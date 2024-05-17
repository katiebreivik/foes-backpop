import numpy as np
import pandas as pd
import emcee

from ctypes import *
from scipy.stats import gaussian_kde
from multiprocessing import Pool
from cosmic import _evolvebin

import os
from scipy.stats import multivariate_normal

libc = cdll.LoadLibrary("./evolve_backpop.so")
np.set_printoptions(suppress=True)

def evolv2(m1, m2, logtb, e, alpha_2, gam, acc_lim_2, vk1, theta1, phi1, omega1):
    m2, m1 = np.sort([m1,m2],axis=0)
    tb = 10**logtb
    z = byref(c_double(0.02))
    gam = byref(c_double(gam))
    zpars = np.zeros(20).ctypes.data_as(POINTER(c_double))
    natal_kick = np.zeros((2,5))
    natal_kick[0,0] = vk1
    natal_kick[0,1] = phi1
    natal_kick[0,2] = theta1
    natal_kick[0,3] = omega1
    natal_kick[0,4] = 3
    natal_kick = natal_kick.T.flatten().ctypes.data_as(POINTER(c_double))
    alpha = np.zeros(2)
    alpha[0] = 0.0
    alpha[1] = alpha_2
    alpha = alpha.flatten().ctypes.data_as(POINTER(c_double))
    
    acc_lim = np.zeros(2)
    acc_lim[0] = 0.0
    acc_lim[1] = acc_lim_2
    acc_lim = acc_lim.flatten().ctypes.data_as(POINTER(c_double))

    
    libc.evolv2_global_(z,zpars,alpha,natal_kick,acc_lim,gam)

    mass = np.array([m1,m2]).ctypes.data_as(POINTER(c_double))
    mass0 = np.array([m1,m2]).ctypes.data_as(POINTER(c_double))
    epoch = np.array([0.0,0.0]).ctypes.data_as(POINTER(c_double))
    ospin = np.array([0.0,0.0]).ctypes.data_as(POINTER(c_double))
    tb = byref(c_double(tb))
    ecc = byref(c_double(e))
    tphysf = byref(c_double(13700.0))
    dtp = byref(c_double(0.0))
    rad = np.array([0.0,0.0]).ctypes.data_as(POINTER(c_double))
    lumin = np.array([0.0,0.0]).ctypes.data_as(POINTER(c_double))
    massc = np.array([0.0,0.0]).ctypes.data_as(POINTER(c_double))
    radc = np.array([0.0,0.0]).ctypes.data_as(POINTER(c_double))
    menv = np.array([0.0,0.0]).ctypes.data_as(POINTER(c_double))
    renv = np.array([0.0,0.0]).ctypes.data_as(POINTER(c_double))
    B_0 = np.array([0.0,0.0]).ctypes.data_as(POINTER(c_double))
    bacc = np.array([0.0,0.0]).ctypes.data_as(POINTER(c_double))
    tacc = np.array([0.0,0.0]).ctypes.data_as(POINTER(c_double))
    tms = np.array([0.0,0.0]).ctypes.data_as(POINTER(c_double))
    bhspin = np.array([0.0,0.0]).ctypes.data_as(POINTER(c_double))
    tphys = byref(c_double(0.0))
    bkick = np.zeros(20).ctypes.data_as(POINTER(c_double))
    kick_info = np.zeros(34).ctypes.data_as(POINTER(c_double)) # Fortran treat n-D array differently than numpy
    bpp_index_out = byref(c_int64(0))
    bcm_index_out = byref(c_int64(0))
    kick_info_out = np.zeros(34).ctypes.data_as(POINTER(c_double))

    bpp_out=np.zeros([1000,43]).flatten().ctypes.data_as(POINTER(c_double))
    kstar = np.array([1,1]).ctypes.data_as(POINTER(c_double))
    libc.evolv2_(kstar,mass,tb,ecc,z,tphysf,
    dtp,mass0,rad,lumin,massc,radc,
    menv,renv,ospin,B_0,bacc,tacc,epoch,tms,
    bhspin,tphys,zpars,bkick,kick_info,
    bpp_index_out,bcm_index_out,kick_info_out,
    bpp_out)
 
    
    bpp = bpp_out._arr.reshape(43,1000)[:,0:bpp_index_out._obj.value].T
    bpp = pd.DataFrame(bpp, 
                       columns=['tphys', 'mass_1', 'mass_2', 'kstar_1', 'kstar_2', 
                                'sep', 'porb', 'ecc', 'RRLO_1', 'RRLO_2', 'evol_type', 
                                'aj_1', 'aj_2', 'tms_1', 'tms_2', 'massc_1', 'massc_2', 
                                'rad_1', 'rad_2', 'mass0_1', 'mass0_2', 'lum_1', 'lum_2', 
                                'teff_1', 'teff_2', 'radc_1', 'radc_2', 'menv_1', 'menv_2', 
                                'renv_1', 'renv_2', 'omega_spin_1', 'omega_spin_2', 'B_1', 
                                'B_2', 'bacc_1', 'bacc_2', 'tacc_1', 'tacc_2', 'epoch_1', 
                                'epoch_2', 'bhspin_1', 'bhspin_2'])
    
    
    # make your selection on kstar types and orbit properties
    bpp = bpp.loc[(bpp.kstar_1 == 13) & (bpp.kstar_2 == 7) & (bpp.porb > 0.0)]
    if len(bpp) >= 1:
        bpp = bpp.iloc[0]
        dat_out = bpp[['mass_1', 'mass_2', 'porb', 'ecc']].to_numpy()
        return dat_out, kick_info_out._arr.reshape(17,2).T,bpp
    else:
        return np.zeros(4),kick_info_out._arr.reshape(17,2).T,np.array([0.0, 0.0, 0.0, 0.0])


#Porb: 4194.7 +/- 112
#Eccentricity: 0.7262  +/- 0.0056 
#Mstar: 0.76  +/- 0.05
#M_BH: 32.8 +/- 0.82
#mean = np.array([1.1, 1.0, 10])
#cov = np.array([[0.1**2, 0, 0], [0, 0.5**2, 0], [0, 0, 1.0**2]])

Porb_m = 4.17 / 24  # days 
Porb_m_err = 0.1 / 24
Ecc_m = 0.01
Ecc_m_err = 0.001
M_ss = 4 
M_ss_err = 0.5
M_ns = 1.3 
M_ns_err = 0.2 
mean = np.array([M_ns, M_ss, Porb_m, Ecc_m])
cov = np.array([[M_ns_err**2, 0, 0, 0], [0, M_ss_err**2, 0, 0], [0, 0, Porb_m_err**2, 0], [0, 0, 0, Ecc_m_err**2]])
rv = multivariate_normal(mean, cov)
#m1, m2, tb, e, alpha, vk1, theta1, phi1, omega1, vk2, theta2, phi2, omega2, acc_lim, qc_kstar2, qc_kstar3, Z
m1lo = 1.0
m2lo = 0.71
tblo = 10.0
elo = 0.0
alphalo_2 = 0.1
vklo = 0.0
thetalo = 0.0
philo = -90.0
omegalo = 0.0
acc_lim_2_lo = 0.0
gamma_lo = 0.0
Zlo = 0.00214

m1hi = 10.0
m2hi = 0.81
tbhi = 15000.0
ehi = 0.9
alphahi_2 = 20.0
vkhi = 300.0
thetahi = 360.0
phihi = 90.0
omegahi = 360
acc_lim_2_hi = 1.0
gamma_hi = 1.0
Zhi = 0.0036

#m1, m2, logtb, e, alpha_1, vk1, theta1, phi1, omega1, logZ
lower_bound = np.array([m1lo, m2lo, np.log10(tblo), elo, alphalo_2, gamma_lo, acc_lim_2_lo, vklo, thetalo, philo, omegalo])
upper_bound = np.array([m1hi, m2hi, np.log10(tbhi), ehi, alphahi_2, gamma_hi, acc_lim_2_hi, vkhi, thetahi, phihi, omegahi])
param_names = ["m1", "m2", "logtb", "e", "alpha1_1", "vk1", "theta1", "phi1", "omega1"]

def likelihood(coord):
    for i in range(len(coord)):
        if (coord[i]<lower_bound[i] or coord[i]>upper_bound[i]):
            return -np.inf
    result = evolv2(*coord)
    fit_coord = result[0]
    if np.any(fit_coord) == 0.0: return -np.inf
    return rv.logpdf(fit_coord) 

n_dim = len(lower_bound)
n_walkers = 1024
p0 = np.random.uniform(lower_bound, upper_bound, size=(n_walkers, len(lower_bound)))

n_steps = 100

params_out = ["m1", "m2", "porb"]
params_in = param_names
sampler = emcee.EnsembleSampler(n_walkers, n_dim, likelihood)#, pool=pool)
sampler.run_mcmc(p0, n_steps, progress=True)

#for p in p0:
#print(p)
#_ = likelihood(p)


np.savez('./test_bl',nwalkers=n_walkers,n_steps=n_steps,chain=sampler.chain)
