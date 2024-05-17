# foes-backpop

to compile the COSMIC shared object library, run the following gfortran commands in the COSMIC repo:

cd COSMIC/cosmic/src

gfortran -shared -fPIC -c assign_remnant.f hrdiag_remnant.f evolv_wrapper.f evolv2.f bpp_array.f checkstate.f comenv.f comprad.f concatkstars.f corerd.f deltat.f dgcore.f gntage.f hrdiag.f instar.f int64.f kick.f mix.f mlwind.f mrenv.f ran3.f rl.f star.f tausworth.f zcnsts.f zfuncs.f

gfortran -shared -o evolve_backpop.so *.o  

cp evolvebin.so /path/to/foes-backpop
