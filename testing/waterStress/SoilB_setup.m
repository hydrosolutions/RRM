prm.modelname = 'SoilB';
prm.seed = 2;
prm.propagate = @StepSoilB;
prm.observe = @SoilB_observe;
prm.n = 1;
prm.p = 1;
prm.collapseratio = 0.005;
prm.m=10000;
prm.method = 'EnKF';
prm.inflation = 1.35;
prm.cycle_nstep = 5;
prm.obs_variance = 2;
prm.nspinup = 5;
prm.ncycle = 25;
prm.nitermax=1;

