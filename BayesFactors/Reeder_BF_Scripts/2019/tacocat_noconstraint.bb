Begin MrBayes;
    log start filename=tacocat.toxsi.log replace;
    execute tacocat.nex;

    prset statefreqpr=dirichlet(1,1,1,1);
    mcmcp ngen=2000000 nruns=2 nchains=4 filename=tacocat.toxsi;
    ss nsteps=50 burninSS=-4;
    log stop;
End;
