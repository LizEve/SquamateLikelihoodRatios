Begin MrBayes;
    log start filename=tacocat.toxsa.log replace;
    execute tacocat.nex;
    Lset  nst=6  rates=gamma;
    Prset statefreqpr=dirichlet(1,1,1,1);
    mcmcp ngen=1000000 nruns=4 nchains=4 filename=tacocat.toxsa;
    constraint Outgroup hard = homo_sapiens chrysemys_picta alligator_mississippiensis gallus_gallus sphenodon_sp;
    constraint Dibamidae hard = anelytropsis_sp dibamus_sp;
    constraint Gekkota hard = strophurus_sp gonatodes_sp saltorius_sp gekko_sp coleonyx_variegatus lialis_sp;
    constraint Scincoidea hard = cordylosaurus_sp cordylus_sp cricosaura_sp lepidophyma_sp plestiodon_fasciatus tiliqua_sp acontias_sp;
    constraint TeiidaeGymno hard = aspidoscelis_tigris tupinambus_sp pholidobollus_sp;
    constraint Amphi hard = bipes_sp rhineura_sp;
    constraint Serpentes hard = typhlops_jamaicensis micrurus_fulvius python_molurus_2;
    constraint Iguania hard = anolis_carolinensis_2 hydrosaurus_sp uta_stansburiana;
    constraint Anguimorpha hard = lanthanotus_borneensis anniella_pulchra heloderma_suspectum varanus_exanthematicus xenosaurus_platyceps;
    constraint tox hard = typhlops_jamaicensis micrurus_fulvius python_molurus_2 anolis_carolinensis_2 hydrosaurus_sp uta_stansburiana lanthanotus_borneensis anniella_pulchra heloderma_suspectum varanus_exanthematicus xenosaurus_platyceps;
    constraint toxai hard = anolis_carolinensis_2 hydrosaurus_sp uta_stansburiana lanthanotus_borneensis anniella_pulchra heloderma_suspectum varanus_exanthematicus xenosaurus_platyceps;
    constraint toxsi hard = anolis_carolinensis_2 hydrosaurus_sp uta_stansburiana typhlops_jamaicensis micrurus_fulvius python_molurus_2;
    constraint toxsa hard = typhlops_jamaicensis micrurus_fulvius python_molurus_2 lanthanotus_borneensis anniella_pulchra heloderma_suspectum varanus_exanthematicus xenosaurus_platyceps;
    constraint scler hard = anelytropsis_sp dibamus_sp strophurus_sp gonatodes_sp saltorius_sp gekko_sp coleonyx_variegatus lialis_sp plestiodon_fasciatus tiliqua_sp acontias_sp cricosaura_sp lepidophyma_sp cordylosaurus_sp cordylus_sp aspidoscelis_tigris tupinambus_sp pholidobollus_sp bipes_sp lacerta_sp rhineura_sp typhlops_jamaicensis micrurus_fulvius python_molurus_2 lanthanotus_borneensis anniella_pulchra heloderma_suspectum varanus_exanthematicus xenosaurus_platyceps;
    prset topologypr = constraints (toxsa, tox, Outgroup, Dibamidae, Gekkota, Scincoidea, TeiidaeGymno, Amphi, Serpentes, Iguania, Anguimorpha);
    ss nsteps=50;
    log stop;

    log start filename=tacocat.toxsi.log replace;
    mcmcp filename=tacocat.toxsi;
    prset topologypr = constraints (toxsi, tox, Outgroup, Dibamidae, Gekkota, Scincoidea, TeiidaeGymno, Amphi, Serpentes, Iguania, Anguimorpha);
    ss nsteps=50;
    log stop;

    log start filename=tacocat.toxai.log replace;
    mcmcp filename=tacocat.toxai;
    prset topologypr = constraints (toxai, tox, Outgroup, Dibamidae, Gekkota, Scincoidea, TeiidaeGymno, Amphi, Serpentes, Iguania, Anguimorpha);
    ss nsteps=50;
    log stop;

    log start filename=tacocat.toxsc.log replace;
    mcmcp filename=tacocat.toxsc;
    prset topologypr = constraints (scler, Outgroup, Dibamidae, Gekkota, Scincoidea, TeiidaeGymno, Amphi, Serpentes, Iguania, Anguimorpha);
    ss nsteps=50;
    log stop;

End;
