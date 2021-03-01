Begin MrBayes;
    log start filename=tacocat.sclero.log replace;
    execute tacocat.nex;
    constraint TeiidaeGymno hard = Ameiva_ameiva Ameivula_mumbuca Colobosaura_modesta Kentropyx_calcarata Micrablepharus_maximiliani;
    constraint Iguania hard = Anolis_brasiliensis Anolis_meridionalis Tropidurus_oreadicus;
    constraint Colubroidea hard = Apostolepis_cearensis Apostolepis_polylepis Chironius_exoletus Erythrolamprus_almadensis Erythrolamprus_poecilogyrus Erythrolamprus_reginae Imantodes_cenchoa Leptodeira_annulata Lygophis_paucidens Oxyrhopus_petolarius Oxyrhopus_trigeminus Philodryas_nattereri Philodryas_olfersii Phimophis_guerini Pseudoboa_neuwiedii Pseudoboa_nigra Psomophis_joberti Sibynomorphus_mikanii Taeniophallus_occipitalis Tantilla_melanocephala Thamnodynastes_hypoconia Xenodon_merremi;
    constraint Viperidae hard = Bothrops_lutzi Bothrops_moojeni Bothrops_pauloensis;
    constraint Scincoidea hard = Brasiliscincus_heathi Copeoglossum_nigropunctatum Notomabuya_frenata;
    constraint Gekkota hard = Gymnodactylus_amarali Hemidactylus_mabouia;
    constraint Outgroup hard = Gallus_gallus;
    constraint Snake hard = Apostolepis_cearensis Apostolepis_polylepis Bothrops_lutzi Bothrops_moojeni Bothrops_pauloensis Chironius_exoletus Corallus_hortulanus Erythrolamprus_almadensis Erythrolamprus_poecilogyrus Erythrolamprus_reginae Imantodes_cenchoa Leptodeira_annulata Liotyphlops_ternetzii Lygophis_paucidens Micrurus_brasiliensis Oxyrhopus_petolarius Oxyrhopus_trigeminus Philodryas_nattereri Philodryas_olfersii Phimophis_guerini Pseudoboa_neuwiedii Pseudoboa_nigra Psomophis_joberti Sibynomorphus_mikanii Taeniophallus_occipitalis Tantilla_melanocephala Thamnodynastes_hypoconia Trilepida_brasiliensis Typhlops_brongersmianus Xenodon_merremi;
    constraint Sclero hard = Ameiva_ameiva Ameivula_mumbuca Amphisbaena_alba Apostolepis_cearensis Apostolepis_polylepis Bothrops_lutzi Bothrops_moojeni Bothrops_pauloensis Brasiliscincus_heathi Chironius_exoletus Colobosaura_modesta Copeoglossum_nigropunctatum Corallus_hortulanus Erythrolamprus_almadensis Erythrolamprus_poecilogyrus Erythrolamprus_reginae Gymnodactylus_amarali Hemidactylus_mabouia Imantodes_cenchoa Kentropyx_calcarata Leptodeira_annulata Liotyphlops_ternetzii Lygophis_paucidens Micrablepharus_maximiliani Micrurus_brasiliensis Notomabuya_frenata Oxyrhopus_petolarius Oxyrhopus_trigeminus Philodryas_nattereri Philodryas_olfersii Phimophis_guerini Pseudoboa_neuwiedii Pseudoboa_nigra Psomophis_joberti Sibynomorphus_mikanii Taeniophallus_occipitalis Tantilla_melanocephala Thamnodynastes_hypoconia Trilepida_brasiliensis Typhlops_brongersmianus Xenodon_merremi;
    constraint ToxPoly hard = Anolis_brasiliensis Anolis_meridionalis Apostolepis_cearensis Apostolepis_polylepis Bothrops_lutzi Bothrops_moojeni Bothrops_pauloensis Chironius_exoletus Corallus_hortulanus Erythrolamprus_almadensis Erythrolamprus_poecilogyrus Erythrolamprus_reginae Imantodes_cenchoa Leptodeira_annulata Liotyphlops_ternetzii Lygophis_paucidens Micrurus_brasiliensis Oxyrhopus_petolarius Oxyrhopus_trigeminus Philodryas_nattereri Philodryas_olfersii Phimophis_guerini Pseudoboa_neuwiedii Pseudoboa_nigra Psomophis_joberti Sibynomorphus_mikanii Taeniophallus_occipitalis Tantilla_melanocephala Thamnodynastes_hypoconia Trilepida_brasiliensis Tropidurus_oreadicus Typhlops_brongersmianus Xenodon_merremi;
    prset topologypr = constraints (Sclero,TeiidaeGymno,Iguania,Colubroidea,Viperidae,Scincoidea,Gekkota,Snake);
    prset statefreqpr=dirichlet(1,1,1,1);
    mcmcp ngen=1000000 nruns=4 nchains=4 filename=tacocat.sclero;
    ss nsteps=50 burninSS=-4;
    log stop;
End;
