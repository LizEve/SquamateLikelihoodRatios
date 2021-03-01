Begin MrBayes;
    log start filename=tacocat.toxai.log replace;
    execute tacocat.nex;
    constraint Outgroup hard = Alligator_mississippiensis Chrysemys_picta Gopherus_evgoodei Sphenodon_punctatus1 Gallus_gallus;
    constraint Iguania hard = Anolis_carolinensis Pogona_vitticeps Anolis_brasiliensis Anolis_meridionalis Tropidurus_oreadicus;
    constraint Viperidae hard = Crotalus_horridus Deinagkistrodon_actutus Protobothrops_mucrosquamatus Bothrops_lutzi Bothrops_moojeni Bothrops_pauloensis;
    constraint Anguimorpha hard = Dopasia_gracilis Shinisaurus_crocodilurus;
    constraint Gekkota hard = Eublepharis_macularius Gekko_japonicus Paroedura_picta Gymnodactylus_amarali Hemidactylus_mabouia;
    constraint Elapoidea hard = Hydrophis_melanocephalus Laticauda_colubrina Ophiophagus_hannah Micrurus_brasiliensis;
    constraint Lacertidae hard = Lacerta_bilineata Lacerta_viridis Podarcis_muralis Zootoca_vivipara;
    constraint TeiidaeGymno hard = Salvator_meriana Ameiva_ameiva Ameivula_mumbuca Colobosaura_modesta Kentropyx_calcarata Micrablepharus_maximiliani;
    constraint Colubroidea hard = Apostolepis_cearensis Apostolepis_polylepis Chironius_exoletus Erythrolamprus_almadensis Erythrolamprus_poecilogyrus Erythrolamprus_reginae Imantodes_cenchoa Leptodeira_annulata Lygophis_paucidens Oxyrhopus_petolarius Oxyrhopus_trigeminus Philodryas_nattereri Philodryas_olfersii Phimophis_guerini Pseudoboa_neuwiedii Pseudoboa_nigra Psomophis_joberti Sibynomorphus_mikanii Taeniophallus_occipitalis Tantilla_melanocephala Thamnodynastes_hypoconia Xenodon_merremi;
    constraint Scincoidea hard = Brasiliscincus_heathi Copeoglossum_nigropunctatum Notomabuya_frenata;
    constraint Snake hard = Crotalus_horridus Deinagkistrodon_actutus Hydrophis_melanocephalus Laticauda_colubrina Ophiophagus_hannah Protobothrops_mucrosquamatus Python_bivittatus Apostolepis_cearensis Apostolepis_polylepis Bothrops_lutzi Bothrops_moojeni Bothrops_pauloensis Chironius_exoletus Corallus_hortulanus Erythrolamprus_almadensis Erythrolamprus_poecilogyrus Erythrolamprus_reginae Imantodes_cenchoa Leptodeira_annulata Liotyphlops_ternetzii Lygophis_paucidens Micrurus_brasiliensis Oxyrhopus_petolarius Oxyrhopus_trigeminus Philodryas_nattereri Philodryas_olfersii Phimophis_guerini Pseudoboa_neuwiedii Pseudoboa_nigra Psomophis_joberti Sibynomorphus_mikanii Taeniophallus_occipitalis Tantilla_melanocephala Thamnodynastes_hypoconia Trilepida_brasiliensis Typhlops_brongersmianus Xenodon_merremi;
    constraint Sclero hard = Crotalus_horridus Deinagkistrodon_actutus Dopasia_gracilis Eublepharis_macularius Gekko_japonicus Hydrophis_melanocephalus Lacerta_bilineata Lacerta_viridis Laticauda_colubrina Ophiophagus_hannah Paroedura_picta Podarcis_muralis Protobothrops_mucrosquamatus Python_bivittatus Salvator_meriana Shinisaurus_crocodilurus Zootoca_vivipara Ameiva_ameiva Ameivula_mumbuca Amphisbaena_alba Apostolepis_cearensis Apostolepis_polylepis Bothrops_lutzi Bothrops_moojeni Bothrops_pauloensis Brasiliscincus_heathi Chironius_exoletus Colobosaura_modesta Copeoglossum_nigropunctatum Corallus_hortulanus Erythrolamprus_almadensis Erythrolamprus_poecilogyrus Erythrolamprus_reginae Gymnodactylus_amarali Hemidactylus_mabouia Imantodes_cenchoa Kentropyx_calcarata Leptodeira_annulata Liotyphlops_ternetzii Lygophis_paucidens Micrablepharus_maximiliani Micrurus_brasiliensis Notomabuya_frenata Oxyrhopus_petolarius Oxyrhopus_trigeminus Philodryas_nattereri Philodryas_olfersii Phimophis_guerini Pseudoboa_neuwiedii Pseudoboa_nigra Psomophis_joberti Sibynomorphus_mikanii Taeniophallus_occipitalis Tantilla_melanocephala Thamnodynastes_hypoconia Trilepida_brasiliensis Typhlops_brongersmianus Xenodon_merremi;
    constraint ToxPoly hard = Anolis_carolinensis Crotalus_horridus Deinagkistrodon_actutus Dopasia_gracilis Hydrophis_melanocephalus Laticauda_colubrina Ophiophagus_hannah Pogona_vitticeps Protobothrops_mucrosquamatus Python_bivittatus Shinisaurus_crocodilurus Anolis_brasiliensis Anolis_meridionalis Apostolepis_cearensis Apostolepis_polylepis Bothrops_lutzi Bothrops_moojeni Bothrops_pauloensis Chironius_exoletus Corallus_hortulanus Erythrolamprus_almadensis Erythrolamprus_poecilogyrus Erythrolamprus_reginae Imantodes_cenchoa Leptodeira_annulata Liotyphlops_ternetzii Lygophis_paucidens Micrurus_brasiliensis Oxyrhopus_petolarius Oxyrhopus_trigeminus Philodryas_nattereri Philodryas_olfersii Phimophis_guerini Pseudoboa_neuwiedii Pseudoboa_nigra Psomophis_joberti Sibynomorphus_mikanii Taeniophallus_occipitalis Tantilla_melanocephala Thamnodynastes_hypoconia Trilepida_brasiliensis Tropidurus_oreadicus Typhlops_brongersmianus Xenodon_merremi;
    constraint ToxAI hard = Anolis_carolinensis Dopasia_gracilis Pogona_vitticeps Shinisaurus_crocodilurus Anolis_brasiliensis Anolis_meridionalis Tropidurus_oreadicus;
    constraint ToxSA hard = Crotalus_horridus Deinagkistrodon_actutus Dopasia_gracilis Hydrophis_melanocephalus Laticauda_colubrina Ophiophagus_hannah Protobothrops_mucrosquamatus Python_bivittatus Shinisaurus_crocodilurus Apostolepis_cearensis Apostolepis_polylepis Bothrops_lutzi Bothrops_moojeni Bothrops_pauloensis Chironius_exoletus Corallus_hortulanus Erythrolamprus_almadensis Erythrolamprus_poecilogyrus Erythrolamprus_reginae Imantodes_cenchoa Leptodeira_annulata Liotyphlops_ternetzii Lygophis_paucidens Micrurus_brasiliensis Oxyrhopus_petolarius Oxyrhopus_trigeminus Philodryas_nattereri Philodryas_olfersii Phimophis_guerini Pseudoboa_neuwiedii Pseudoboa_nigra Psomophis_joberti Sibynomorphus_mikanii Taeniophallus_occipitalis Tantilla_melanocephala Thamnodynastes_hypoconia Trilepida_brasiliensis Typhlops_brongersmianus Xenodon_merremi;
    constraint ToxSI hard = Anolis_carolinensis Crotalus_horridus Deinagkistrodon_actutus Hydrophis_melanocephalus Laticauda_colubrina Ophiophagus_hannah Pogona_vitticeps Protobothrops_mucrosquamatus Python_bivittatus Anolis_brasiliensis Anolis_meridionalis Apostolepis_cearensis Apostolepis_polylepis Bothrops_lutzi Bothrops_moojeni Bothrops_pauloensis Chironius_exoletus Corallus_hortulanus Erythrolamprus_almadensis Erythrolamprus_poecilogyrus Erythrolamprus_reginae Imantodes_cenchoa Leptodeira_annulata Liotyphlops_ternetzii Lygophis_paucidens Micrurus_brasiliensis Oxyrhopus_petolarius Oxyrhopus_trigeminus Philodryas_nattereri Philodryas_olfersii Phimophis_guerini Pseudoboa_neuwiedii Pseudoboa_nigra Psomophis_joberti Sibynomorphus_mikanii Taeniophallus_occipitalis Tantilla_melanocephala Thamnodynastes_hypoconia Trilepida_brasiliensis Tropidurus_oreadicus Typhlops_brongersmianus Xenodon_merremi;
    prset topologypr = constraints (Outgroup,Iguania,Viperidae,Anguimorpha,Gekkota,Elapoidea,Lacertidae,TeiidaeGymno,Colubroidea,Scincoidea,Snake,ToxPoly,ToxAI);
    prset statefreqpr=dirichlet(1,1,1,1);
    mcmcp ngen=2000000 nruns=2 nchains=4 filename=tacocat.toxai;
    ss nsteps=50 burninSS=-4;
    log stop;
End;
