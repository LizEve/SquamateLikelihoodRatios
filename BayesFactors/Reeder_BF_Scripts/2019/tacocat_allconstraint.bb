Begin MrBayes;
    log start filename=tacocat.XXX.log replace;
    execute tacocat.nex;
    constraint Amphisbaenia hard = Amphisbaena_fuliginosa Bipes_biporus Bipes_canaliculatus Diplometopon_zarudnyi Geocalamus_acutus Rhineura_floridana Trogonophis_wiegmanni;
	constraint Anguimorpha hard = Anniella_pulchra Celestus_enneagrammus Elgaria_multicarinata Heloderma_horridum Heloderma_suspectum Lanthanotus_borneensis Shinisaurus_crocodilurus Varanus_acanthurus Varanus_exanthematicus Varanus_salvator Xenosaurus_grandis Xenosaurus_platyceps Pseudopus_apodus;
	constraint Dibamidae hard = Anelytropsis_papillosus Dibamus_novaeguineae;
	constraint Gekkota hard = Aeluroscalobates_felinus Coleonyx_variegatus Delma_borea Eublepharis_macularius Gekko_gecko Gonatodes_albogularis Lialis_burtonis Phelsuma_lineata Rhacodactylus_auriculatus Saltuarius_cornutus Strophurus_ciliaris Teratoscincus;
	constraint TeiidaeGymno hard = Colobosaura_modesta Pholidobolus Aspidoscelis_tigris Callopistes_maculatus Teius_teyou Tupinambis_teguixin;
	constraint Iguania hard = Agama_agama Anolis_carolinensis Basiliscus_basiliscus Brachylophus_fasciatus Brookesia_brygooi Calotes_emma Chalarodon_madagascariensis Chamaeleo Corytophanes_cristatus Crotaphytus_collaris Dipsosaurus_dorsalis Enyalioides_laticeps Gambelia_wislizenii Hydrosaurus Leiocephalus_barahonensis Leiolepis_belliana Leiosaurus_catamarcensis Liolaemus_bellii Morunasaurus_annularis Oplurus_cyclurus Petrosaurus_mearnsi Phrynosoma_platyrhinos Phymaturus_palluma Physignathus_cocincinus Pogona_vitticeps Polychrus_marmoratus Pristidactylus_torquatus Sauromalus_ater Sceloporus_variabilis Stenocercus_guentheri Uma_scoparia Uranoscodon_superciliosus Uromastyx_aegyptus Urostrophus_vautieri Uta_stansburiana Plica_plica;
	constraint Lacertidae hard = Lacerta_viridis Takydromus_ocellatus;
	constraint Scincoidea hard = Acontias Amphiglossus_splendidus Brachymeles_gracilis Cordylosaurus_subtesselatus Cordylus_mossambicus Cricosaura_typica Eugongylus_rufescens Feylinia_polylepis Lepidophyma_flavimaculatu Platysaurus Plestiodon_fasciatus Scincus Sphenomorphus_solomonis Tiliqua_scincoides Trachylepis_quinquetaeniata Xantusia_vigilis Zonosaurus_ornatus;
	constraint Outgroup hard = Alligator Chelydra Crocodylus Dromaius Gallus Mus Podocnemis Sphenodon_punctatus Tachyglos Homo;
	constraint Boinae hard = Boa_constrictor Calabaria_reinhardtii Epicrates_striatus Eryx_colubrinus Exiliboa_placata Lichanura_trivirgata Ungaliophis_continentalis;
	constraint Colubroidea hard = Afronatrix_anoscopus Amphiesma_stolata Coluber_constrictor Diadophis_punctatus Homalopsis_buccata Lampropeltis_getula Natrix_natrix Pareas_hamptoni Thamnophis_marcianus Xenochrophis_piscator Xenodermus_javanicus;
	constraint Elapoidea hard = Aparallactus_werneri Atractaspis_irregularis Laticauda_colubrina Lycophidion_capense Micrurus_fulvius Naja Notechis_scutatus;
	constraint Pythonidae hard = Aspidites_melanocephalus Python_molurus;
	constraint Tropidophiidae hard = Trachyboa_boulengeri Tropidophis_haetianus;
	constraint Viperidae hard = Agkistrodon_contortrix Azemiops_feae Bothrops_asper Causus Daboia_russelli Lachesis_muta;
	constraint Snake hard = Acrochordus_granulatus Anilius_scytale Liotyphlops_albirostris Boa_constrictor Calabaria_reinhardtii Epicrates_striatus Eryx_colubrinus Exiliboa_placata Lichanura_trivirgata Ungaliophis_continentalis Casarea_dussumieri Afronatrix_anoscopus Amphiesma_stolata Coluber_constrictor Diadophis_punctatus Homalopsis_buccata Lampropeltis_getula Natrix_natrix Pareas_hamptoni Thamnophis_marcianus Xenochrophis_piscator Xenodermus_javanicus Cylindrophis_rufus Aparallactus_werneri Atractaspis_irregularis Laticauda_colubrina Lycophidion_capense Micrurus_fulvius Naja Notechis_scutatus Leptotyphlops Loxocemus_bicolor Aspidites_melanocephalus Python_molurus Trachyboa_boulengeri Tropidophis_haetianus Typhlops_jamaicensis Uropeltis_melanogaster Agkistrodon_contortrix Azemiops_feae Bothrops_asper Causus Daboia_russelli Lachesis_muta Xenopeltis_unicolor;
	constraint Sclero hard = Amphisbaena_fuliginosa Bipes_biporus Bipes_canaliculatus Diplometopon_zarudnyi Geocalamus_acutus Rhineura_floridana Trogonophis_wiegmanni Anniella_pulchra Celestus_enneagrammus Elgaria_multicarinata Heloderma_horridum Heloderma_suspectum Lanthanotus_borneensis Shinisaurus_crocodilurus Varanus_acanthurus Varanus_exanthematicus Varanus_salvator Xenosaurus_grandis Xenosaurus_platyceps Pseudopus_apodus Anelytropsis_papillosus Dibamus_novaeguineae Aeluroscalobates_felinus Coleonyx_variegatus Delma_borea Eublepharis_macularius Gekko_gecko Gonatodes_albogularis Lialis_burtonis Phelsuma_lineata Rhacodactylus_auriculatus Saltuarius_cornutus Strophurus_ciliaris Teratoscincus Colobosaura_modesta Pholidobolus Lacerta_viridis Takydromus_ocellatus Acontias Amphiglossus_splendidus Brachymeles_gracilis Cordylosaurus_subtesselatus Cordylus_mossambicus Cricosaura_typica Eugongylus_rufescens Feylinia_polylepis Lepidophyma_flavimaculatu Platysaurus Plestiodon_fasciatus Scincus Sphenomorphus_solomonis Tiliqua_scincoides Trachylepis_quinquetaeniata Xantusia_vigilis Zonosaurus_ornatus Aspidoscelis_tigris Callopistes_maculatus Teius_teyou Tupinambis_teguixin Acrochordus_granulatus Anilius_scytale Liotyphlops_albirostris Boa_constrictor Calabaria_reinhardtii Epicrates_striatus Eryx_colubrinus Exiliboa_placata Lichanura_trivirgata Ungaliophis_continentalis Casarea_dussumieri Afronatrix_anoscopus Amphiesma_stolata Coluber_constrictor Diadophis_punctatus Homalopsis_buccata Lampropeltis_getula Natrix_natrix Pareas_hamptoni Thamnophis_marcianus Xenochrophis_piscator Xenodermus_javanicus Cylindrophis_rufus Aparallactus_werneri Atractaspis_irregularis Laticauda_colubrina Lycophidion_capense Micrurus_fulvius Naja Notechis_scutatus Leptotyphlops Loxocemus_bicolor Aspidites_melanocephalus Python_molurus Trachyboa_boulengeri Tropidophis_haetianus Typhlops_jamaicensis Uropeltis_melanogaster Agkistrodon_contortrix Azemiops_feae Bothrops_asper Causus Daboia_russelli Lachesis_muta Xenopeltis_unicolor;
	constraint ToxPoly hard = Anniella_pulchra Celestus_enneagrammus Elgaria_multicarinata Heloderma_horridum Heloderma_suspectum Lanthanotus_borneensis Shinisaurus_crocodilurus Varanus_acanthurus Varanus_exanthematicus Varanus_salvator Xenosaurus_grandis Xenosaurus_platyceps Pseudopus_apodus Agama_agama Anolis_carolinensis Basiliscus_basiliscus Brachylophus_fasciatus Brookesia_brygooi Calotes_emma Chalarodon_madagascariensis Chamaeleo Corytophanes_cristatus Crotaphytus_collaris Dipsosaurus_dorsalis Enyalioides_laticeps Gambelia_wislizenii Hydrosaurus Leiocephalus_barahonensis Leiolepis_belliana Leiosaurus_catamarcensis Liolaemus_bellii Morunasaurus_annularis Oplurus_cyclurus Petrosaurus_mearnsi Phrynosoma_platyrhinos Phymaturus_palluma Physignathus_cocincinus Pogona_vitticeps Polychrus_marmoratus Pristidactylus_torquatus Sauromalus_ater Sceloporus_variabilis Stenocercus_guentheri Uma_scoparia Uranoscodon_superciliosus Uromastyx_aegyptus Urostrophus_vautieri Uta_stansburiana Plica_plica Acrochordus_granulatus Anilius_scytale Liotyphlops_albirostris Boa_constrictor Calabaria_reinhardtii Epicrates_striatus Eryx_colubrinus Exiliboa_placata Lichanura_trivirgata Ungaliophis_continentalis Casarea_dussumieri Afronatrix_anoscopus Amphiesma_stolata Coluber_constrictor Diadophis_punctatus Homalopsis_buccata Lampropeltis_getula Natrix_natrix Pareas_hamptoni Thamnophis_marcianus Xenochrophis_piscator Xenodermus_javanicus Cylindrophis_rufus Aparallactus_werneri Atractaspis_irregularis Laticauda_colubrina Lycophidion_capense Micrurus_fulvius Naja Notechis_scutatus Leptotyphlops Loxocemus_bicolor Aspidites_melanocephalus Python_molurus Trachyboa_boulengeri Tropidophis_haetianus Typhlops_jamaicensis Uropeltis_melanogaster Agkistrodon_contortrix Azemiops_feae Bothrops_asper Causus Daboia_russelli Lachesis_muta Xenopeltis_unicolor;
	constraint ToxAI hard = Anniella_pulchra Celestus_enneagrammus Elgaria_multicarinata Heloderma_horridum Heloderma_suspectum Lanthanotus_borneensis Shinisaurus_crocodilurus Varanus_acanthurus Varanus_exanthematicus Varanus_salvator Xenosaurus_grandis Xenosaurus_platyceps Pseudopus_apodus Agama_agama Anolis_carolinensis Basiliscus_basiliscus Brachylophus_fasciatus Brookesia_brygooi Calotes_emma Chalarodon_madagascariensis Chamaeleo Corytophanes_cristatus Crotaphytus_collaris Dipsosaurus_dorsalis Enyalioides_laticeps Gambelia_wislizenii Hydrosaurus Leiocephalus_barahonensis Leiolepis_belliana Leiosaurus_catamarcensis Liolaemus_bellii Morunasaurus_annularis Oplurus_cyclurus Petrosaurus_mearnsi Phrynosoma_platyrhinos Phymaturus_palluma Physignathus_cocincinus Pogona_vitticeps Polychrus_marmoratus Pristidactylus_torquatus Sauromalus_ater Sceloporus_variabilis Stenocercus_guentheri Uma_scoparia Uranoscodon_superciliosus Uromastyx_aegyptus Urostrophus_vautieri Uta_stansburiana Plica_plica;
	constraint ToxSA hard = Anniella_pulchra Celestus_enneagrammus Elgaria_multicarinata Heloderma_horridum Heloderma_suspectum Lanthanotus_borneensis Shinisaurus_crocodilurus Varanus_acanthurus Varanus_exanthematicus Varanus_salvator Xenosaurus_grandis Xenosaurus_platyceps Pseudopus_apodus Acrochordus_granulatus Anilius_scytale Liotyphlops_albirostris Boa_constrictor Calabaria_reinhardtii Epicrates_striatus Eryx_colubrinus Exiliboa_placata Lichanura_trivirgata Ungaliophis_continentalis Casarea_dussumieri Afronatrix_anoscopus Amphiesma_stolata Coluber_constrictor Diadophis_punctatus Homalopsis_buccata Lampropeltis_getula Natrix_natrix Pareas_hamptoni Thamnophis_marcianus Xenochrophis_piscator Xenodermus_javanicus Cylindrophis_rufus Aparallactus_werneri Atractaspis_irregularis Laticauda_colubrina Lycophidion_capense Micrurus_fulvius Naja Notechis_scutatus Leptotyphlops Loxocemus_bicolor Aspidites_melanocephalus Python_molurus Trachyboa_boulengeri Tropidophis_haetianus Typhlops_jamaicensis Uropeltis_melanogaster Agkistrodon_contortrix Azemiops_feae Bothrops_asper Causus Daboia_russelli Lachesis_muta Xenopeltis_unicolor;
	constraint ToxSI hard = Agama_agama Anolis_carolinensis Basiliscus_basiliscus Brachylophus_fasciatus Brookesia_brygooi Calotes_emma Chalarodon_madagascariensis Chamaeleo Corytophanes_cristatus Crotaphytus_collaris Dipsosaurus_dorsalis Enyalioides_laticeps Gambelia_wislizenii Hydrosaurus Leiocephalus_barahonensis Leiolepis_belliana Leiosaurus_catamarcensis Liolaemus_bellii Morunasaurus_annularis Oplurus_cyclurus Petrosaurus_mearnsi Phrynosoma_platyrhinos Phymaturus_palluma Physignathus_cocincinus Pogona_vitticeps Polychrus_marmoratus Pristidactylus_torquatus Sauromalus_ater Sceloporus_variabilis Stenocercus_guentheri Uma_scoparia Uranoscodon_superciliosus Uromastyx_aegyptus Urostrophus_vautieri Uta_stansburiana Plica_plica Acrochordus_granulatus Anilius_scytale Liotyphlops_albirostris Boa_constrictor Calabaria_reinhardtii Epicrates_striatus Eryx_colubrinus Exiliboa_placata Lichanura_trivirgata Ungaliophis_continentalis Casarea_dussumieri Afronatrix_anoscopus Amphiesma_stolata Coluber_constrictor Diadophis_punctatus Homalopsis_buccata Lampropeltis_getula Natrix_natrix Pareas_hamptoni Thamnophis_marcianus Xenochrophis_piscator Xenodermus_javanicus Cylindrophis_rufus Aparallactus_werneri Atractaspis_irregularis Laticauda_colubrina Lycophidion_capense Micrurus_fulvius Naja Notechis_scutatus Leptotyphlops Loxocemus_bicolor Aspidites_melanocephalus Python_molurus Trachyboa_boulengeri Tropidophis_haetianus Typhlops_jamaicensis Uropeltis_melanogaster Agkistrodon_contortrix Azemiops_feae Bothrops_asper Causus Daboia_russelli Lachesis_muta Xenopeltis_unicolor;
	prset topologypr = constraints (Sclero,ToxPoly,ToxAI,ToxSA,ToxSI,Amphisbaenia,Anguimorpha,Dibamidae,Gekkota,TeiidaeGymno,Iguania,Lacertidae,Scincoidea,Outgroup,Boinae,Colubroidea,Elapoidea,Pythonidae,Tropidophiidae,Viperidae,Snake);
    prset statefreqpr=dirichlet(1,1,1,1);
    mcmcp ngen=1000000 nruns=4 nchains=4 filename=tacocat.XXX;
    ss nsteps=50 burninSS=-4;
    log stop;
End;