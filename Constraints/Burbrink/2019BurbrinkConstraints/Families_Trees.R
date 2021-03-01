library(ape)

x <- "((chrysemys_picta,alligator_mississippiensis,gallus_gallus),((strophurus_sp,gonatodes_sp,saltorius_sp,lialis_sp),(cordylosaurus_sp,cordylus_sp,cricosaura_sp,lepidophyma_sp,plestiodon_fasciatus,tiliqua_sp),(aspidoscelis_tigris,tupinambus_sp,pholidobollus_sp),(bipes_sp,rhineura_sp),micrurus_fulvius,(anniella_pulchra,heloderma_suspectum,xenosaurus_platyceps)));"
y <- "((homo_sapiens,chrysemys_picta,alligator_mississippiensis,gallus_gallus,sphenodon_sp),(anolis_carolinensis_2,hydrosaurus_sp),(dibamus_sp,(strophurus_sp,gonatodes_sp,saltorius_sp,gekko_sp,coleonyx_variegatus,lialis_sp),(cordylus_sp,plestiodon_fasciatus,tiliqua_sp),(aspidoscelis_tigris,tupinambus_sp),rhineura_sp,(typhlops_jamaicensis,micrurus_fulvius,python_molurus_2),(lanthanotus_borneensis,anniella_pulchra,xenosaurus_platyceps)));"
z <- "(Kentropyx_calcarata:0.0109893696,((((((Typhlops_brongersmianus:0.0627900103,(Bothrops_moojeni:0.0195861198,Tantilla_melanocephala:0.0209819421)100_100
:0.0383223313)100_100:0.0368754732,((Anolis_meridionalis:0.0010448896,Anolis_brasiliensis:0.0006709941)100_100:0.0370710223,Tropidurus_oreadicus:0.0308
643818)100_100:0.0240763187)100_100:0.0015088974,Ophisaurus_gracilis:0.7153691699)100_100:0.0026577726,(((Gymnodactylus_amarali:0.0251972113,Hemidactyl
us_mabouia:0.0312956575)100_100:0.0318486216,Outgroup_outgroup:0.1736770065)100_100:0.0043661536,(Brasiliscincus_heathi:0.0039027396,(Copeoglossum_nigr
opunctatum:0.0029673186,Notomabuya_frenata:0.0031261982)100_100:0.0006163549)100_100:0.0449982559)100_100:0.0066837867)100_100:0.0050993046,Amphisbaena
_alba:0.0618370986)100_100:0.0400059809,(Colobosaura_modesta:0.0216100830,Micrablepharus_maximiliani:0.0277099824)100_100:0.0264979275)100_100:0.018433
6724,Ameivula_mumbuca:0.0127196352);"
t <- read.tree(text=z)
x <- root(t,outgroup="Outgroup_outgroup",resolve.root = TRUE)
y <- ladderize(x, right = TRUE)
plot(y)

AI <- "((homo_sapiens, chrysemys_picta, alligator_mississippiensis, gallus_gallus, sphenodon_sp), (anelytropsis_sp, dibamus_sp), (strophurus_sp, gonatodes_sp, saltorius_sp, gekko_sp, coleonyx_variegatus, lialis_sp), (cordylosaurus_sp, cordylus_sp, cricosaura_sp, lepidophyma_sp, plestiodon_fasciatus, tiliqua_sp, acontias_sp), (aspidoscelis_tigris, tupinambus_sp, pholidobollus_sp), (bipes_sp, rhineura_sp), ((typhlops_jamaicensis, micrurus_fulvius, python_molurus_2), ((anolis_carolinensis_2, hydrosaurus_sp, uta_stansburiana), (lanthanotus_borneensis, anniella_pulchra, heloderma_suspectum, varanus_exanthematicus, xenosaurus_platyceps))));"
SA <- "((homo_sapiens, chrysemys_picta, alligator_mississippiensis, gallus_gallus, sphenodon_sp), (anelytropsis_sp, dibamus_sp), (strophurus_sp, gonatodes_sp, saltorius_sp, gekko_sp, coleonyx_variegatus, lialis_sp), (cordylosaurus_sp, cordylus_sp, cricosaura_sp, lepidophyma_sp, plestiodon_fasciatus, tiliqua_sp, acontias_sp), (aspidoscelis_tigris, tupinambus_sp, pholidobollus_sp), (bipes_sp, rhineura_sp), (((typhlops_jamaicensis, micrurus_fulvius, python_molurus_2), (lanthanotus_borneensis, anniella_pulchra, heloderma_suspectum, varanus_exanthematicus, xenosaurus_platyceps)), (anolis_carolinensis_2, hydrosaurus_sp, uta_stansburiana)));"
SI <- "((homo_sapiens, chrysemys_picta, alligator_mississippiensis, gallus_gallus, sphenodon_sp), (anelytropsis_sp, dibamus_sp), (strophurus_sp, gonatodes_sp, saltorius_sp, gekko_sp, coleonyx_variegatus, lialis_sp), (cordylosaurus_sp, cordylus_sp, cricosaura_sp, lepidophyma_sp, plestiodon_fasciatus, tiliqua_sp, acontias_sp), (aspidoscelis_tigris, tupinambus_sp, pholidobollus_sp), (bipes_sp, rhineura_sp), (((typhlops_jamaicensis, micrurus_fulvius, python_molurus_2), (anolis_carolinensis_2, hydrosaurus_sp, uta_stansburiana)), (lanthanotus_borneensis, anniella_pulchra, heloderma_suspectum, varanus_exanthematicus, xenosaurus_platyceps)));" 
SC <- "((homo_sapiens, chrysemys_picta, alligator_mississippiensis, gallus_gallus, sphenodon_sp), (anolis_carolinensis_2, hydrosaurus_sp, uta_stansburiana), ((anelytropsis_sp, dibamus_sp), (strophurus_sp, gonatodes_sp, saltorius_sp, gekko_sp, coleonyx_variegatus, lialis_sp), (cordylosaurus_sp, cordylus_sp, cricosaura_sp, lepidophyma_sp, plestiodon_fasciatus, tiliqua_sp, acontias_sp), (aspidoscelis_tigris, tupinambus_sp, pholidobollus_sp), (bipes_sp, rhineura_sp), (typhlops_jamaicensis, micrurus_fulvius, python_molurus_2), (lanthanotus_borneensis, anniella_pulchra, heloderma_suspectum, varanus_exanthematicus, xenosaurus_platyceps)));"

t <- read.tree(text=SC)

plot(t)


Outgroup<- c("Gallus_gallus")
Gekkota<- c("Hemidactylus_mabouia","Gymnodactylus_amarali")
Scincoidea<- c("Brasiliscincus_heathi", "Copeoglossum_nigropunctatum", "Notomabuya_frenata")
Gymnothalmidae <- c("Colobosaura_modesta", "Micrablepharus_maximiliani")
Teioidea<- c("Ameiva_ameiva", "Ameivula_mumbuca", "Kentropyx_calcarata")
Amphisbaenia<-c("Amphisbaena_alba")
Iguania<- c("Anolis_brasiliensis", "Anolis_meridionalis", "Tropidurus_oreadicus")
Serpentes<- c("Corallus_hortulanus", "Chironius_exoletus", "Tantilla_melanocephala", "Sibynomorphus_mikanii", "Philodryas_olfersii", "Bothrops_lutzi", "Bothrops_moojeni", "Bothrops_pauloensis", "Liotyphlops_ternetzii", "Micrurus_brasiliensis", "Trilepida_brasiliensis", "Typhlops_brongersmianus")
Boidae <- c("Corallus_hortulanus")
Colubridae <- c("Chironius_exoletus", "Tantilla_melanocephala")
Dipsadidae <- c("Sibynomorphus_mikanii", "Philodryas_olfersii")
Viperidae <- c("Bothrops_lutzi", "Bothrops_moojeni", "Bothrops_pauloensis")
Anomalepididae <- c("Liotyphlops_ternetzii")
Elapidae <- c("Micrurus_brasiliensis")
Leptotyphlopidae <- c("Trilepida_brasiliensis")
Typhlopidae <- c("Typhlops_brongersmianus")
Toxicofera <- c("Anolis_brasiliensis", "Anolis_meridionalis", "Tropidurus_oreadicus", "Corallus_hortulanus", "Chironius_exoletus", "Tantilla_melanocephala", "Sibynomorphus_mikanii", "Philodryas_olfersii", "Bothrops_lutzi", "Bothrops_moojeni", "Bothrops_pauloensis", "Liotyphlops_ternetzii", "Micrurus_brasiliensis", "Trilepida_brasiliensis", "Typhlops_brongersmianus")


