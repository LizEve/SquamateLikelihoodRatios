#! /usr/bin/env python
from Bio import AlignIO
from Bio import SeqIO
from Bio.Alphabet import generic_dna
import os
import glob
import shutil 
import itertools
import re 
import pandas as pd


# Make folders for each locus and put nexus file into folder
def makeFolders(n, mainDir):

	# Get locus name 
	locus=os.path.splitext(n)[0]

	# Create paths 
	dirPath = os.path.join(mainDir,locus)

	# Make directories
	if not os.path.exists(dirPath):
		os.mkdir(dirPath)	

	# Copy nexus to data folder
	os.system("cp %s %s" % (n,dirPath))

	# Gut check 
	print("File "+ n + " moved to "+dirPath)
	return dirPath,locus

# List missing taxa from MSA file 
def missingTax(nex,allTaxa):

	# Open MSA file and get list of taxa
	alignment = AlignIO.read(open(nex),"nexus")
	locusTaxa=[]
	for record in alignment:
		locusTaxa.append(str(record.id))

	# Get list of all taxa not in locus 
	missingTaxa=list(allTaxa.difference(set(locusTaxa)))
	return missingTaxa 

def editBayesBlock(suffix,fName,missingTaxa,dumpFolder,dirPath,mainDir):
	# Keep track of which constraints are removed for locus 
	removedConstraints=[]

	# Open template bb file and new locus specific bb file 
	bbin=open('tacocat'+suffix)
	bbout=open(fName+suffix,"w+")

	# For each line in bayes block
	for line in bbin:
		# Replace tacocat with locus name 
		line=line.replace("tacocat",fName)

		# Remove missing taxa, replace with nothing, then remove double commas
		for t in missingTaxa:
			line=line.replace(" "+t,"",20)
			line=line.replace(",,",",",20)
		
		# Flag all empty constraints. and store in list 
		if "=;" in line:
			
			# Grab name of constraint to remove from prset topology constraints 
			emptyConstraint=line.split(" ")[-3]
			removedConstraints.append(emptyConstraint)
		else:
			# If empty constraint is in prset line, replace with nothing
			for c in removedConstraints:
				line=line.replace(c,"")

			# clean up prset and other lines that might have extra commas
			line = line.replace(", ,",",",20)
			line = line.replace(",,",",",20)
			line = line.replace(",)",")",20)
			line = line.replace(", )",")",20)
			line = line.replace("(, ","(",20)
			line = line.replace("(,","(",20)

			# Write line to new file. 
			bbout.write(line)
	bbin.close()
	
	return removedConstraints


def setup(mainDir,allTaxa,bbFiles,outFile,reqClades,folders = True):
	
	os.chdir(mainDir)

	# Make folder to dump nexus files with missing clades. Delete manually later, didn't want to write a delete folder call into a loop in case things go wrong.  
	dumpFolder=os.path.join(mainDir,"missingCladeDir")
	if not os.path.exists(dumpFolder):
		os.mkdir(dumpFolder)	

	# Keep track of which clades are missing for each locus
	nexMissing={}

	# Keep track of which loci are kept 
	nexList=[]

	# For each nexus file
	for nex in glob.glob('*nex'):

		# Make folder for locus 
		if folders == True:
			
			# Make a folder and move MSA file into folder 
			dirPath,fName = makeFolders(nex, mainDir)

		else:
			# Get locus name 
			fName=os.path.split(n)[0]
			# Create paths 
			dirPath = os.path.join(mainDir,fName)

		print(dirPath)
		# Move into folder and grab bb file 
		os.chdir(dirPath)

		os.system("cp ../*.bb .")
		
		# Open MSA and get taxa that are missing
		missingTaxa = missingTax(nex,allTaxa)
		
		# counter
		x=0

		# keep track of required clades that are missing 
		reqMissingClades=[]

		# Iterate through bayesblock files 
		for suffix in bbFiles:

			# For the first bayesblock, check if any required clades are missing
			if x == 0:

				# Edit BayesBlock and return list of clades that have no taxa in them. 
				missingClades=editBayesBlock(suffix,fName,missingTaxa,dumpFolder,dirPath,mainDir)
				
				# Get the set of requires clades that are missing
				bad = set(missingClades).intersection(set(reqClades))

				# If there are clades missing, add 10 to x to pass down to a section to move files around 
				if not len(bad) == 0:
					# Get which clades are missing for locus
					reqMissingClades=list(bad)
					x += 10

				# If there are no missing clades, add one to counter and keep going with bayesblocks
				elif len(bad) == 0:
					
					# Edit BayesBlock and return list of clades that have no taxa in them. 
					missingClades=editBayesBlock(suffix,fName,missingTaxa,dumpFolder,dirPath,mainDir)
					nexList.append(nex)
					x += 1
			
			# For following bb files when there are no missing clades
			elif 0 < x < 10:

				# Edit BayesBlock and return list of clades that have no taxa in them. 
				missingClades=editBayesBlock(suffix,fName,missingTaxa,dumpFolder,dirPath,mainDir)

			# Move whole folder and nexus file 
			elif x >= 10:
				# Move entire folder to compost file for manual deletion
					try:
						shutil.move(dirPath,dumpFolder)
						shutil.move(os.path.join(mainDir,nex),dumpFolder)
					except:
						print(NameError)
		
		# When done with locus 
		os.system("rm ./tacocat*")
		os.chdir(mainDir)
			
		# Save in dictionary
		if len(reqMissingClades) > 0:
			nexMissing[nex]=reqMissingClades

	with open(os.path.join(mainDir,"nexFilesSetup.txt"), 'w') as f:
			for z in nexList:
				f.write("%s \n" % z)

	# print dictionary to file 
	pd.DataFrame.from_dict(data=nexMissing, orient='index').to_csv(outFile, header=False)

def main():
	# Assums you have a folder full of nexus files, and template bb files. 
	mainDir='/Users/ChatNoir/Projects/Squam/Burbrink/BurbrinkSS/'
	bbFiles = ["_sclero.bb","_toxpoly.bb","_toxai.bb","_toxsa.bb","_toxsi.bb"]
	outFile = '/Users/ChatNoir/Projects/Squam/Burbrink/missingClades.csv'
	# set of all taxa, needs to be a set and not a list or tuple 
	allTaxa={'Amphisbaenidae_Amphisbaena_fuliginosa_I1051', 'Amphisbaenidae_Geocalamus_acutus_I8584', 'Bipedidae_Bipes_biporus_I6291', 'Bipedidae_Bipes_canaliculatus_I7475', 'Trogonophiidae_Trogonophis_wiegmanni_I8588', 'Anguidae_Elgaria_multicarinata_I7480', 'Anguidae_Ophisaurus_apodus_I8102', 'Anniellidae_Anniella_pulchra_I8586', 'Anguidae_Celestus_enneagrammus_I8580', 'Helodermatidae_Heloderma_horridum_I1046', 'Helodermatidae_Heloderma_suspectum_I1047', 'Varanidae_Varanus_acanthurus_I1007', 'Varanidae_Varanus_bengalensis_I0970', 'Varanidae_Varanus_exanthematicus_I1827', 'Varanidae_Varanus_griseus_I0971', 'Varanidae_Varanus_salvator_I1914', 'Xenosauridae_Xenosaurus_grandis_I3684', 'Xenosauridae_Xenosaurus_platyceps_I1021', 'Carphodactylidae_Saltuarius_cornutus_I20515', 'Diplodactylidae_Rachodactylus_leachianus_I20516', 'Diplodactylidae_Strophurus_ciliaris_I20517', 'Eublepharidae_Coleonyx_variegatus_I1048', 'Eublepharidae_Eublepharis_macularius_I1042', 'Gekkonidae_Cnemaspis_psychedelica_I3182', 'Gekkonidae_Cyrtodactylus_pulchellus_I3192', 'Gekkonidae_Gekko_gecko_I1041', 'Gekkonidae_Phelsuma_lineata_I7483', 'Pygopodidae_Delma_borea_I6290', 'Pygopodidae_Lialis_burtonis_I8124', 'Sphaerodactylidae_Gonatodes_albogularis_I1037', 'Sphaerodactylidae_Teratoscincus_sp_I8838', 'Gymnophthalmidae_Colobosaura_modesta_I1010', 'Gymnophthalmidae_Pholidobolus_montium_I8840', 'Agamidae_Acanthosaura_crucigera_I3191', 'Agamidae_Agama_agama_I1038', 'Agamidae_Calotes_emma_I1012', 'Agamidae_Gonocephalus_zamboanga_I3172', 'Agamidae_Leiolepis_belliana_I3180', 'Agamidae_Physignathus_cocincinus_I1015', 'Agamidae_Pogona_vitticeps_I1018', 'Agamidae_Uromastyx_sp_I1049', 'Chamaeleonidae_Brookesia_brygooi_I7477', 'Chamaeleonidae_Chamaeleo_calyptratus_I1026', 'Corytophanidae_Basiliscus_basiliscus_I8591', 'Corytophanidae_Basiliscus_plumifrons_I3905', 'Corytophanidae_Corytophanes_cristatus_I12759', 'Crotaphytidae_Crotaphytus_collaris_I6766', 'Dactyloidae_Anolis_carolinensis_I8301', 'Hoplocercidae_Enyalioides_sp_I8839', 'Hoplocercidae_Hoplocercus_spinosus_I12271', 'Iguanidae_Brachylophus_fasciatus_I7476', 'Iguanidae_Dipsosaurus_dorsalis_I1039', 'Iguanidae_Sauromalus_ater_I8350', 'Leiocephalidae_Leiocephalus_barahonensis_I8842', 'Liolaemidae_Ctenoblepharis_adspersa_I3171', 'Liolaemidae_Liolaemus_signifer_I3170', 'Liolaemidae_Phymaturus_palluma_I3130', 'Opluridae_Chalarodon_madagascariensis_I7478', 'Opluridae_Oplurus_cyclurus_I7482', 'Phrynosomatidae_Petrosaurus_mearnsi_I8348', 'Phrynosomatidae_Phrynosoma_cornutum_I1020', 'Phrynosomatidae_Sceloporus_magister_I7356', 'Phrynosomatidae_Uma_scoparia_I8104', 'Phrynosomatidae_Uta_stansburiana_I7484', 'Polychrotidae_Polychrus_marmoratus_I1011', 'Tropiduridae_Plica_plica_I1034', 'Tropiduridae_Stenocercus_guentheri_I3907', 'Tropiduridae_Uranoscodon_superciliosus_I1036', 'Lacertidae_Takydromus_ocellatus_I8176', 'Lacertidaea_Lacerta_agilis_I1016', 'Leiosauridae_Leiosaurus_bellii_I12269', 'Leiosauridae_Pristidactylus_scapulatus_I8583', 'Leiosauridae_Urostrophus_vautieri_I1033', 'Cordylidae_Cordylus_sp_I1043', 'Cordylidae_Platysaurus_imperator_I8339', 'Gerrhosauridae_Zonosaurus_ornatus_I7486', 'Scincidae_Acontias_percivali_I8346', 'Scincidae_Amphiglossus_splendidus_I7474', 'Scincidae_Brachymeles_gracilis_I8837', 'Scincidae_Eugongylus_rufescens_I8126', 'Scincidae_Eutropis_multifasciata_I3183', 'Scincidae_Feylinia_polylepis_I8099', 'Scincidae_Mabuya_quinquetaeniata_I7481', 'Scincidae_Scincus_scincus_I7357', 'Scincidae_Sphenomorphus_solomonis_I8125', 'Scincidae_Tiliqua_rugosa_I5874', 'Xantusiidae_Lepidophyma_flavimaculatum_I1031', 'Xantusiidae_Xantusia_henshawi_I8351', 'Teiidae_Ameiva_ameiva_I0533', 'Teiidae_Aspidoscelis_deppii_I3146', 'Teiidae_Callopistes_maculatus_I8343', 'Teiidae_Cnemidophorus_lemniscatus_I0581', 'Teiidae_Cnemidophorus_ocellifer_I0588', 'Teiidae_Dicrodon_heterolepis_I3160', 'Teiidae_Dracaena_guianensis_I0594', 'Teiidae_Kentropyx_calcarata_I0599', 'Teiidae_Teius_teyou_I3135', 'Teiidae_Tupinambis_teguixin_I0524', 'Tropidophiidae_Trachyboa_boulengeri_I1341', 'Tropidophiidae_Tropidophis_feicki_I7426', 'Tropidophiidae_Tropidophis_haetianus_I9380', 'Tropidophiidae_Tropidophis_maculatus_I7430', 'Tropidophiidae_Tropidophis_paucisquamis_I1308', 'Dibamidae_Dibamus_novaeguineae_I6292', 'Shinisauridae_Shinisaurus_crocodilurus_I1045', 'Sphenodontidae_Sphenodon_punctatus_I5870', 'Acrochordidae_Acrochordus_granulatus_I7464', 'Acrochordidae_Acrochordus_javanicus_I1320', 'Anomalepididae_Liotyphlops_albirostris_I8450', 'Anomalepididae_Liotyphlops_beui_I1305', 'Anomalepididae_Typhlophis_squamosus_I1350', 'Boidae_Boa_constrictor_I1307', 'Boidae_Corallus_caninus_I1306', 'Boidae_Corallus_cropanii_I1303', 'Boidae_Corallus_hortulanus_I1358', 'Boidae_Corallus_ruschembergerii_I1005', 'Boidae_Epicrates_alvarezi_I7462', 'Boidae_Epicrates_angulifer_I4004', 'Boidae_Epicrates_assisi_I1348', 'Boidae_Epicrates_cenchria_I7458', 'Boidae_Epicrates_crassus_I7463', 'Boidae_Epicrates_striatus_I13179', 'Boidae_Eunectes_murinus_I1374', 'Erycidae_Calabaria_reinhardtii_I4001', 'Boidae_Candoia_aspera_I1854', 'Boidae_Candoia_bibroni_I1312', 'Erycidae_Charina_bottae_I4003', 'Erycidae_Lichanura_trivirgata_I8100', 'Erycidae_Eryx_colubrinus_I1326', 'Erycidae_Eryx_conicus_I1323', 'Erycidae_Eryx_jaculus_I4009', 'Erycidae_Eryx_jayakari_I1324', 'Erycidae_Eryx_johnii_I1325', 'Erycidae_Eryx_miliaris_I4005', 'Erycidae_Eryx_muelleri_I4006', 'Erycidae_Eryx_tataricus_I4007', 'Boidae_Acrantophis_dumerilii_I3999', 'Boidae_Sanzinia_madagascariensis_I1328', 'Ungaliophiidae_Exiliboa_placata_I9382', 'Ungaliophiidae_Ungaliophis_panamensis_I1345', 'Calamaridae_Calamaria_parvimentata_I20524', 'Calamaridae_Calamaria_septentrionalis_I20519', 'Calamaridae_Calamaria_suluensis_I20518', 'Calamaridae_Calamaria_yunnanensis_I20526', 'Colubridae_Ahaetulla_prasina_I12180', 'Colubridae_Boiga_angulata_I12181', 'Colubridae_Coluber_constrictor_I7469', 'Colubridae_Dasypeltis_scabra_I12189', 'Colubridae_Dendrelaphis_pictus_I12190', 'Colubridae_Gonyosoma_oxycephalum_I12196', 'Colubridae_Lampropeltis_getula_I0101', 'Colubridae_Lycodon_aulicus_I12201', 'Colubridae_Oligodon_cyclurus_I12204', 'Colubridae_Oxybelis_aeneus_I12205', 'Colubridae_Rhynchocalamus_melanocephalus_I12208', 'Colubridae_Simophis_rhinostoma_I12212', 'Dipsadidae_Alsophis_cantherigerus_I7405', 'Dipsadidae_Clelia_equatoriana_I1006', 'Dipsadidae_Diadophis_punctatus_I1126', 'Dipsadidae_Dipsas_catesbyi_I1078', 'Dipsadidae_Elapomorphus_quinquelineatus_I1088', 'Dipsadidae_Erythrolamprus_bizona_I7533', 'Dipsadidae_Farancia_erytrogramma_I1125', 'Dipsadidae_Geophis_godmani_I7535', 'Dipsadidae_Heterodon_nasicus_I1224', 'Dipsadidae_Hydrops_triangularis_I1061', 'Dipsadidae_Hypsiglena_torquata_I7540', 'Dipsadidae_Thermophis_baileyi_I1241', 'Grayiidae_Grayia_ornata_I12197', 'Homalopsidae_Bitia_hydroides_I12239', 'Homalopsidae_Cantoria_violacea_I12240', 'Homalopsidae_Enhydris_enhydris_I12241', 'Homalopsidae_Fordonia_leucobalia_I12242', 'Homalopsidae_Gerarda_prevostiana_I12243', 'Homalopsidae_Homalopsis_buccata_I12244', 'Homalopsidae_Myrrophis_chinensis_I12245', 'Natricidae_Afronatrix_anoscopus_I12267', 'Natricidae_Amphiesma_stolata_I8096', 'Natricidae_Macropisthodon_rudis_I12255', 'Natricidae_Natrix_natrix_I6769', 'Natricidae_Seminatrix_pygaea_I12264', 'Natricidae_Sinonatrix_percarinata_I12257', 'Natricidae_Storeria_dekayi_I12258', 'Natricidae_Thamnophis_marcianus_I7473', 'Natricidae_Tropidoclonion_lineatum_I12259', 'Natricidae_Xenochrophis_psicator_I6770', 'Pareatidae_Pareas_margaritophorus_I8103', 'Pareatidae_Pareas_nuchalis_I7472', 'Pseudoxenodontidae_Pseudoxenodon_bambusicola_I20522', 'Pseudoxenodontidae_Pseudoxenodon_karlschmidti_I20525', 'Pseudoxenodontidae_Pseudoxenodon_macrops_I20523', 'Sibynophiidae_Scaphiodontophis_annulatus_I12210', 'Sibynophiidae_Sibynophis_collaris_I12211', 'Xenodermatidae_Achalinus_rufescens_I6190', 'Xenodermatidae_Xenodermus_javanicus_I8177', 'Cylindrophiidae_Cylindrophis_maculatus_I1318', 'Cylindrophiidae_Cylindrophis_ruffus_I1330', 'Atractaspididae_Aparallactus_capensis_I12176', 'Atractaspididae_Aparallactus_werneri_I7466', 'Atractaspididae_Atractaspis_bibronii_I7467', 'Atractaspididae_Atractaspis_engaddensis_I12177', 'Atractaspididae_Homoroselaps_lacteus_I12178', 'Atractaspididae_Micrelaps_muelleri_I12179', 'Colubroides_Cyclocorus_lineatus_I12220', 'Elapidae_Acanthophis_pyrrhus_I12221', 'Elapidae_Bungarus_multicinctus_I12224', 'Elapidae_Calliophis_maculiceps_I12225', 'Elapidae_Dendroaspis_viridis_I12226', 'Elapidae_Laticauda_colubrina_I5831', 'Elapidae_Micrurus_fulvius_I8101', 'Elapidae_Naja_naja_I5829', 'Elapidae_Notechis_scutatus_I5840', 'Elapidae_Ophiophagus_hannah_I12231', 'Elapidae_Pelamis_platurus_I12233', 'Elapidae_Pseudonaja_nuchalis_I12235', 'Elapidae_Sinomicrurus_macclellandi_I12236', 'Elapoidea_Psammodynastes_pulverulentus_I12238', 'Lamprophiidae_Bothrophthalmus_sp_I12246', 'Lamprophiidae_Gonionotophis_brussauxi_I12249', 'Lamprophiidae_Gonionotophis_poensis_I12175', 'Lamprophiidae_Hormonotus_modestus_I12250', 'Lamprophiidae_Mehelya_guirali_I6187', 'Psammophiidae_Mimophis_mahfalensis_I12261', 'Psammophiidae_Rhamphiophis_rubripunctatus_I0181', 'Pseudoxyrhophiidae_Dromicodryas_bernieri_I12247', 'Pseudoxyrhophiidae_Duberria_lutrix_I12248', 'Pseudoxyrhophiidae_Duberria_variegata_I17939', 'Pseudoxyrhophiidae_Ithycyphus_miniatus_I12251', 'Pseudoxyrhophiidae_Liophidium_torquatum_I12252', 'Pseudoxyrhophiidae_Liopholidophis_rhadinaea_I17931', 'Pseudoxyrhophiidae_Lycodryas_inornatus_I17932', 'Pseudoxyrhophiidae_Madagascarophis_colubrinus_I12254', 'Pseudoxyrhophiidae_Parastenophis_betsileanus_I17930', 'Pseudoxyrhophiidae_Phisalixella_sp_I17938', 'Gerrhopilidae_Gerrhopilus_hades_I8123', 'Gerrhopilidae_Gerrhopilus_inornatus_I8122', 'Leptotyphlopidae_Epictia_tenella_I1349', 'Leptotyphlopidae_Rena_dulcis_I16786', 'Leptotyphlopidae_Siagonodon_septemstriatus_I12268', 'Leptotyphlopidae_Trilepida_fuliginosa_I1354', 'Leptotyphlopidae_Trilepida_koppesi_I1309', 'Leptotyphlopidae_Trilepida_macrolepis_I1304', 'Leptotyphlopidae_Trilepida_salgueiroi_I1311', 'Pythonidae_Antaresia_childreni_I1839', 'Pythonidae_Antaresia_maculosa_I4000', 'Pythonidae_Antaresia_perthensis_I1868', 'Pythonidae_Antaresia_stimsoni_I1866', 'Pythonidae_Apodora_papuana_I3929', 'Pythonidae_Aspidites_melanocephalus_I1843', 'Pythonidae_Aspidites_ramsayi_I1332', 'Pythonidae_Bothrochilus_boa_I1313', 'Pythonidae_Leiopython_hoserae_I1901', 'Pythonidae_Liasis_fuscus_I2447', 'Pythonidae_Liasis_mackloti_I4008', 'Pythonidae_Malayopython_reticulatus_I1917', 'Pythonidae_Morelia_amethistina_I2434', 'Pythonidae_Morelia_bredli_I2445', 'Pythonidae_Morelia_oenpelliensis_I1844', 'Pythonidae_Morelia_spilota_I5869', 'Pythonidae_Morelia_viridis_I3914', 'Pythonidae_Python_bivittatus_I12758', 'Pythonidae_Python_curtus_I1316', 'Pythonidae_Python_regius_I1352', 'Pythonidae_Python_sebae_I2437', 'Typhlopidae_Amerotyphlops_amoipira_I1781', 'Typhlopidae_Amerotyphlops_brongersmianus_I1776', 'Typhlopidae_Amerotyphlops_reticulatus_I1782', 'Typhlopidae_Indiotyphlops_braminus_I1343', 'Typhlopidae_Indotyphlops_pushpakumara_I0833', 'Typhlopidae_Rhinotyphlops_lalandei_I1342', 'Typhlopidae_Typhlops_biminiensis_I7445', 'Typhlopidae_Typhlops_lumbricalis_I7446', 'Uropeltidae_Melanophidium_punctatum_I1346', 'Uropeltidae_Pseudotyphlops_philippinus_I0822', 'Uropeltidae_Rhinophis_blythii_I0823', 'Uropeltidae_Rhinophis_homolepis_I0824', 'Uropeltidae_Uropeltis_melanogaster_I8448', 'Viperidae_Agkistrodon_contortrix_I2093', 'Viperidae_Atropoides_nummifer_I1984', 'Viperidae_Azemiops_feae_I1982', 'Viperidae_Bothriechis_nigroviridis_I1980', 'Viperidae_Bothrops_asper_I1960', 'Viperidae_Causus_maculatus_I2001', 'Viperidae_Daboia_russellii_I0796', 'Viperidae_Lachesis_muta_I1956', 'Aniliidae_Anilius_scytale_I1310', 'Bolyeriidae_Casarea_dussumieri_I8815', 'Loxocemidae_Loxocemus_bicolor_I1331', 'Xenopeltidae_Xenopeltis_unicolor_I3916'}
	# There isnt an outgroup constraint in the bb, because it is only one taxa, but this means no genes got kicked out when they should have. 
	reqClades=['Iguania','Anguimorpha','Snake','Outgroup']
	setup(mainDir,allTaxa,bbFiles,outFile,reqClades)

if __name__ == "__main__":
    main()
