#! /usr/bin/env python
from Bio import AlignIO
from Bio import SeqIO
from Bio.Alphabet import generic_dna
import os
import glob
import shutil 
import itertools
import re 

# Burbrink Version
# requires a main folder with: nexus files, and 5 template bayesblocks 

# Make folders for each locus and put nexus file into folder
def makeFolders(n, suffix, mainDir):

	# get file name 
	nexus=str(n)
	nexusIndex = nexus.find(suffix)
	fName = nexus[:nexusIndex]

	# Create paths 
	dirPath = os.path.join(mainDir,fName)

	# Make directories
	if not os.path.exists(dirPath):
		os.mkdir(dirPath)	

	# Copy nexus to data folder
	os.system("cp %s %s" % (fName+".nex",dirPath))

	# Gut check 
	print("File "+ n + " moved")
	return dirPath,fName


def editFile(fName,missingTaxa,bbFiles):
	for suffix in bbFiles:
		# Open template bb file and new locus specific bb file 
		bbin=open('tacocat'+suffix)
		bbout=open(fName+suffix,"w+")
		removedConstraints=[]
		# For each line in template, remove missing taxa, replace tacocat with locus
		for line in bbin:
			# Replace tacocat with locus name 
			line=line.replace("tacocat",fName)
			# Remove all missing taxa, replace with nothing, then remove double commas
			for t in missingTaxa:
				line=line.replace(" "+t,"",20)
				line=line.replace(",,",",",20)
			# Remove all empty constraints. and store in list 
			if "=;" in line:
				print(line)
				# Grab name of constraint to remove from prset topology constraints 
				emptyConstraint=line.split(" ")[-3]
				removedConstraints.append(emptyConstraint)
				print(emptyConstraint)
			else:
				# if constraint in prset line, replace with nothing
				for c in removedConstraints:
					line=line.replace(c,"")
				# clean up prset and other lines that might have extra commas
				line = line.replace(", ,",",",20)
				line = line.replace(",,",",",20)
				line = line.replace(",)",")",20)
				line = line.replace(", )",")",20)
				bbout.write(line)
		
		bbin.close()
		bbout.close()


def setup(mainDir,suffix,allTaxa,bbFiles,folders = True):
	os.chdir(mainDir)
	# For each nexus file, make a folder, and get a path to folder
	for nex in glob.glob('*%s' % suffix):
		if folders == True:
			# Make a folder and move MSA file into folder 
			dirPath,fName = makeFolders(nex, suffix, mainDir)
		else:
			nexus=str(nex)
			nexusIndex = nexus.find(suffix)
			fName = nexus[:nexusIndex]
			dirPath = os.path.join(mainDir,fName)
		# Move into folder and grab bb file 
		os.chdir(dirPath)
		os.system("cp ../*.bb .")
		# Open MSA file and get list of missing taxa for locus
		alignment = AlignIO.read(open(nex),"nexus")
		locusTaxa=[]
		for record in alignment:
			locusTaxa.append(str(record.id))
		# Get list of all taxa not in locus 
		missingTaxa=list(allTaxa.difference(set(locusTaxa)))
		editFile(fName,missingTaxa,bbFiles)
		os.system("rm tacocat*.bb")
		os.chdir(mainDir)


suffix='.nex'
mainDir='/Users/ChatNoir/Projects/Squam/Reeder/ReederData'
#mainDir='/work/gmount/ReederSS'
# set of all taxa 
Taxa={'Notechis_scutatus', 'Feylinia_polylepis', 'Lanthanotus_borneensis', 'Homalopsis_buccata', 'Celestus_enneagrammus', 'Cricosaura_typica', 'Rhineura_floridana', 'Acrochordus_granulatus', 'Sphenomorphus_solomonis', 'Zonosaurus_ornatus', 'Eugongylus_rufescens', 'Heloderma_horridum', 'Bipes_biporus', 'Aeluroscalobates_felinus', 'Pseudopus_apodus', 'Gekko_gecko', 'Diadophis_punctatus', 'Morunasaurus_annularis', 'Tiliqua_scincoides', 'Cylindrophis_rufus', 'Tupinambis_teguixin', 'Xantusia_vigilis', 'Plica_plica', 'Chamaeleo', 'Rhacodactylus_auriculatus', 'Xenosaurus_grandis', 'Boa_constrictor', 'Amphiesma_stolata', 'Thamnophis_marcianus', 'Corytophanes_cristatus', 'Strophurus_ciliaris', 'Diplometopon_zarudnyi', 'Saltuarius_cornutus', 'Trachyboa_boulengeri', 'Pholidobolus', 'Oplurus_cyclurus', 'Calotes_emma', 'Xenodermus_javanicus', 'Alligator', 'Trachylepis_quinquetaeniata', 'Xenochrophis_piscator', 'Dromaius', 'Amphisbaena_fuliginosa', 'Agkistrodon_contortrix', 'Takydromus_ocellatus', 'Dibamus_novaeguineae', 'Basiliscus_basiliscus', 'Anniella_pulchra', 'Varanus_acanthurus', 'Anolis_carolinensis', 'Gambelia_wislizenii', 'Enyalioides_laticeps', 'Lialis_burtonis', 'Anelytropsis_papillosus', 'Lichanura_trivirgata', 'Natrix_natrix', 'Brachymeles_gracilis', 'Callopistes_maculatus', 'Loxocemus_bicolor', 'Typhlops_jamaicensis', 'Podocnemis', 'Azemiops_feae', 'Ungaliophis_continentalis', 'Stenocercus_guentheri', 'Colobosaura_modesta', 'Leptotyphlops', 'Teius_teyou', 'Tropidophis_haetianus', 'Atractaspis_irregularis', 'Crotaphytus_collaris', 'Eryx_colubrinus', 'Naja', 'Uma_scoparia', 'Shinisaurus_crocodilurus', 'Liolaemus_bellii', 'Uta_stansburiana', 'Eublepharis_macularius', 'Elgaria_multicarinata', 'Phrynosoma_platyrhinos', 'Xenopeltis_unicolor', 'Phymaturus_palluma', 'Aspidites_melanocephalus', 'Uranoscodon_superciliosus', 'Leiocephalus_barahonensis', 'Heloderma_suspectum', 'Chelydra', 'Physignathus_cocincinus', 'Bothrops_asper', 'Petrosaurus_mearnsi', 'Laticauda_colubrina', 'Brachylophus_fasciatus', 'Platysaurus', 'Acontias', 'Pareas_hamptoni', 'Homo', 'Coleonyx_variegatus', 'Afronatrix_anoscopus', 'Micrurus_fulvius', 'Bipes_canaliculatus', 'Agama_agama', 'Trogonophis_wiegmanni', 'Mus', 'Leiosaurus_catamarcensis', 'Varanus_exanthematicus', 'Anilius_scytale', 'Coluber_constrictor', 'Lachesis_muta', 'Causus', 'Gallus', 'Teratoscincus', 'Pogona_vitticeps', 'Calabaria_reinhardtii', 'Geocalamus_acutus', 'Polychrus_marmoratus', 'Plestiodon_fasciatus', 'Brookesia_brygooi', 'Lycophidion_capense', 'Liotyphlops_albirostris', 'Casarea_dussumieri', 'Python_molurus', 'Urostrophus_vautieri', 'Sceloporus_variabilis', 'Amphiglossus_splendidus', 'Cordylus_mossambicus', 'Hydrosaurus', 'Xenosaurus_platyceps', 'Epicrates_striatus', 'Chalarodon_madagascariensis', 'Uropeltis_melanogaster', 'Gonatodes_albogularis', 'Varanus_salvator', 'Tachyglos', 'Lepidophyma_flavimaculatu', 'Lampropeltis_getula', 'Scincus', 'Daboia_russelli', 'Pristidactylus_torquatus', 'Leiolepis_belliana', 'Cordylosaurus_subtesselatus', 'Exiliboa_placata', 'Lacerta_viridis', 'Crocodylus', 'Sphenodon_punctatus', 'Delma_borea', 'Aparallactus_werneri', 'Aspidoscelis_tigris', 'Sauromalus_ater', 'Phelsuma_lineata', 'Uromastyx_aegyptus', 'Dipsosaurus_dorsalis'}
bbFiles = ["_toxpoly.bb","_sclero.bb","_toxai.bb","_toxsa.bb","_toxsi.bb"]
setup(mainDir,suffix,Taxa,bbFiles)