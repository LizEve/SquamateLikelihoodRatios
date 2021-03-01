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
mainDir='/Users/ChatNoir/Projects/Squam/Singhal/SinghalData'
# set of all taxa 
Taxa={'Leptodeira_annulata', 'Erythrolamprus_reginae', 'Notomabuya_frenata', 'Erythrolamprus_almadensis', 'Lygophis_paucidens', 'Bothrops_pauloensis', 'Copeoglossum_nigropunctatum', 'Psomophis_joberti', 'Philodryas_olfersii', 'Ameiva_ameiva', 'Amphisbaena_alba', 'Sphenodon_punctatus2', 'Sibynomorphus_mikanii', 'Oxyrhopus_trigeminus', 'Trilepida_brasiliensis', 'Salvator_meriana', 'Imantodes_cenchoa', 'Philodryas_nattereri', 'Alligator_mississippiensis', 'Bothrops_lutzi', 'Taeniophallus_occipitalis', 'Crotalus_horridus', 'Gymnodactylus_amarali', 'Brasiliscincus_heathi', 'Pseudoboa_neuwiedii', 'Zootoca_vivipara', 'Gopherus_evgoodei', 'Pogona_vitticeps', 'Erythrolamprus_poecilogyrus', 'Lacerta_bilineata', 'Chironius_exoletus', 'Anolis_carolinensis', 'Podarcis_muralis', 'Typhlops_brongersmianus', 'Laticauda_colubrina', 'Eublepharis_macularius', 'Thamnodynastes_hypoconia', 'Varanus_komodoensis', 'Anolis_meridionalis', 'Pseudoboa_nigra', 'Hemidactylus_mabouia', 'Xenodon_merremi', 'Phimophis_guerini', 'Tropidurus_oreadicus', 'Apostolepis_polylepis', 'Shinisaurus_crocodilurus', 'Ophiophagus_hannah', 'Colobosaura_modesta', 'Sphenodon_punctatus1', 'Apostolepis_cearensis', 'Corallus_hortulanus', 'Hydrophis_melanocephalus', 'Liotyphlops_ternetzii', 'Micrurus_brasiliensis', 'Gallus_gallus', 'Protobothrops_mucrosquamatus', 'Micrablepharus_maximiliani', 'Oxyrhopus_petolarius', 'Lacerta_viridis', 'Paroedura_picta', 'Kentropyx_calcarata', 'Deinagkistrodon_actutus', 'Python_bivittatus', 'Dopasia_gracilis', 'Ameivula_mumbuca', 'Tantilla_melanocephala', 'Anolis_brasiliensis', 'Chrysemys_picta', 'Bothrops_moojeni', 'Gekko_japonicus'}
bbFiles = ["_toxpoly.bb","_sclero.bb","_toxai.bb","_toxsa.bb","_toxsi.bb"]
setup(mainDir,suffix,Taxa,bbFiles)