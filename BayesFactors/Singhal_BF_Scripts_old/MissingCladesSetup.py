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
	try:
		alignment = AlignIO.read(open(nex),"nexus")
		locusTaxa=[]
		for record in alignment:
			locusTaxa.append(str(record.id))

		# Get list of all taxa not in locus 
		missingTaxa=list(allTaxa.difference(set(locusTaxa)))
		return missingTaxa 
	except:
		w = 10
		return w

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
		
		if missingTaxa == 10:
			x=10
		else:
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
	mainDir='/Users/ChatNoir/Projects/Squam/Singhal/SinghalSS/'
	bbFiles = ["_sclero.bb","_toxpoly.bb","_toxai.bb","_toxsa.bb","_toxsi.bb"]
	outFile = '/Users/ChatNoir/Projects/Squam/Singhal/missingClades.csv'
	# set of all taxa, needs to be a set and not a list or tuple 
	allTaxa={'Leptodeira_annulata', 'Erythrolamprus_reginae', 'Notomabuya_frenata', 'Erythrolamprus_almadensis', 'Lygophis_paucidens', 'Bothrops_pauloensis', 'Copeoglossum_nigropunctatum', 'Psomophis_joberti', 'Philodryas_olfersii', 'Ameiva_ameiva', 'Amphisbaena_alba', 'Sphenodon_punctatus2', 'Sibynomorphus_mikanii', 'Oxyrhopus_trigeminus', 'Trilepida_brasiliensis', 'Salvator_meriana', 'Imantodes_cenchoa', 'Philodryas_nattereri', 'Alligator_mississippiensis', 'Bothrops_lutzi', 'Taeniophallus_occipitalis', 'Crotalus_horridus', 'Gymnodactylus_amarali', 'Brasiliscincus_heathi', 'Pseudoboa_neuwiedii', 'Zootoca_vivipara', 'Gopherus_evgoodei', 'Pogona_vitticeps', 'Erythrolamprus_poecilogyrus', 'Lacerta_bilineata', 'Chironius_exoletus', 'Anolis_carolinensis', 'Podarcis_muralis', 'Typhlops_brongersmianus', 'Laticauda_colubrina', 'Eublepharis_macularius', 'Thamnodynastes_hypoconia', 'Varanus_komodoensis', 'Anolis_meridionalis', 'Pseudoboa_nigra', 'Hemidactylus_mabouia', 'Xenodon_merremi', 'Phimophis_guerini', 'Tropidurus_oreadicus', 'Apostolepis_polylepis', 'Shinisaurus_crocodilurus', 'Ophiophagus_hannah', 'Colobosaura_modesta', 'Sphenodon_punctatus1', 'Apostolepis_cearensis', 'Corallus_hortulanus', 'Hydrophis_melanocephalus', 'Liotyphlops_ternetzii', 'Micrurus_brasiliensis', 'Gallus_gallus', 'Protobothrops_mucrosquamatus', 'Micrablepharus_maximiliani', 'Oxyrhopus_petolarius', 'Lacerta_viridis', 'Paroedura_picta', 'Kentropyx_calcarata', 'Deinagkistrodon_actutus', 'Python_bivittatus', 'Dopasia_gracilis', 'Ameivula_mumbuca', 'Tantilla_melanocephala', 'Anolis_brasiliensis', 'Chrysemys_picta', 'Bothrops_moojeni', 'Gekko_japonicus'}
	reqClades=['Iguania','Anguimorpha','Snake','Outgroup']
	setup(mainDir,allTaxa,bbFiles,outFile,reqClades)

if __name__ == "__main__":
    main()
