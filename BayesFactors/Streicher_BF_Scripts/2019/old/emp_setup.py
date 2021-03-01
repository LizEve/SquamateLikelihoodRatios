#! /usr/bin/env python
from Bio import AlignIO
from Bio import SeqIO
from Bio.Alphabet import generic_dna
import os
import glob
import shutil 
import itertools
import re 

# Streicher version 
# requires a main folder with: phylip files, and template bayesblock 

# Make folders for each locus and put nexus file into folder
def makeFolders(n, suffix, mainDir):
	nexus=str(n)
	nexusIndex = nexus.find(suffix)
	fName = nexus[:nexusIndex]
	dirPath = os.path.join(mainDir,fName)
	
	# Create paths 
	dirPath = os.path.join(mainDir,fName)

	# Make directories
	if not os.path.exists(dirPath):
		os.mkdir(dirPath)	

	# Copy nexus to data folder
	print(n)
	SeqIO.convert(n,"phylip-relaxed",fName+".nex","nexus",alphabet=generic_dna)
	os.system("mv %s %s" % (fName+".nex",dirPath))
	return dirPath,fName


def editFile(fName,missingTaxa):
	# Open template bb file and new locus specific bb file 
	bbin=open('tacocat.bb')
	bbout=open(fName+'.bb',"w+")
	# For each line in template, remove missing taxa, replace tacocat with locus
	for line in bbin:
		for t in missingTaxa:
			line=line.replace(" "+t,"")
			line=line.replace("tacocat",fName)
			line=line.replace(",,",",")
		if "=;" in line:
			emptyConstraint=line.split(" ")[-3]
			missingTaxa.append(emptyConstraint)
		else:
			bbout.write(line)
	bbin.close()
	bbout.close()


def setup(mainDir,suffix,allTaxa,folders = True):
	os.chdir(mainDir)
	# For each nexus file, make a folder, and get a path to folder
	for phy in glob.glob('*%s' % suffix):
		if folders == True:
			# Make a folder and move MSA file into folder 
			dirPath,fName = makeFolders(phy, suffix, mainDir)
		else:
			nexus=str(phy)
			nexusIndex = nexus.find(suffix)
			fName = nexus[:nexusIndex]
			dirPath = os.path.join(mainDir,fName)
		n=fName+'.nex'
		# Move into folder and grab bb file 
		os.chdir(dirPath)
		os.system("cp ../tacocat.bb .")
		# Open MSA file and get list of missing taxa for locus
		alignment = AlignIO.read(open(n),"nexus")
		locusTaxa=[]
		for record in alignment:
			locusTaxa.append(str(record.id))
		# Get list of all taxa not in locus 
		missingTaxa=list(allTaxa.difference(set(locusTaxa)))
		editFile(fName,missingTaxa)
		os.system("rm tacocat.bb")
		#os.system("rm *.phy")
		os.chdir(mainDir)

#mainDir=os.getcwd()
suffix='.nexus.phy'
mainDir='/Users/ChatNoir/Projects/Squam/Streicher/ASTRAL_37_taxa_4178_alignments'
allTaxa={'homo_sapiens', 'anolis_carolinensis_2', 'bipes_sp', 'chrysemys_picta', 'aspidoscelis_tigris', 'cricosaura_sp', 'tupinambus_sp', 'plestiodon_fasciatus', 'hydrosaurus_sp', 'lacerta_sp', 'lepidophyma_sp', 'typhlops_jamaicensis', 'cordylosaurus_sp', 'lanthanotus_borneensis', 'anelytropsis_sp', 'strophurus_sp', 'gonatodes_sp', 'cordylus_sp', 'tiliqua_sp', 'rhineura_sp', 'saltorius_sp', 'anniella_pulchra', 'dibamus_sp', 'heloderma_suspectum', 'gekko_sp', 'alligator_mississippiensis', 'acontias_sp', 'lialis_sp', 'varanus_exanthematicus', 'micrurus_fulvius', 'pholidobollus_sp', 'coleonyx_variegatus', 'gallus_gallus', 'xenosaurus_platyceps', 'sphenodon_sp', 'uta_stansburiana', 'python_molurus_2'}

setup(mainDir,suffix,allTaxa)