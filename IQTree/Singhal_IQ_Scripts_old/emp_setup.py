#! /usr/bin/env python
from Bio import AlignIO
from Bio import SeqIO
from Bio.Alphabet import generic_dna
import os
import glob
import shutil 
import itertools
import re 
import dendropy
import copy 

# Streicher version IQ Tree
# requires a main folder with: nexus files, and 5 template constraint files



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



def editFile(n,fName,c,locusList):
	outCon = fName+"."+c
	# Read in tree file 
	allTaxa = dendropy.TaxonNamespace()
	constraint = dendropy.Tree.get(path=c, schema="newick", taxon_namespace=allTaxa)
	# Keep taxa from locus 
	constraint.retain_taxa_with_labels(locusList)
	# WRite new locus specific constraint file 
	constraint.write(path=outCon, schema="newick")
	# Remove template constraint file 
	os.system("rm %s" % c)


def setup(mainDir,conFiles,folders = True):

	os.chdir(mainDir)
	
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

		# Move into folder and grab constraint files
		os.chdir(dirPath)
		os.system("cp ../*.constraint .")

		# Open locus sequence alignment, get taxa names
		locusTaxa = dendropy.TaxonNamespace()
		alignment = dendropy.DnaCharacterMatrix.get(path=nex, schema="nexus", taxon_namespace=locusTaxa)
		locusList = locusTaxa.labels()

		for f in conFiles:
			editFile(nex,fName,f,locusList)
		# When done with locus 
		os.chdir(mainDir)

def main():
	dataSet = "Singhal"
	conFiles=["scleroglossa.constraint", "toxicoferaP.constraint", "toxicoferaAI.constraint", "toxicoferaSA.constraint", "toxicoferaSI.constraint"]
	mainDir = os.path.join("/Users/ChatNoir/Projects/Squam/",dataSet,dataSet+"IQ")

	setup(mainDir,conFiles)


if __name__ == "__main__":
    main()