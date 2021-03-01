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

# Burbrink version IQ Tree
# requires a main folder with: nexus files, and 5 template constraint files


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


def setup(mainDir,suffix,folders = True):
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
		# Move into folder and grab constraint files
		os.chdir(dirPath)
		os.system("cp ../*.constraint .")
		# Open locus sequence alignment, get taxa names
		locusTaxa = dendropy.TaxonNamespace()
		alignment = dendropy.DnaCharacterMatrix.get(path=nex, schema="nexus", taxon_namespace=locusTaxa)
		locusList = locusTaxa.labels()
		editFile(nex,fName,"scleroglossa.constraint",locusList)
		editFile(nex,fName,"toxicoferaP.constraint",locusList)
		editFile(nex,fName,"toxicoferaAI.constraint",locusList)
		editFile(nex,fName,"toxicoferaSA.constraint",locusList)
		editFile(nex,fName,"toxicoferaSI.constraint",locusList)
		os.chdir(mainDir)

suffix='.nex'
mainDir='/Users/ChatNoir/Projects/Squam/Burbrink/BurbrinkData'

setup(mainDir,suffix)