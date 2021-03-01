#! /usr/bin/env python

import os
import glob
import shutil 
import itertools
import re 

# requires a main folder with:
# nexusFiles - folder containing all nexus files to run

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
	os.system("cp %s %s" % (n,dirPath))
	return dirPath,fName



def setup(mainDir,suffix,folders = False):
	os.chdir(mainDir)
	for n in glob.glob('*%s' % suffix):
		if folders == True:
			x = makeFolders(n, suffix, mainDir)
			dirPath = x[0]
			fName = x[1] 
		else:
			nexus=str(n)
			nexusIndex = nexus.find(suffix)
			fName = nexus[:nexusIndex]
			dirPath = os.path.join(mainDir,fName)
		os.chdir(mainDir)



mainDir=os.getcwd()
setup(mainDir,".ntg.nex",True)

