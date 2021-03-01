#!/bin/bash
import glob
import os
import zipfile

def grabZip(zipFile,outFolder,suffix):
	zipZ = zipfile.ZipFile(zipFile, 'r')
	for file in zipZ.namelist():
		if file.endswith(suffix):
			zipZ.extract(file, outFolder)

def grabFiles(zipDir,outFolder,globber):
	os.chdir(zipDir)
	for f in glob.glob(globber):
		#grabZip(f,outFolder,'posterior.trees')
		#grabZip(f,outFolder,'posterior.var')
		#grabZip(f,outFolder,'posterior.log')
		grabZip(f,outFolder,'.nex')
		grabZip(f,outFolder,'1.log')
		grabZip(f,outFolder,'2.log')
		grabZip(f,outFolder,'3.log')
		grabZip(f,outFolder,'4.log')
		grabZip(f,outFolder,'1.trees')
		grabZip(f,outFolder,'2.trees')
		grabZip(f,outFolder,'3.trees')
		grabZip(f,outFolder,'4.trees')
		grabZip(f,outFolder,'1.var')
		grabZip(f,outFolder,'2.var')
		grabZip(f,outFolder,'3.var')
		grabZip(f,outFolder,'4.var')



# Need to edit to run in current directory and have specific zips files as user intput and output as automatic from zip folder name

zipDir='/Users/ChatNoir/Projects/Squam/Chapter_Squamate'
outFolder='constraintsByFamSmall_500k_all'
globber='constraintsByFamSmall_500k_prtl.zip'
grabFiles(zipDir,outFolder,globber)