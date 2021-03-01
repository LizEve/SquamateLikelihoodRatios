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

def setup(mainDir,suffix,allTaxa,bbFile,outFile):
	os.chdir(mainDir)
	nexMissing={}
	# For each nexus file
	for nex in glob.glob('*%s' % suffix):
		removedConstraints=[]
		# Open MSA file and get list of missing taxa for locus
		alignment = AlignIO.read(open(nex),"nexus")
		locusTaxa=[]
		for record in alignment:
			locusTaxa.append(str(record.id))
		# Get list of all taxa not in locus 
		missingTaxa=list(allTaxa.difference(set(locusTaxa)))
		bbin=open(bbFile)
		# For each line in template, remove missing taxa
		for line in bbin:
			# Remove all missing taxa, replace with nothing, then remove double commas
			for t in missingTaxa:
				line=line.replace(" "+t,"",20)
				line=line.replace(",,",",",20)
			# Flag all empty constraints. and store in list 
			if "=;" in line:
				# Grab name of constraint to remove from prset topology constraints 
				emptyConstraint=line.split(" ")[-3]
				removedConstraints.append(emptyConstraint)
		bbin.close()
		# save in dictionary
		nexMissing[nex]=removedConstraints
		#print(nex,removedConstraints)
	# print dictionary to file 
	print(nexMissing)
	pd.DataFrame.from_dict(data=nexMissing, orient='index').to_csv(outFile, header=False)

def main():
	suffix='.nex'
	mainDir='/Users/ChatNoir/Projects/Squam/Streicher/StreicherNexus'
	bbFile = "/Users/ChatNoir/Projects/Squam/scripts/BayesFactors/Streicher_BF_Scripts/tacocat_allconstraint.bb"
	outFile = '/Users/ChatNoir/Projects/Squam/Streicher/missingClades.csv'
	# set of all taxa 
	Taxa={'chrysemys_picta', 'cricosaura_sp', 'varanus_exanthematicus', 'pholidobollus_sp', 'dibamus_sp', 'python_molurus_2', 'tupinambus_sp', 'anniella_pulchra', 'lialis_sp', 'hydrosaurus_sp', 'lepidophyma_sp', 'xenosaurus_platyceps', 'gekko_sp', 'bipes_sp', 'typhlops_jamaicensis', 'cordylosaurus_sp', 'cordylus_sp', 'homo_sapiens', 'uta_stansburiana', 'sphenodon_sp', 'acontias_sp', 'tiliqua_sp', 'anolis_carolinensis_2', 'coleonyx_variegatus', 'micrurus_fulvius', 'saltorius_sp', 'plestiodon_fasciatus', 'gonatodes_sp', 'lanthanotus_borneensis', 'anelytropsis_sp', 'heloderma_suspectum', 'lacerta_sp', 'rhineura_sp', 'aspidoscelis_tigris', 'gallus_gallus', 'strophurus_sp', 'alligator_mississippiensis'}
	setup(mainDir,suffix,Taxa,bbFile,outFile)

if __name__ == "__main__":
    main()
