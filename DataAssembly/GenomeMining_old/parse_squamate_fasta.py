#! /usr/bin/env python
import glob
from Bio import SeqIO
import os
import time


def convertLoci(files) :
	loci = []
	# store all loci names in list. this will include duplicates
	for f in files :
		sequences = list(SeqIO.parse(open(f, 'r'), 'fasta'))
		for s in sequences :
			loci.append(s.id)
	mainDir = os.getcwd()
	if not os.path.exists(mainDir+"/loci_files/"):
		os.mkdir(mainDir+"/loci_files/")	
	count = 1
	for locus in set(loci) :
		print "Processing locus %s / %s." %(str(count), len(set(loci)))
		outName = mainDir+"/loci_files/" + locus + ".fasta"
		new_sequences = []
		for f in files :
			taxon = f.split(".fasta")[0]
			sequences = list(SeqIO.parse(open(f, 'r'), 'fasta'))
			for s in sequences :
				if locus == s.id :
					new_seq = s
					new_seq.id = taxon.split("/")[-1] + "_" + s.id.split(" ")[-1]
					new_seq.name = ""
					new_seq.description = ""
					new_sequences.append(new_seq)
					
		SeqIO.write(new_sequences, outName, "fasta")			
		count += 1


start = time.time()

files = glob.glob("*.fasta")
## convert from taxon-based fasta to loci-based fasta for alignments
convertLoci(files)

end = time.time()
print(end-start)
print((end-start)/60)
print((end-start)/(60*60))