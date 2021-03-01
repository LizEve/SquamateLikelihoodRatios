#! /usr/bin/env python
import glob
import os
import sys
from Bio import AlignIO
from Bio import Alphabet
from Bio import SeqIO
from Bio.Nexus import Nexus
import time

def remove_locus_name(files):
	for f in files:
		new_file = f.split(".")[0]+".aln.temp"
		newrecords = []
		for record in SeqIO.parse(open(f),'fasta'):
			newHeader = record.id.split("_")[0]+"_"+record.id.split("_")[1]
			record.id=newHeader
			newrecords.append(record)
		SeqIO.write(newrecords, new_file, 'fasta')


def convert_to_nexus(files):
	for t in files:
		nexus_file = t.split(".")[0]+".sixg.nex"
		try:
				AlignIO.convert(t, "fasta", nexus_file, "nexus", alphabet=Alphabet.generic_dna)
		except ValueError:
				os.system("rm %s" % nexus_file)
				print("%s no sequence left after trimming" % t)
				sys.exc_clear()

start = time.time()

files = glob.glob("*.fasta.aln-gb")
remove_locus_name(files)
temp_files = glob.glob("*.aln.temp")
convert_to_nexus(temp_files)

end = time.time()
print(end-start)
print((end-start)/60)
print((end-start)/(60*60))