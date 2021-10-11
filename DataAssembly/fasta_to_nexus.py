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
	short = []
	for f in files:
		new_file = f.split(".")[0]+".aln.temp"
		newrecords = []
		for record in SeqIO.parse(open(f),'fasta'):
			if len(record) >= 100:
				rev = record.id.split('_R_') # mafft reversed and renamed some seqs
				if len(rev) == 1:
					newHeader = record.id.split("_")[0]+"_"+record.id.split("_")[1]
				elif len(rev) == 2:
					newHeader = record.id.split("_")[2]+"_"+record.id.split("_")[3]
				if 'Ophiosaurus' not in newHeader:
					record.id=newHeader
					newrecords.append(record)
				else:
					pass
			else:
				short.append(f)
		SeqIO.write(newrecords, new_file, 'fasta')
	print(str(set(short)))


def convert_to_nexus(files):
	for t in files:
		nexus_file = t.split(".")[0]+".nex"
		try:
				AlignIO.convert(t, "fasta", nexus_file, "nexus", alphabet=Alphabet.generic_dna)
		except ValueError:
				os.system("rm %s" % nexus_file)
				print("%s no sequence left after trimming" % t)


start = time.time()

files = glob.glob("*.fasta.aln-gb")
remove_locus_name(files)
temp_files = glob.glob("*.aln.temp")
convert_to_nexus(temp_files)

end = time.time()
print(end-start)
print((end-start)/60)
print((end-start)/(60*60))