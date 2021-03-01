#! /usr/bin/env python
import glob
import os
import time

def align_mafft(files):
	for f in files:
		print(f)
		os.system("mafft --maxiterate 1000 --globalpair --adjustdirection %s > %s.aln" % (f,f))

start = time.time()

files = glob.glob("*.fasta")
align_mafft(files)

end = time.time()
print(end-start)
print((end-start)/60)
print((end-start)/(60*60))