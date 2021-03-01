#! /usr/bin/env python
import glob
import os
import time
# Runs gblocks with default parameters


def trim_gblocks(files,gblocks_path):
	for f in files:
		print(f)
		os.system("%s %s -t=d -b1=14 -b2=14 -b3=8 -b4=5" % (gblocks_path,f))

start = time.time()

files = glob.glob("*.aln")
gblocks_path = "/home/gmount/Gblocks_0.91b/Gblocks"
trim_gblocks(files, gblocks_path)


end = time.time()
print(end-start)
print((end-start)/60)
print((end-start)/(60*60))