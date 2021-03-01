#! /usr/bin/env python
from Bio.Nexus import Nexus
from os import walk
import glob
import os
import re

mainDir = '/home/gmount/Sq_CL/ConcatTreeThreeGen'
outFile = '/home/gmount/Sq_CL/ConcatTreeThreeGen/ThreeGenConcat2.nex'
globber= "*threeg.nex"

file_list = []
os.chdir(mainDir)
for file in glob.glob(globber):
    with open(file) as fp:
        for i,line in enumerate(fp):
            if i==2:
                nchar=line.split('=')[2].split(";")[0]
                if int(nchar) >= 100:
                     file_list.append(file)
            elif i > 2:
                break
    fp.close()


nexi =  [(fname, Nexus.Nexus(fname)) for fname in file_list]
combined = Nexus.combine(nexi)
combined.write_nexus_data(filename=open(outFile, 'w'))