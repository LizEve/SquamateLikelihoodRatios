#! /usr/bin/env python
from Bio import SeqIO
import os
import glob
import sys

'''
Copy nexus files from one folder and paste them to a different folder as fasta files. 
'''

def fastaToNexus(nexusDir,HoTdir):
    for n in glob.glob('*.nex'):
        # Get locus name 
        fName=os.path.splitext(n)[0]
        # Read in nexus file 
        records = SeqIO.parse(n, "nexus")
        # Convert to fasta 
        count = SeqIO.write(records, os.path.join(HoTdir,fName+".fasta"), "fasta")
        print("Converted %i records" % count)


def main():
    
    dataSet = sys.argv[1]#"SinghalOG"
    
    # Create path variables 
    nexusDir=os.path.join("/work/gmount",dataSet+'SS')
    HoTDir=os.path.join("/work/gmount",dataSet+'_HoT')

    # Make directory if needed 
    if not os.path.exists(HoTDir):
        os.mkdir(HoTDir)

    # Move into dir with nexus files	
    os.chdir(nexusDir)
    fastaToNexus(nexusDir,HoTDir)

if __name__ == "__main__":
    main()