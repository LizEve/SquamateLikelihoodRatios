#!/bin/bash
import os
import glob
import re
from Bio import SeqIO

original_file = "SqCL_probes_v1.fasta"

with open(original_file, 'r') as original:
    # start list of all loci
    all_loci=[]
    for record in SeqIO.parse(original_file, 'fasta'):
        n = re.split("_|-",record.id)
        # n[1] gives name for AHE and genes, n[2] works for uces
        if n[1] != "uce":
            all_loci.append(n[1])
        elif n[1] == "uce":
            all_loci.append(n[2])
        else:
            print(str(n[1])+": does not conform to record id patterns accounted for in this script")    
    loci_set = set(all_loci)
print(len(all_loci),len(loci_set))

    
original_file = "SqCL_probes_v1.fasta"
corrected_file = "probes_sqc.fasta"

# This takes a while 20 mins or so I think. 
with open(original_file, 'r') as original, open(corrected_file, 'w') as corrected:
    l_counter = 1
    for locus in loci_set:
        p_counter = 1
        for record in SeqIO.parse(original_file, 'fasta'):
            n = re.split("_|-",record.id)
            # n[1] gives name for AHE and genes, n[2] works for uces
            if n[1] != "uce":
                if n[1] == locus:
                    record.description = "|source: " + record.id
                    # edit new id
                    newid = "sqc-"+str(l_counter)+"_p"+str(p_counter)
                    print(record.id,newid)
                    # replace new id
                    record.id = newid
                    SeqIO.write(record, corrected, 'fasta')
                    p_counter += 1
            elif n[1] == "uce":
                if n[2] == locus:
                    record.description = "|source: " + record.id
                    # edit new id
                    newid = "sqc-"+str(l_counter)+"_p"+str(p_counter)
                    print(record.id,newid)
                    # replace new id
                    record.id = newid
                    SeqIO.write(record, corrected, 'fasta')
                    p_counter += 1
            else:
                print(str(n[1])+": does not conform to record id patterns accounted for in this script")    
        l_counter +=1
print(loci_set)
print(5)
