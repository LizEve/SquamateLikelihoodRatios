
import re
import os
from Bio import SeqIO

def renameFasta(original_file,corrected_file,ident,probe_file):
    sqcProbeDict={}
    sqcList=[]
    # Create dictionary of temporary sqc name to original probe name. 
    for record in SeqIO.parse(probe_file, 'fasta'):
        d=re.split("\|",record.description)
        ##['sqc-1_p4 ', 'source: Anolis-uce-1869-p1-UCE_60']
        sqc = re.split("_",d[0])[0]
        sqcList.append(sqc)
        longName = re.split("_|-",d[1])
        ##['source: Gekko', 'uce', '3922', 'p1', 'UCE', '73']
        ##['source: anoCar2', 'L118', 'AHE', '0']
        ##['source: Gekko', 'BACH1', 'genes', '126']
        if longName[2] == "AHE":
            newName=longName[2]+"-"+str(longName[1])   
        elif longName[2] == "genes":    
            newName="gene-"+str(longName[1]) 
        elif longName[1] == "uce":
            newName=longName[1]+"-"+str(longName[2])  
        sqcProbeDict[sqc]=newName
    print(len(sqcProbeDict))

    with open(corrected_file, 'w') as corrected:
        for key, value in sqcProbeDict.items():
            for record in SeqIO.parse(original_file, 'fasta'):
                # save original record id before altering it 
                oldHeader = record.id
                n = re.split("\|",oldHeader)
                sqcName = re.split(":",n[3])[1]
                contig = str(ident)+"_contig"+re.split(":",n[1])[1]
                print(sqcName,key)
                if sqcName == key:
                    record.id = value
                    record.description = contig
                    print(sqcName)
                    print(record.id)
                    SeqIO.write(record, corrected, 'fasta')

probe_file = "/home/gmount/Sq_CL/GrabGenomeSeqs/uce-genome/probes_sqc.fasta"
mainDir = "/home/gmount/Sq_CL/GrabGenomeSeqs/uce-genome/sixGenomes-genome-fasta"
os.chdir(mainDir)
renameFasta("allmis2.fasta","Alligator_mississippiensis.fasta","ASM28112v4",probe_file)
renameFasta("anocar2.fasta","Anolis_carolinensis.fasta","AnoCar2.0",probe_file) 
renameFasta("chrpic3.fasta","Chrysemys_picta.fasta","Chrysemyspictabellii3.0.3",probe_file) 
renameFasta("gekjap1.fasta","Gekko_japonicus.fasta","GekkojaponicusV1.1",probe_file)
renameFasta("ophgra1.fasta","Ophisaurus_gracilis.fasta","100119",probe_file)
renameFasta("pytbiv5.fasta","Python_bivittatus.fasta","Pythonmolurusbivittatus5.0.2",probe_file)




