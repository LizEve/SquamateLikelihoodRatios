
import re
import os
from Bio import SeqIO


def getProbeDict(probe_file):
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
    return sqcProbeDict


def renameFasta(original_file,corrected_file,ident,sqcProbeDict):
    with open(corrected_file, 'w') as corrected:
            for record in SeqIO.parse(original_file, 'fasta'):
                # save original record id before altering it 
                oldHeader = record.id
                n = re.split("\|",oldHeader)
                sqcName = re.split(":",n[3])[1]
                contig = str(ident)+"_contig"+re.split(":",n[1])[1]
                record.id = sqcProbeDict[sqcName]
                record.description = contig
                SeqIO.write(record, corrected, 'fasta')

probe_file='/home/gmount/Squamates/GrabGenomeSeqs/uce-genome/probes_sqc.fasta'
sqcProbeDict=getProbeDict(probe_file)
mainDir = "/home/gmount/Squamates/GrabGenomeSeqs/rename-fastas"
os.chdir(mainDir)
renameFasta("allmis2.fasta","Alligator_mississippiensis.fasta","Alligator_mississippiensis",sqcProbeDict)
renameFasta("anocar2.fasta","Anolis_carolinensis.fasta","Anolis_carolinensis",sqcProbeDict)
renameFasta("chrpic3.fasta","Chrysemys_picta.fasta","Chrysemys_picta",sqcProbeDict)
renameFasta("crohor1.fasta","Crotalus_horridus.fasta","Crotalus_horridus",sqcProbeDict)
renameFasta("deiacu1.fasta","Deinagkistrodon_actutus.fasta","Deinagkistrodon_actutus",sqcProbeDict)
renameFasta("eubmac1.fasta","Eublepharis_macularius.fasta","Eublepharis_macularius",sqcProbeDict)
renameFasta("gekjap1.fasta","Gekko_japonicus.fasta","Gekko_japonicus",sqcProbeDict)
renameFasta("gopevg1.fasta","Gopherus_evgoodei.fasta","Gopherus_evgoodei",sqcProbeDict)
renameFasta("hydmel1.fasta","Hydrophis_melanocephalus.fasta","Hydrophis_melanocephalus",sqcProbeDict)
renameFasta("lacbil1.fasta","Lacerta_bilineata.fasta","Lacerta_bilineata",sqcProbeDict)
renameFasta("lacvir1.fasta","Lacerta_viridis.fasta","Lacerta_viridis",sqcProbeDict)
renameFasta("latcor1.fasta","Laticauda_colubrina.fasta","Laticauda_colubrina",sqcProbeDict)
renameFasta("ophhan1.fasta","Ophiophagus_hannah.fasta","Ophiophagus_hannah",sqcProbeDict)
renameFasta("ophgra1.fasta","Ophisaurus gracilis.fasta","Ophisaurus gracilis",sqcProbeDict)
renameFasta("parpic1.fasta","Paroedura_picta.fasta","Paroedura_picta",sqcProbeDict)
renameFasta("podmur1.fasta","Podarcis_muralis.fasta","Podarcis_muralis",sqcProbeDict)
renameFasta("pogvit1.fasta","Pogona_vitticeps.fasta","Pogona_vitticeps",sqcProbeDict)
renameFasta("promuc1.fasta","Protobothrops_mucrosquamatus.fasta","Protobothrops_mucrosquamatus",sqcProbeDict)
renameFasta("pytbiv5.fasta","Python_bivittatus.fasta","Python_bivittatus",sqcProbeDict)
renameFasta("salmer1.fasta","Salvator_meriana.fasta","Salvator_meriana",sqcProbeDict)
renameFasta("shicro1.fasta","Shinisaurus_crocodilurus.fasta","Shinisaurus_crocodilurus",sqcProbeDict)
renameFasta("sphpun1.fasta","Sphenodon_punctatus1.fasta","Sphenodon_punctatus1",sqcProbeDict)
renameFasta("sphpun2.fasta","Sphenodon_punctatus2.fasta","Sphenodon_punctatus2",sqcProbeDict)
renameFasta("varkom1.fasta","Varanus_komodoensis.fasta","Varanus_komodoensis",sqcProbeDict)
renameFasta("zooviv1.fasta","Zootoca_vivipara.fasta","Zootoca_vivipara",sqcProbeDict)


