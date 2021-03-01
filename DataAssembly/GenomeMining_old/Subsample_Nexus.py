#! /usr/bin/env python
import glob
import os
from Bio import AlignIO
from Bio import Alphabet
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment

def removeTaxa(infolder,outfolder,required,optional):
        os.chdir(infolder)
        # Make a folder
        if not os.path.exists(outfolder):
                os.mkdir(outfolder)
        # Iterate nexus files
        for f in glob.glob("*.nex"):
                # List for new nexus file, and file path
                newrecords = []
                new_file =os.path.join(outfolder,f)
                # Open alignment nexus
                alignment = AlignIO.read(open(f),'nexus')
                # Make list of all taxa that are required in your alignment
                for record in alignment:
                        if record.id in required:
                                z=record.seq.count("t")+record.seq.count("g")+record.seq.count("c")+record.seq.count("a")
                                if z > 10:
                                        newrecords.append(record)
                # For all alignments that have all taxa, continue
                if len(newrecords) == len(required):
                # Add option taxa to new alignment
                        for record in alignment:
                                if record.id in optional:
                                        newrecords.append(record)
                        # change format and write out to new file
                        nr = MultipleSeqAlignment(newrecords)
                        AlignIO.write(nr, new_file, 'nexus')
                        # Print all files that only have the required amount of taxa
                        if len(newrecords) < len(required):
                                print(f,len(newrecords))
                #elif len(newrecords) < len(required):
			#print("%s does not have all the taxa" % f)

required = ["Anolis_brasiliensis", "Anolis_meridionalis", "Tropidurus_oreadicus", "Hemidactylus_mabouia", "Gymnodactylus_amarali", "Colobosaura_modesta", "Micrablepharus_maximiliani", "Brasiliscincus_heathi", "Copeoglossum_nigropunctatum", "Notomabuya_frenata", "Ameivula_mumbuca", "Kentropyx_calcarata", "Amphisbaena_alba",  "Tantilla_melanocephala", "Bothrops_moojeni", "Typhlops_brongersmianus",  "Ophisaurus_gracilis"]
optional = ["Chrysemys_picta","Alligator_mississippiensis","Gallus_gallus"]
infolder = "/home/gmount/Sq_CL/NTG/nexus_files"
outfolder = "/home/gmount/Sq_CL/NTG/nexus_files_reducedtx"
removeTaxa(infolder,outfolder,required,optional)