#! /usr/bin/env python
import glob
import os
from Bio import AlignIO
from Bio import Alphabet
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment



def removeTaxa(infolder,outfolder,required):
        os.chdir(infolder)
        # Make a folder
        if not os.path.exists(outfolder):
                os.mkdir(outfolder)
        o = open(os.path.join(outfolder,"outgroup.info"),"w")
        q = open(os.path.join(outfolder,"excludedLoci.info"),"w")
        # Iterate nexus files
        for f in glob.glob("*.nex"):
                # List for new nexus file, and file path
                newrecords = []
                all_taxa = []
                new_file =os.path.join(outfolder,f)
                # Open alignment nexus
                alignment = AlignIO.read(open(f),'nexus')
                # Make list of all taxa that are required in your alignment
                for record in alignment:
                        all_taxa.append(record.id)
                        if record.id in required:
                                newrecords.append(record)
                # For all alignments that have all taxa, continue
                if len(newrecords) == len(required):
                # Add option taxa to new alignment
                        if "Gallus_gallus" in all_taxa:
                                newRecords=matchRename("Gallus_gallus","Outgroup_outgroup",alignment,newrecords)
                                nr = MultipleSeqAlignment(newrecords)
                                AlignIO.write(nr, new_file, 'nexus')
                                o.write("%s,%s\n" % (f,"Gallus_gallus"))
                        elif "Alligator_mississippiensis" in all_taxa:
                                newRecords=matchRename("Alligator_mississippiensis","Outgroup_outgroup",alignment,newrecords)
                                nr = MultipleSeqAlignment(newrecords)
                                AlignIO.write(nr, new_file, 'nexus')
                                o.write("%s,%s\n" % (f,"Alligator_mississippiensis"))
                        elif "Chrysemys_picta" in all_taxa:
                                newRecords=matchRename("Chrysemys_picta","Outgroup_outgroup",alignment,newrecords)
                                nr = MultipleSeqAlignment(newrecords)
                                AlignIO.write(nr, new_file, 'nexus')
                                o.write("%s,%s\n" % (f,"Chrysemys_picta"))
                        else: 
                                q.write("%s,%s\n" % (f,"noOutgroup"))
                elif len(newrecords) < len(required):
                        q.write("%s,%s\n" % (f,"tooFewTaxa"))

#required = ["Anolis_brasiliensis", "Anolis_meridionalis", "Tropidurus_oreadicus", "Hemidactylus_mabouia", "Gymnodactylus_amarali", "Colobosaura_modesta", "Micrablepharus_maximiliani", "Brasiliscincus_heathi", "Copeoglossum_nigropunctatum", "Notomabuya_frenata", "Ameivula_mumbuca", "Kentropyx_calcarata", "Amphisbaena_alba",  "Tantilla_melanocephala", "Bothrops_moojeni", "Typhlops_brongersmianus",  "Ophisaurus_gracilis"]
required = ["Anolis_brasiliensis", "Anolis_meridionalis", "Tropidurus_oreadicus", "Hemidactylus_mabouia", "Gymnodactylus_amarali", "Colobosaura_modesta", "Micrablepharus_maximiliani", "Brasiliscincus_heathi", "Copeoglossum_nigropunctatum", "Notomabuya_frenata", "Ameivula_mumbuca", "Kentropyx_calcarata", "Amphisbaena_alba",  "Tantilla_melanocephala", "Bothrops_moojeni", "Typhlops_brongersmianus"]
infolder = os.getcwd()
outfolder = os.path.join(infolder,"nexus_files_subsampled")
removeTaxa(infolder,outfolder,required)