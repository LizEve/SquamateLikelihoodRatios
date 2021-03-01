#! /usr/bin/env python
import glob
import os
from Bio import AlignIO
from Bio import Alphabet
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment

'''
Select loci that contain all ingroupt taxa, choose one outgroup and rename to "Outgroup_outgroup". 
Counts taxa as missing if sequence does not have at least 5 of each bp.
Outputs info on why loci were excluded. 
`outgroup.info` - locus,outgroup
`excludedLoci.info` - locus,reason for exclusion
'''


def matchRename(oldName,newName,alignment,newrecordsList):
        for record in alignment:
                if record.id == oldName:
                        record.id = newName
                        newrecordsList.append(record)
        return newrecordsList

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
                                # Count instances of atgc, require each to occur 5 times
                                atgc=[record.seq.count("t"),record.seq.count("g"),record.seq.count("c"),record.seq.count("a")]
                                if all(i >= 5 for i in atgc):
                                        newrecords.append(record)
                                else:
                                        q.write("%s,%s\n" % (f,"tooShortSeq"))
                # For all alignments that have all taxa, continue
                if len(newrecords) == len(required):
                # Add option taxa to new alignment
                        if "Gallus_gallus" in all_taxa:
                                newRecords=matchRename("Gallus_gallus","Outgroup_outgroup",alignment,newrecords)
                                # change format and write out to new file
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
                        AlignIO.write(nr, new_file, 'nexus')
                # Not enough taxa 
                elif len(newrecords) < len(required):
                        q.write("%s,%s\n" % (f,"tooFewTaxa"))

#required = ["Anolis_brasiliensis", "Anolis_meridionalis", "Tropidurus_oreadicus", "Hemidactylus_mabouia", "Gymnodactylus_amarali", "Colobosaura_modesta", "Micrablepharus_maximiliani", "Brasiliscincus_heathi", "Copeoglossum_nigropunctatum", "Notomabuya_frenata", "Ameivula_mumbuca", "Kentropyx_calcarata", "Amphisbaena_alba",  "Tantilla_melanocephala", "Bothrops_moojeni", "Typhlops_brongersmianus",  "Ophisaurus_gracilis"]
required = ["Anolis_brasiliensis", "Anolis_meridionalis", "Tropidurus_oreadicus", "Hemidactylus_mabouia", "Gymnodactylus_amarali", "Colobosaura_modesta", "Micrablepharus_maximiliani", "Brasiliscincus_heathi", "Copeoglossum_nigropunctatum", "Notomabuya_frenata", "Ameivula_mumbuca", "Kentropyx_calcarata", "Amphisbaena_alba",  "Tantilla_melanocephala", "Bothrops_moojeni", "Typhlops_brongersmianus"]
#infolder = "/home/gmount/Sq_CL/NTG/nexus_files"
#outfolder = "/home/gmount/Sq_CL/NTG/nexus_files_reducedtx"
infolder = os.getcwd()
outfolder = os.path.join(infolder,"nexus_files_subsampled")
removeTaxa(infolder,outfolder,required)