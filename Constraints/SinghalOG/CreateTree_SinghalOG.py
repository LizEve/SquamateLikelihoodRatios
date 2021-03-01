import re
import os
import sys
from tempfile import mkstemp
from shutil import copyfile
import dendropy 
import ete3
from collections import OrderedDict
from operator import itemgetter   
from Bio import AlignIO
import glob
import pickle


def getStarTree(allTaxa):
    """Takes
    Returns"""
    # Create empty tree
    t1 = ete3.Tree()
    
    # Iteratively add all taxa
    for tx in allTaxa:
        t1.add_child(name=tx) 
        
    return t1

def addBipart(t,allTaxa,newBipart):
    """ Takes a tree, a bipartition, and a list of all taxa in the tree
    Returns new tree
    Has bug - say if AB is new bipart, but is already in clade ABCD, AB will be added as new bipartition, and the ABCD will not be preserved. 
    Solution - add biparts from fewest to most taxa. 
    This failed for singhal dataset, fixing those trees by hand. 
    """
    
    # Get list of all leaf names minus bipartition    
    save = list(set(allTaxa)-set(newBipart))

    #print(len(allTaxa),len(newBipart),len(save))

    # get mrca of taxa in bipart
    mrca=t.get_common_ancestor(newBipart)

    # Separate and copy bipart as subtree
    subtree=mrca.copy()

    # Prune off any taxa not in bipartition from subtree, in case of polytomy
    subtree.prune(newBipart)

    # save all of tree except bipartition
    t.prune(save)

    # Add new bipartition subtree to original tree as new node.  
    t.add_child(subtree)

    return(t)

def makeTree(fileName,t,outTree,allTaxa):
    # Create list with biparts added 
    addedBipartsT1 = []
    addedBipartsNeither = []
    lenDict = {}

    # Open file 
    for line in open(fileName):
        line = line.strip("\n")
        #print(line)
        # Make into list 
        bipart = line.split(" ")
        # Clean up spaces in names 
        bipart = [x.replace(' ','') for x in bipart]
        # remove empty items in list 
        bipart = list(filter(None, bipart))
        #print(bipart)
        # Write bipart and number of taxa in bipart 
        lenDict[str(bipart)]=[int(len(bipart)),bipart]

        # Order biparts by number of taxa
        # Have to do this because of bipart adding fuction, more details in notes for that function. 
        byTaxaDict = OrderedDict(sorted(lenDict.items(), key = itemgetter(1), reverse = False))

        # For each bipartition in dictionary, add to tree
        for key,value in byTaxaDict.items():
            #print(value[0])
            bipart = value[1]
            #print(t.get_ascii())
            newBipart = bipart
            # Get list of all leaf names minus bipartition    
            save = list(set(allTaxa)-set(newBipart))
            # get mrca of taxa in bipart
            mrca=t.get_common_ancestor(newBipart)
            # Separate and copy bipart as subtree
            subtree=mrca.copy()
            #print(subtree.get_ascii())
            # Prune off any taxa not in bipartition from subtree, in case of polytomy
            subtree.prune(newBipart)
            # save all of tree except bipartition
            t.prune(save)
            #print(t.get_ascii())
            # Add new bipartition subtree to original tree as new node.  
            t.add_child(subtree)
    print(t.get_ascii())
    # Convert to dendropy tree list to write to file
    # Ete3 newick formatting: 5 = internal and leaf branches + leaf names, 9 = leaf names
    dt = dendropy.Tree.get(data=t.write(format = 9),schema = 'newick')

    # Output tree
    out=open(outTree, "w")
    dt.write(file=out, schema='newick')
    out.close()


def main():
    dataSet = "SinghalOG"
    # Get list of all taxa 
    allTaxa={'Ameiva_ameiva', 'Ameivula_mumbuca', 'Amphisbaena_alba', 'Anolis_brasiliensis', 'Anolis_meridionalis', 'Apostolepis_cearensis', 'Apostolepis_polylepis', 'Bothrops_lutzi', 'Bothrops_moojeni', 'Bothrops_pauloensis', 'Brasiliscincus_heathi', 'Chironius_exoletus', 'Colobosaura_modesta', 'Copeoglossum_nigropunctatum', 'Corallus_hortulanus', 'Erythrolamprus_almadensis', 'Erythrolamprus_poecilogyrus', 'Erythrolamprus_reginae', 'Gallus_gallus', 'Gymnodactylus_amarali', 'Hemidactylus_mabouia', 'Imantodes_cenchoa', 'Kentropyx_calcarata', 'Leptodeira_annulata', 'Liotyphlops_ternetzii', 'Lygophis_paucidens', 'Micrablepharus_maximiliani', 'Micrurus_brasiliensis', 'Notomabuya_frenata', 'Oxyrhopus_petolarius', 'Oxyrhopus_trigeminus', 'Philodryas_nattereri', 'Philodryas_olfersii', 'Phimophis_guerini', 'Pseudoboa_neuwiedii', 'Pseudoboa_nigra', 'Psomophis_joberti', 'Sibynomorphus_mikanii', 'Taeniophallus_occipitalis', 'Tantilla_melanocephala', 'Thamnodynastes_hypoconia', 'Trilepida_brasiliensis', 'Tropidurus_oreadicus', 'Typhlops_brongersmianus', 'Xenodon_merremi'}
    # define folders and files 
    mainFolder = os.path.join("/Users/ChatNoir/Projects/Squam/scripts/Constraints",dataSet)
    scler = os.path.join(mainFolder,dataSet+"_Sclero_iq.txt")
    toxpoly = os.path.join(mainFolder,dataSet+"_ToxPoly_iq.txt")
    # out file names 
    sOut = os.path.join(mainFolder,"scleroglossa.constraint")
    tpOut = os.path.join(mainFolder,"toxicoferaP.constraint")
    inFiles=(scler,toxpoly)
    outFiles=(sOut,tpOut)
    # Get initial star tree
    t = getStarTree(allTaxa)
    for x in range(0,2):
        makeTree(inFiles[x],t,outFiles[x],allTaxa)

if __name__=='__main__':
    main()