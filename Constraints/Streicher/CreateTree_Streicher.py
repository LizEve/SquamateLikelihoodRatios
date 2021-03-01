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
    dataSet = "Streicher"
    # Get list of all taxa 
    allTaxa={'chrysemys_picta', 'cricosaura_sp', 'varanus_exanthematicus', 'pholidobollus_sp', 'dibamus_sp', 'python_molurus_2', 'tupinambus_sp', 'anniella_pulchra', 'lialis_sp', 'hydrosaurus_sp', 'lepidophyma_sp', 'xenosaurus_platyceps', 'gekko_sp', 'bipes_sp', 'typhlops_jamaicensis', 'cordylosaurus_sp', 'cordylus_sp', 'homo_sapiens', 'uta_stansburiana', 'sphenodon_sp', 'acontias_sp', 'tiliqua_sp', 'anolis_carolinensis_2', 'coleonyx_variegatus', 'micrurus_fulvius', 'saltorius_sp', 'plestiodon_fasciatus', 'gonatodes_sp', 'lanthanotus_borneensis', 'anelytropsis_sp', 'heloderma_suspectum', 'lacerta_sp', 'rhineura_sp', 'aspidoscelis_tigris', 'gallus_gallus', 'strophurus_sp', 'alligator_mississippiensis'}
    # define folders and files 
    mainFolder = os.path.join("/Users/ChatNoir/Projects/Squam/scripts/Constraints",dataSet)
    scler = os.path.join(mainFolder,dataSet+"_Sclero_iq.txt")
    toxpoly = os.path.join(mainFolder,dataSet+"_ToxPoly_iq.txt")
    toxai = os.path.join(mainFolder,dataSet+"_ToxAI_iq.txt")
    toxsa = os.path.join(mainFolder,dataSet+"_ToxSA_iq.txt")
    toxsi = os.path.join(mainFolder,dataSet+"_ToxSI_iq.txt")
    # out file names 
    sOut = os.path.join(mainFolder,"scleroglossa.constraint")
    tpOut = os.path.join(mainFolder,"toxicoferaP.constraint")
    taiOut = os.path.join(mainFolder,"toxicoferaAI.constraint")
    tsaOut = os.path.join(mainFolder,"toxicoferaSA.constraint")
    tsiOut = os.path.join(mainFolder,"toxicoferaSI.constraint")
    inFiles=(scler,toxpoly,toxai,toxsa,toxsi)
    outFiles=(sOut,tpOut,taiOut,tsaOut,tsiOut)
    # Get initial star tree
    t = getStarTree(allTaxa)

    for x in range(0,5):
        makeTree(inFiles[x],t,outFiles[x],allTaxa)

if __name__=='__main__':
    main()