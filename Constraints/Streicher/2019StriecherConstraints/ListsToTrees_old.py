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

def writeTaxa(globber,folder):
    cwd = os.getcwd()
    os.chdir(folder)
    allTaxa=[]
    for f in glob.glob(globber):
        #print(f)
        alignment = AlignIO.read(open(f),"phylip-relaxed")
        for record in alignment:
            allTaxa.append(str(record.id))
    os.chdir(cwd)
    return list(set(allTaxa))



def getStarTree(allTaxa):
    """Takes
    Returns"""
    # Create empty tree
    t1 = ete3.Tree()
    
    # Iteratively add all taxa
    for tx in allTaxa:
        t1.add_child(name=tx) 
        
    return t1


def checkConflict(nextbipart,addedBiparts):
    """ Takes list of bipartitions already in tree, and new bipartition to be added
    Returns number of bipartitions in tree that conflict with new bipartition.
    Returns number of bipartitions that include the new bipart.
    """
    # Counter for number of conflicts
    conflicts=0

    # Iterate through current bipartitions
    for bipart in addedBiparts:
        
        # Get number of items that overlap between bipartitions
        overlap = len(set(nextbipart)&set(bipart))
        
        # Get length of smaller bipartition
        smaller = min([len(nextbipart),len(bipart)])
        
        # Check if one bipartition is contained in the other, ie the overlap == smaller bipartition
        if overlap == smaller:
            conflicts+=0
            
        # No taxa are shared
        elif overlap == 0:
            conflicts+=0
            
        # Some but not all taxa are shared, indicating conflict
        elif overlap < smaller:
            conflicts+=1
            
        else:
            print("bad, wrong, uh oh.")
            
    return conflicts

def addBipart(t,allTaxa,newBipart):
    """ Takes a tree, a bipartition, and a list of all taxa in the tree
    Returns new tree
    Has bug - say if AB is new bipart, but is already in clade ABCD, AB will be added as new bipartition, and the ABCD will not be preserved. 
    Solution - add biparts from fewest to most taxa. 
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

def pickTreeAddBipart(newBipart,t1,addedBipartsT1,addedBipartsNeither,allTaxa):
    """ Takes bipartition to be added, trees, lists of biparts already added. 
    Returns updated trees and lists
    Uses defs: checkConflict, addBipart
    """

    # Check for conflicts for first tree
    conflicts = checkConflict(newBipart,addedBipartsT1)
    
    # If conflict for first tree
    if conflicts > 0:
        addedBipartsNeither.append(newBipart)
            
    # No conflicts, add to first tree
    else:
        # Add to bipart list
        addedBipartsT1.append(newBipart)
        
        # Add to tree
        t1=addBipart(t1,allTaxa,newBipart)
        
    return t1,addedBipartsT1,addedBipartsNeither

def fileToTree(fileName,allTaxa,t1,constraintsTree,outLog):

    # Create list with biparts added 
    addedBipartsT1 = []
    addedBipartsNeither = []
    lenDict = {}

    # Open file 
    for line in open(fileName):
        # A bunch of cleanup that I will do in R or textwrangler if I have to redo this
        line = line.strip("\n")
        line = line.strip(")")
        line = line.strip("(")
        # Make into list 
        bipart = line.split(",")
        # Clean up spaces in names 
        bipart = [x.replace(' ','') for x in bipart]
        # Write bipart and number of taxa in bipart 
        lenDict[str(bipart)]=[int(len(bipart)),bipart]

    # Order biparts by number of taxa
    # Have to do this because of bipart adding fuction, more details in notes for that function. 
    byTaxaDict = OrderedDict(sorted(lenDict.items(), key = itemgetter(1), reverse = False))
    # For each bipartition in dictionary, add to tree
    for key,value in byTaxaDict.items():
        bipart = value[1]
        t1,addedBipartsT1,addedBipartsNeither=pickTreeAddBipart(bipart,t1,addedBipartsT1,addedBipartsNeither,allTaxa)

    # Convert to dendropy tree list to write to file
    # Ete3 newick formatting: 5 = internal and leaf branches + leaf names, 9 = leaf names
    dt = dendropy.Tree.get(data=t1.write(format = 9),schema = 'newick')

    # Output tree
    out=open(constraintsTree, "w")
    dt.write(file=out, schema='newick')
    out.close()

    # Write to main log file
    
    outLog.write(str(fileName)+" "+str(len(addedBipartsT1))+" biparts"+"\n")
    outLog.write(dt.__str__()+"\n\n")
    outLog.write("Community extra biparts,"+str(len(addedBipartsNeither))+":\n")
    outLog.write(str(addedBipartsNeither)+"\n\n\n")


def main():
    nexFolder = "/Users/ChatNoir/Projects/Squam/Streicher/StreicherData"
    constraintsFile = "/Users/ChatNoir/Projects/Squam/scripts/Constraints/Streicher_IQ_sclero.txt" 
    constraintsTree = '/Users/ChatNoir/Projects/Squam/scripts/Constraints/Streicher_IQ_sclero_tree.txt'
    constraintsFile1 = "/Users/ChatNoir/Projects/Squam/scripts/Constraints/Streicher_IQ_toxpoly.txt" 
    constraintsTree1 = '/Users/ChatNoir/Projects/Squam/scripts/Constraints/Streicher_IQ_toxpoly_tree.txt'
    constraintsFile2 = "/Users/ChatNoir/Projects/Squam/scripts/Constraints/Streicher_IQ_toxai.txt" 
    constraintsTree2 = '/Users/ChatNoir/Projects/Squam/scripts/Constraints/Streicher_IQ_toxai_tree.txt'
    constraintsFile3 = "/Users/ChatNoir/Projects/Squam/scripts/Constraints/Streicher_IQ_toxsa.txt" 
    constraintsTree3 = '/Users/ChatNoir/Projects/Squam/scripts/Constraints/Streicher_IQ_toxsa_tree.txt'
    constraintsFile4 = "/Users/ChatNoir/Projects/Squam/scripts/Constraints/Streicher_IQ_toxsi.txt" 
    constraintsTree4 = '/Users/ChatNoir/Projects/Squam/scripts/Constraints/Streicher_IQ_toxsi_tree.txt'

    # Get list of all taxa 
    aT = writeTaxa("*phy",nexFolder)
    print(len(aT))
    # Get initial star tree
    t = getStarTree(aT)
    outLog = open("log.log","w+")
    fileToTree(constraintsFile,aT,t,constraintsTree,outLog)
    fileToTree(constraintsFile1,aT,t,constraintsTree1,outLog)
    fileToTree(constraintsFile2,aT,t,constraintsTree2,outLog)
    fileToTree(constraintsFile3,aT,t,constraintsTree3,outLog)
    fileToTree(constraintsFile4,aT,t,constraintsTree4,outLog)
    outLog.close()


if __name__=='__main__':
    main()