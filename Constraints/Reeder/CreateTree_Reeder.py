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
    dataSet = "Reeder"
    # Get list of all taxa 
    allTaxa={'Notechis_scutatus', 'Feylinia_polylepis', 'Lanthanotus_borneensis', 'Homalopsis_buccata', 'Celestus_enneagrammus', 'Cricosaura_typica', 'Rhineura_floridana', 'Acrochordus_granulatus', 'Sphenomorphus_solomonis', 'Zonosaurus_ornatus', 'Eugongylus_rufescens', 'Heloderma_horridum', 'Bipes_biporus', 'Aeluroscalobates_felinus', 'Pseudopus_apodus', 'Gekko_gecko', 'Diadophis_punctatus', 'Morunasaurus_annularis', 'Tiliqua_scincoides', 'Cylindrophis_rufus', 'Tupinambis_teguixin', 'Xantusia_vigilis', 'Plica_plica', 'Chamaeleo', 'Rhacodactylus_auriculatus', 'Xenosaurus_grandis', 'Boa_constrictor', 'Amphiesma_stolata', 'Thamnophis_marcianus', 'Corytophanes_cristatus', 'Strophurus_ciliaris', 'Diplometopon_zarudnyi', 'Saltuarius_cornutus', 'Trachyboa_boulengeri', 'Pholidobolus', 'Oplurus_cyclurus', 'Calotes_emma', 'Xenodermus_javanicus', 'Alligator', 'Trachylepis_quinquetaeniata', 'Xenochrophis_piscator', 'Dromaius', 'Amphisbaena_fuliginosa', 'Agkistrodon_contortrix', 'Takydromus_ocellatus', 'Dibamus_novaeguineae', 'Basiliscus_basiliscus', 'Anniella_pulchra', 'Varanus_acanthurus', 'Anolis_carolinensis', 'Gambelia_wislizenii', 'Enyalioides_laticeps', 'Lialis_burtonis', 'Anelytropsis_papillosus', 'Lichanura_trivirgata', 'Natrix_natrix', 'Brachymeles_gracilis', 'Callopistes_maculatus', 'Loxocemus_bicolor', 'Typhlops_jamaicensis', 'Podocnemis', 'Azemiops_feae', 'Ungaliophis_continentalis', 'Stenocercus_guentheri', 'Colobosaura_modesta', 'Leptotyphlops', 'Teius_teyou', 'Tropidophis_haetianus', 'Atractaspis_irregularis', 'Crotaphytus_collaris', 'Eryx_colubrinus', 'Naja', 'Uma_scoparia', 'Shinisaurus_crocodilurus', 'Liolaemus_bellii', 'Uta_stansburiana', 'Eublepharis_macularius', 'Elgaria_multicarinata', 'Phrynosoma_platyrhinos', 'Xenopeltis_unicolor', 'Phymaturus_palluma', 'Aspidites_melanocephalus', 'Uranoscodon_superciliosus', 'Leiocephalus_barahonensis', 'Heloderma_suspectum', 'Chelydra', 'Physignathus_cocincinus', 'Bothrops_asper', 'Petrosaurus_mearnsi', 'Laticauda_colubrina', 'Brachylophus_fasciatus', 'Platysaurus', 'Acontias', 'Pareas_hamptoni', 'Homo', 'Coleonyx_variegatus', 'Afronatrix_anoscopus', 'Micrurus_fulvius', 'Bipes_canaliculatus', 'Agama_agama', 'Trogonophis_wiegmanni', 'Mus', 'Leiosaurus_catamarcensis', 'Varanus_exanthematicus', 'Anilius_scytale', 'Coluber_constrictor', 'Lachesis_muta', 'Causus', 'Gallus', 'Teratoscincus', 'Pogona_vitticeps', 'Calabaria_reinhardtii', 'Geocalamus_acutus', 'Polychrus_marmoratus', 'Plestiodon_fasciatus', 'Brookesia_brygooi', 'Lycophidion_capense', 'Liotyphlops_albirostris', 'Casarea_dussumieri', 'Python_molurus', 'Urostrophus_vautieri', 'Sceloporus_variabilis', 'Amphiglossus_splendidus', 'Cordylus_mossambicus', 'Hydrosaurus', 'Xenosaurus_platyceps', 'Epicrates_striatus', 'Chalarodon_madagascariensis', 'Uropeltis_melanogaster', 'Gonatodes_albogularis', 'Varanus_salvator', 'Tachyglos', 'Lepidophyma_flavimaculatu', 'Lampropeltis_getula', 'Scincus', 'Daboia_russelli', 'Pristidactylus_torquatus', 'Leiolepis_belliana', 'Cordylosaurus_subtesselatus', 'Exiliboa_placata', 'Lacerta_viridis', 'Crocodylus', 'Sphenodon_punctatus', 'Delma_borea', 'Aparallactus_werneri', 'Aspidoscelis_tigris', 'Sauromalus_ater', 'Phelsuma_lineata', 'Uromastyx_aegyptus', 'Dipsosaurus_dorsalis'}
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