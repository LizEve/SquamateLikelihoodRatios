#! /usr/bin/env python
import glob
import os
import time 

def countTaxa(taxaList,files,outfolder):
	mainDir = os.getcwd()
	if not os.path.exists(mainDir+outfolder):
		os.mkdir(mainDir+outfolder)	
	for f in files:
		x = []
		missing = []
		for taxa in taxaList:
			if taxa in open(f).read():
				x.append(1)
			else:
				x.append(0)
				missing.append(taxa)
		print(len(taxaList),sum(x))
		#print(f)
		print(missing)
		if sum(x) == len(taxaList):
			os.system("cp %s .%s" % (f,outfolder))
		#else:
			#print("%s does not have all the taxa" % f)

start = time.time()

files = glob.glob("*.fasta")
taxaList = ["Python_bivittatus","Anolis_brasiliensis", "Anolis_meridionalis", "Tropidurus_oreadicus", "Gekko_japonicus", "Hemidactylus_mabouia", "Gymnodactylus_amarali", "Colobosaura_modesta", "Micrablepharus_maximiliani", "Brasiliscincus_heathi", "Copeoglossum_nigropunctatum", "Notomabuya_frenata", "Ameiva_ameiva", "Ameivula_mumbuca", "Kentropyx_calcarata", "Amphisbaena_alba",  "Gallus_gallus", "Tantilla_melanocephala", "Bothrops_moojeni", "Typhlops_brongersmianus",  "Ophisaurus_gracilis"]
outfolder = "/completeLociFiles"
countTaxa(taxaList,files,outfolder)

end = time.time()
print(end-start)
print((end-start)/60)
print((end-start)/(60*60))