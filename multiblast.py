import os

genome = "/home/barbara/data/genomes/Clectularius/Clec_Bbug02212013.genome.fa"
#genome = "/home/barbara/data/genomes/Ofasciatus/Ofas.scaffolds.fa"
species = "clec"
#species = "ofas"
blastfolder = "/home/barbara/Dropbox/oncopeltus/annotations/"
filelist = ['dmel_odd', 'tcas_odd', 'dmel_bowl', 'tcas_bowl', 'dmel_krup', 'tcas_krup', 'tcas_krup2', 'dmel_sob', 'tcas_sob']
add = 1000

for i in filelist:
	out = i + "-" + species
	blast = "tblastn -db %s -query %s%s -out %s%s" %(genome,blastfolder,i,blastfolder,out)
	os.system(blast)
	readblast = "python readblast.py %s %s%s %s %s" %(genome,blastfolder,out,add,add)
	#os.system(readblast)
	


