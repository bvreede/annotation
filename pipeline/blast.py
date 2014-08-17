'''
This script can be used in an automated pipeline to blast all sequences
in a folder to a genome. The genome needs to be converted to a blastable
database before starting, and its location needs to be specified in this
script.
The script works with NCBI blast 2.2.29+ in python 2.7.
Author: Barbara Vreede
Contact: b.vreede@gmail.com
Date: 17 August 2014
'''

import os, sys

'''
CUSTOMIZE THE FOLLOWING PATHS!
To allow easy switching to different genomes, add any new genomes used
to the 'genomedict' dictionary instead of replacing the existing information.
'''
genomedir = "/home/barbara/data/genomes/" #directory where genomes are located
genomedict = {"Ofas":"Ofasciatus/Ofas.scaffolds.fa", "Clec":"Clectularius/Clec_Bbug02212013.genome.fa"} #name of blastable genome dbs


'''
INPUT PARAMETERS
taken from the command line and subsequently modified
'''
if len(sys.argv) <= 3:
	sys.exit("USAGE: python blast.py species path/to/inputfolder blasttype")

species = sys.argv[1]
infolder = sys.argv[2]
blasttype = sys.argv[3]

if infolder[-1:] == "/": #to prevent mishaps when folder name is specified with / at the end
	infolder = infolder[:-1]

outfolder = "%s_blastout" %(infolder)
genome = genomedir + genomedict[species]

if os.path.exists(outfolder):
	try:
		input("Adding (and possibly overwriting) content to %s. Press enter to confirm..." %outfolder)
	except SyntaxError:
		pass
else:
	os.system("mkdir %s" %outfolder)

'''
Specifying files to blast.
'''
filelist = []
for filename in os.listdir(infolder):
	if filename == ".DS_Store": #this file sometimes gives errors when running on mac
		continue
	filelist.append(filename)

'''
Executing the blast command.
'''
for i in filelist:
	out = "%s_blast2%s" %(i,species)
	blast = "%s -db %s -query %s/%s -out %s/%s" %(blasttype,genome,infolder,i,outfolder,out)
	os.system(blast)

