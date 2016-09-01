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
genomedict = {"Ofas":"Ofasciatus/Ofas.scaffolds.fa", "Clec":"Clectularius/Clec_Bbug02212013.genome.fa", "Ccap":"Ccapitata/Ccap01172013-genome.fa", "Otau":"Otaurus/Otaur.scaffolds.fa"} #name of blastable genome dbs


'''
INPUT PARAMETERS
taken from the command line and subsequently modified
'''
if len(sys.argv) <= 3:
	sys.exit("USAGE: python blast.py species path/to/inputfolder blasttype")

species = sys.argv[1]
infolder = sys.argv[2]
blasttype = sys.argv[3]

'''
Verifying input parameters
'''
if species in genomedict:
	genome = genomedir + genomedict[species]
else: 
	sys.exit("Species genome not found. Add the path to the getseq.py script, or check your spelling. (e.g. 'Ofas'.)")
if os.path.exists(genome):
	readgenome = open(genome)
else:
	sys.exit("Could not find genome directory. Verify path in getseq.py code.")

if infolder[-1:] == "/": #to prevent mishaps when folder name is specified with / at the end
	infolder = infolder[:-1]
if os.path.exists(infolder):
	pass
else:
	sys.exit("Could not find input directory. Verify path in input parameters.")

outfolder = "%s_blast2%s" %(infolder,species)
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
	if filename[0] == ".": #So not to run on hidden files...
		continue
	filelist.append(filename)

'''
Executing the blast command.
'''
for i in filelist:
	out = "%s_blast2%s.txt" %(i,species)
	blast = "%s -db %s -query %s/%s -out %s/%s" %(blasttype,genome,infolder,i,outfolder,out)
	os.system(blast)
