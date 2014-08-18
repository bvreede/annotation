'''
This script can be used in a pipeline to connect the parsed info
from a blast output with the sequence retriever (getseq.py).
The script requires customization of the genome file location.
The script works with NCBI blast 2.2.29+ in python 2.7.
Author: Barbara Vreede
Contact: b.vreede@gmail.com
Date: 18 August 2014
'''

import sys,os

'''
CUSTOMIZE THE FOLLOWING PATHS!
To allow easy switching to different genomes, add any new genomes used
to the 'genomedict' dictionary instead of replacing the existing information.
'''
genomedir = "/home/barbara/data/genomes/" #directory where genomes are located
genomedict = {"Ofas":"Ofasciatus/Ofas.scaffolds.fa", "Clec":"Clectularius/Clec_Bbug02212013.genome.fa"} #name of blastable genome dbs

'''
Specifying input parameters
'''
if len(sys.argv) <= 3:
	sys.exit("USAGE: python parse_seq.py species infolder extra_nt")

species = sys.argv[1]
infolder = sys.argv[2]
extra = sys.argv[3]

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
try:
	int(extra)
except ValueError:
	sys.exit("Last argument ('extrant') should be a number: how many extra nucleotides to add per sequence extraction.")

'''
run the getseq.py script on each line
'''
filelist = []
for filename in os.listdir(infolder):
	if filename[-4:] == ".csv":
		filelist.append(filename)

print filelist
