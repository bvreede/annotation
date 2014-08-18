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
	sys.exit("USAGE: python parse_seq.py species infolder extra_nt max_scaffolds")

species = sys.argv[1]
infolder = sys.argv[2]
extra = int(sys.argv[3])
max_scafcount = int(sys.argv[4])

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
translate scaffold, start, and end site into a system command
that retrieves the corresponding sequence and puts it in the appropriate
output folder.
'''
scafcount = 0

def getseq(scaf_e,start_e,end_e,outfolder):
	start_b = start_e - extra
	end_b = end_e + extra
	scaf_b = scaf_e[8:]
	command = "python getseq.py %s %s %s %s %s" %(species,scaf_b,start_b,end_b,outfolder)
	global scafcount #variable that counts whenever a 
	scafcount += 1
	if scafcount <= max_scafcount:
		os.system(command)

'''
run the getseq.py script on each line
'''
filelist = []
for filename in os.listdir(infolder):
	if filename[-4:] == ".csv":
		filelist.append(filename)


for k in filelist:
	infile = open("%s/%s" %(infolder,k))
	outfolder = "%s/%s" %(infolder,k[:-4])
	if os.path.exists(outfolder):
		pass
	else:
		os.system("mkdir %s" %outfolder)
	count = 0
	for line in infile:
		info = line.split(',')
		if len(info) <= 4 or info[0] == "Scaffold": # this indicates the headers and meta info in the file
			continue
		#only actual hits from this point on!
		scaffold = info[0] #saves scaffold number
		start = int(info[5]) #saves start number
		end = int(info[6]) #saves end number
		if count != 0: #count = 0 only in the first hit, where scaf_e etc. are not yet specified.
			if scaffold == scaf_e:
				start_e = min(start_e,start)
				end_e = max(end_e,end)
			else:
				getseq(scaf_e,start_e,end_e,outfolder)
				scaf_e = scaffold
				start_e = start
				end_e = end	
		else: #this condition should only be passed in the first hit.
			scaf_e = scaffold
			start_e = start
			end_e = end
		count += 1 #increases with every line
	getseq(scaf_e,start_e,end_e,outfolder)
	scafcount = 0
