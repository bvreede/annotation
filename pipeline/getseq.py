'''
This script can be used in an automated pipeline to retrieve
a fragment of a genome sequence (in a fasta file) and save it
separately.
Requires customization of directory paths in the script: location
of genome fasta files and the output folder.
The script works with python 2.7.
Author: Barbara Vreede
Contact: b.vreede@gmail.com
Date: 17 August 2014
'''

import sys,os

if len(sys.argv) <= 4:
	sys.exit("USAGE: python getseq.py species scaffoldnr start end [optional:outfolder]")

species = sys.argv[1]
scaf = sys.argv[2]
start=int(sys.argv[3])
end=int(sys.argv[4])

'''
Customize the following paths:
'''
try:
	dir_out = sys.argv[5]
except IndexError: #if there are only 4 arguments then the outfolder will be default.
	dir_out = "/home/barbara/Dropbox/oncopeltus/annotations"
genomedir = "/home/barbara/data/genomes/"
genomedict = {"Ofas":"Ofasciatus/Ofas.scaffolds.fa", "Clec":"Clectularius/Clec_Bbug02212013.genome.fa", "Ccap":"Ccapitata/Ccap01172013-genome.fa"}

'''
Verifying input data, open genome file and create output folder.
'''
if species in genomedict:
	genome = genomedir + genomedict[species]
else: 
	sys.exit("Species genome not found. Add the path to the getseq.py script, or check your spelling. (e.g. 'Ofas'.)")
if os.path.exists(genome):
	readgenome = open(genome)
else:
	sys.exit("Could not find genome directory. Verify path in getseq.py code.")
for i in range(len(scaf)):
	scaftest = scaf[i:]
	try:
		int(scaftest) #checks if the 'scaffold' input is a number. Removes characters until it is.
		scaf = scaftest
		break
	except ValueError:
		continue
try: #last check: if after all the character-removing the scaffold *still* isn't a number...
	int(scaf)
except ValueError:
	sys.exit("USAGE: python getseq.py species scaffoldnr start end [optional:outfolder]. Scaffold not a number.")
if os.path.exists(dir_out):
	pass
else:
	os.system("mkdir %s" %dir_out)

'''
Find the appropriate scaffold in the fasta file, and save
the entire scaffold sequence to local memory.
'''
def scaffoldextract(readgenome,scaffold):
	scafcollect = ""
	marker = 0 # is "1" when reading and saving the correct scaffold sequence
	for line in readgenome:
		if line[0] == ">":
			header = line[1:].strip()
			if header == scaffold:
				marker = 1
				print "Found " + scaffold + ", saving sequence..."
			else:
				marker = 0
		elif marker == 1:
			scafcollect += line.strip()
	return scafcollect

'''
Open the outputfile, call scaffoldextract, save the appropriate
part of the scaffold to the outputfile.
'''
scaffold = "Scaffold" + scaf
scafcollect = scaffoldextract(readgenome,scaffold)
if start <= 0:
	start = 1
if end >= len(scafcollect):
	end = len(scafcollect)
output = "%s/%s_Scaffold%s:%s..%s" %(dir_out,species,scaf,start,end)
outputfile = open(output,"w")
query = scafcollect[start-1:end]
outputfile.write(query)
outputfile.close()
