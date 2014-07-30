import sys

if len(sys.argv) <= 4:
	sys.exit("USAGE: python getseq.py [species (4 letter abbreviation)] [scaffold (number only)] [start] [end]")

species = sys.argv[1]
scaf = sys.argv[2]
start=int(sys.argv[3])
end=int(sys.argv[4])

dir_out = "/home/barbara/Dropbox/oncopeltus/annotations"
genomedir = "/home/barbara/data/genomes/"
genomedict = {"Ofas":"Ofasciatus/Ofas.scaffolds.fa", "Clec":"Clectularius/Clec_Bbug02212013.genome.fa"}

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

output = "%s/%s_Scaffold%s:%s..%s.txt" %(dir_out,species,scaf,start,end)
if species in genomedict:
	genome = genomedir + genomedict[species]
else: 
	sys.exit("Species genome not found. Add the path to the python script, or check your spelling.")
readgenome = open(genome)
outputfile = open(output,"w")
scaffold = "Scaffold" + scaf
scafcollect = scaffoldextract(readgenome,scaffold)
query = scafcollect[start-1:end]
outputfile.write(query)
outputfile.close()
