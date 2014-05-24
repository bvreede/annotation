import csv, sys, os


if len(sys.argv) <= 1:
	sys.exit("USAGE: python blastalign.py path/to/inputfile (result from exonblast.py)")

genome = "/home/barbara/data/genomes/Ofasciatus/Ofas.scaffolds.fa"
inputdb = sys.argv[1] # input file
blastlist = csv.reader(open(inputdb))
outputfolder = "alignments"
scaffolder = genome.split('/')[-1][:-3]

output = open("x.txt","w") # decoy outputfile so that close statement can be used

if os.path.exists(outputfolder):
	pass
else:
	os.mkdir(outputfolder)
if os.path.exists(scaffolder):
	pass
else:
	os.mkdir(scaffolder)

'''
MODULE TAKEN FROM READBLAST.PY (and modified):
Goes into genome fasta file and collects the appropriate 
scaffold; saves as separate fasta files.
'''
def scaffoldextract(scaffold):
	marker = 0 # is "1" when reading and saving scaffold sequence
	readgenome = open(genome)
	for line in readgenome:
		if line[0] == ">":
			scaf = line[1:].strip()
			if scaf == scaffold:
				marker = 1
				outputfile = open("%s/%s.fa" %(scaffolder,scaffold), "w")
				outputfile.write(line)
			else:
				if marker == 1:
					outputfile.close()
					readgenome.close()
					break
		elif marker == 1:
			outputfile.write(line.strip())

def scaffoldfind(blastresults,output):
	scaffolds = []
	scafdict = {}
	newres = []
	for r in blastresults:
		if len(r) > 4:
			scaffolds.append(r[4])
			if type(r[2]) == list: 
				r[2].append(r[4])
			else:
				r[2] = [r[4]]
			print r[2]
	for s in scaffolds:
		n = scaffolds.count(s)
		scafdict[s] = n
	v = list(scafdict.values())
	k = list(scafdict.keys())
	scaffold = k[v.index(max(v))]
	scaffoldextract(scaffold)
	'''
	exon = ""
	scafcheck = []
	for t in newres:
		if t[2] != exon:
			if scaffold not in scafcheck and len(scafcheck) > 0:
				
			exon = t[2]
			scafcheck = []
			scafcheck.append(t[4])
		else:
			scafcheck.append(t[4])
	if scaffold not in scafcheck and len(scafcheck) >0:
	'''	
	return scaffold
	

fbid = ""
blastresults = []
for line in blastlist:
	if line[0] != fbid:
		if len(blastresults) != 0:
			scaffold = scaffoldfind(blastresults,output)
		blastresults = []
		output.close()
		fbid = line[0]
		output = open("%s/%s-align.csv" %(outputfolder,fbid),"w")
		blastresults.append(line)
	else:
		blastresults.append(line)
scaffold = scaffoldfind(blastresults,output)

output.close()
os.remove("x.txt") # remove decoy outputfile
