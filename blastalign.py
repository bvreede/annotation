import csv, sys, os


if len(sys.argv) <= 1:
	sys.exit("USAGE: python blastalign.py path/to/inputfile (result from exonblast.py)")

genome = "/home/barbara/data/genomes/Ofasciatus/Ofas.scaffolds.fa"
inputdb = sys.argv[1] # input file
blastlist = csv.reader(open(inputdb))
outputfolder = "alignments"
scaffolder = genome.split('/')[-1][:-3] # name of genome fasta file minus .fa

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
scaffold; saves the scaffold sequence as separate fasta file.
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

'''
Takes a matrix of blastresults from one gene (per result:
[fbid,genename,exon,resultno,scaffold,direction,start,end])
and finds the most common scaffold.
Also calls the scaffoldextract module and saves the scaffold
as a separate file.
'''
def scaffoldfind(blastresults):
	scaffolds = []
	scafdict = {}
	newres = []
	genemeta = open("blastmeta/%s_meta.txt" %(blastresults[0][0]), "a")
	genemeta.write("\n\n++++SCAFFOLDS ISOLATED:++++\n\n")
	for r in blastresults:
		if len(r) > 4:
			scaffolds.append(r[4])
			newres.append(r) # ensures that only completed blasts are included in 'newres' matrix
	for s in scaffolds: # put scaffold name and their occurrence count in a dictionary
		n = scaffolds.count(s)
		scafdict[s] = n
	v = list(scafdict.values()) # list of scaffold occurrences (order as in dictionary)
	k = list(scafdict.keys()) # list of scaffold names (order as in dictionary)
	scaffold = k[v.index(max(v))] # finds the FIRST scaffold with the highest occurrence
	scaffoldextract(scaffold)
	genemeta.write("Most common scaffold: %s\n" %(scaffold))
	'''
	FOLLOWING BIT WRITTEN TO ENSURE ALL EXONS ARE INCLUDED IN SCAFFOLD
	BUT MAY BE SUPERFLUOUS
	v.remove(max(v))
	k.remove(scaffold)
	scaffold2 = k[v.index(max(v))]
	scaffoldextract(scaffold2)
	'''
	exon = ""
	scafcheck = []
	for t in newres:
		if t[2] != exon:
			if scaffold not in scafcheck and len(scafcheck) > 0:
				genemeta.write("Exon %s blasts against %s\n" %(exon,scafcheck))
			exon = t[2]
			scafcheck = []
			scafcheck.append(t[4])
		else:
			scafcheck.append(t[4])
	if scaffold not in scafcheck and len(scafcheck) >0:
		genemeta.write("Exon %s blasts against %s\n" %(exon,scafcheck))
	genemeta.close()
	return scaffold
	

fbid = ""
blastresults = []
for line in blastlist:
	if line[0] != fbid:
		if len(blastresults) != 0:
			scaffold = scaffoldfind(blastresults)
		blastresults = []
		output.close()
		fbid = line[0]
		output = open("%s/%s-align.csv" %(outputfolder,fbid),"w")
		blastresults.append(line)
	else:
		blastresults.append(line)
scaffold = scaffoldfind(blastresults)

output.close()
os.remove("x.txt") # remove decoy outputfile
