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


def align(to_align,scaffold,direction,genemeta)


'''
From matrix of results (curated blastresults: only those
that are actual blast hits -- 'newres' from the scaffoldfind
module -- are used) and the most common scaffold it isolates
only those hits which are against this scaffold.
Also checks directions, removes errors in directions, and posts
info to the meta file.
'''
def isolateresults(blastresults,scaffold,genemeta):
	newres2 = []
	dirx = []
	for r in blastresults:
		if r[4] == scaffold:
			newres2.append(r)
			dirx.append(r[5][0])
	if dirx.count('+') == (0 or len(dirx)):
		direction = dirx[0]
		pass
	else:
		if dirx.count('+') > dirx.count('-'):
			direction = '+'
		elif dirx.count('-') > dirx.count('+'):
			direction = '-'
		else:
			direction = 'X'
	newres2.append(direction)
	remove = []
	for i in range(len(dirx)):
		if dirx[i] != direction:
			genemeta.write("%s hit %s: reverse orientation, removed.\n" %(newres2[i][2], newres2[i][3]))
			remove.append(i)
	for i in remove[::-1]:
		newres2.pop(i)
	return newres2


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
	blasted = isolateresults(newres,scaffold,genemeta)
	direction = blasted[-1]
	to_align = blasted[:-1]
	align(to_align,scaffold,direction,genemeta)
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
