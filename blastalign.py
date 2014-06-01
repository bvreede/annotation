import csv, sys, os


if len(sys.argv) <= 2:
	sys.exit("USAGE: python blastalign.py path/to/inputfile (result from exonblast.py) path/to/genome")

inputdb = sys.argv[1] # input file
genome = sys.argv[2] # genome (eg: "/home/barbara/data/genomes/Ofasciatus/Ofas.scaffolds.fa"
blastlist = csv.reader(open(inputdb))
outputfolder = "results"
scaffolder = genome.split('/')[-1][:-3] # name of genome fasta file minus .fa
extra = 500 # n basepairs to add to scaffold sequence (NB! the last line will not be printed so the last alignment will be less than n basepairs from the end!)
lline = 90 # adjustable number for linelength in resultfile

if os.path.exists(outputfolder):
	pass
else:
	os.mkdir(outputfolder)
if os.path.exists(scaffolder):
	pass
else:
	os.mkdir(scaffolder)

'''
Alignment module: takes all selected blast hits
and the scaffold ID, and aligns hits to the scaffold
sequence.
Adjustable variables: additional scaffold sequence (before 
and after the first/last alignment) as 'extra'; length
of each line as 'lline'. Both are found in the definitions
on the top of this script.

Main input database is called 'to_align' and consists of:
[0] = fbid
[1] = genename
[2] = exonID
[3] = exonnr
[4] = scaffold
[5] = frame
[6] = start in genome
[7] = end in genome
[8] = start in exon (query)
[9] = end in exon (query)
[10] = length of exon (query)
(except when no hits are found, then it's only 0-2)
'''
def align(to_align,scaffold,genemeta):
	fbid = to_align[0][0]
	genename = to_align[0][1].replace(' ','_')
	genemeta.write("\n\n++++ALIGNMENTS:++++\n\n")
	genemeta.write("exon (no.)\t\tframe\tstart\tend\n")
	for l in to_align:
		genemeta.write("%s (%s)\t%s\t%s\t%s\n" %(l[2],l[3],l[5],l[6],l[7]))
	genemeta.close()
	output = open("%s/%s-align.txt" %(outputfolder,genename),"w")
	scaffile = open("%s/%s.fa" %(scaffolder,scaffold))
	scafseq = ''
	for line in scaffile:
		if line[0] != '>':
			scafseq += line.strip()
	startend = []
	for p in to_align:
		startend.append(int(p[6]))
		startend.append(int(p[7]))
	start = min(startend) - extra
	end = max(startend) + extra
	scafseq2 = scafseq[(start):end]
	for i in range(len(scafseq2)/lline):
		startline = i*lline + start # absolute sequence start site (in scaffold)
		endline = (i+1)*lline + start
		for e in to_align:
			emin = min(int(e[6]),int(e[7]))
			emax = max(int(e[6]),int(e[7]))
			tag = "%s-%s (%s)" %(e[2],e[3],e[5])
			space1=(28 - len(tag))*' ' # adds spaces to 'tag' to a total of 28
			if emin >= startline and emin < endline:
				space2=(emin-startline)*' '
				align=(lline-len(space2))*'+'
				output.write('%s%s %s%s\n' %(tag,space1,space2,align))
			elif emax >= startline and emax < endline:
				align=(emax-startline)*'+'
				output.write('%s%s %s\n' %(tag,space1,align))
			elif emin <= startline and emax >= endline:
				align=lline*'+'
				output.write('%s%s %s\n' %(tag,space1,align))
		space = (28 - len(str(startline)))*' '# adds spaces to 'startline' to a total of 28
		output.write('\n%s%s %s %s\n' %(space, startline, scafseq2[i*lline:(i+1)*lline], endline-1))
	scaffile.close()
	output.close()

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
	if dirx.count('+') == (0 or len(dirx)) and len(dirx) > 0:
		direction = dirx[0]
		pass
	else:
		if dirx.count('+') > dirx.count('-'):
			direction = '+'
		elif dirx.count('-') > dirx.count('+'):
			direction = '-'
		else:
			direction = 'X'
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
and finds the most common scaffold and gene orientation.
Also calls the scaffoldextract module and saves the scaffold
as a separate file.
This is the main module of the script: it calls both the
scaffold extractor (scaffoldextract), blastresult isolator
(isolateresults) and the alignment module (align).
'''
def scaffoldfind(blastresults):
	scaffolds = []
	scafdict = {}
	newres = []
	genename = blastresults[0][1].replace(' ','_')
	print "Aligning blast hits for gene '%s'..." %(genename)
	genemeta = open("%s/%s-meta.txt" %(outputfolder,genename), "a")
	genemeta.write("\n\n++++RAW HITS:++++\n\n")
	genemeta.write("exon ID\t\t\tlength\thit no.\tstart\tend\tscaffold\tframe\tstart\tend\n")
	for hit in blastresults:#adds meta-info about each hit to the meta file.
		if len(hit) > 4:
			genemeta.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(hit[2],hit[10],hit[3],hit[8],hit[9],hit[4],hit[5],hit[6],hit[7]))
		else:
			genemeta.write("None.")
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
	if len(newres) == 0:
		genemeta.write("No blastresults, no alignments!\n\n")
		scaffold = "No-scaffold"
		print "No blastresults to align for gene '%s'.\nContinuing..." %(genename)
	else:
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
			scafcheck.append(t[4]) # makes list of all scaffolds to check that they are the same
		else:
			scafcheck.append(t[4])
	if scaffold not in scafcheck and len(scafcheck) >0:
		genemeta.write("Exon %s blasts against %s\n" %(exon,scafcheck))
	to_align = isolateresults(newres,scaffold,genemeta)
	if len(to_align) > 0:
		align(to_align,scaffold,genemeta)
	else:
		genemeta.write("ERROR determining scaffold and results to add in blast! Aborted alignment.")


fbid = ""
blastresults = []
for line in blastlist: # the following loop ensures that scaffoldfind is called per gene, not per entry.
	if line[0] != fbid:
		if len(blastresults) != 0:
			scaffoldfind(blastresults)
		blastresults = []
		fbid = line[0]
		blastresults.append(line)
	else:
		blastresults.append(line)
scaffoldfind(blastresults)

