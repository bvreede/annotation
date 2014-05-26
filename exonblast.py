import csv, os, time, sys

if len(sys.argv) <= 2:
	sys.exit("USAGE: python exonblast.py path/to/inputfile (result from exonfinder_prot.py) path/to/genome")

inputdb = sys.argv[1] # input file
inputname = inputdb.split("/")[-1]

genome = sys.argv[2] # genome eg. "/home/barbara/data/genomes/Ofasciatus/Ofas.scaffolds.fa"
blasttype = "tblastn"
csvfolder = "csv"
metafolder = "results"
blastoutput = "blastoutput"
output = csvfolder + "/" + inputname[:-10] + "-blast.csv"
exonlist = csv.reader(open(inputdb))
ex_in_scaf = open(output,"w")

if os.path.exists(blastoutput):
	pass
else:
	os.mkdir(blastoutput)
if os.path.exists(metafolder):
	pass
else:
	os.mkdir(metafolder)

'''
MODULE TAKEN FROM READBLAST.PY:
reads the output file and returns info needed to extract from the
genome fasta file:
* Scaffold number
* frame (not necessary to extract, but useful info for outputfile)
* Start site
* End site
'''
def blastreader(blast):
	mainlist = [] #list of lists: scaffolds and frames with their corresponding start-end sites
	numlist = [] #list of all start-end sites per scaffold and frame
	scaf = ""
	frame = ""
	for line in blast:
		if line[0] == ">": # indicates the start of a new scaffold
			if len(numlist) != 0:
				start = numlist[0]
				end = numlist[-1]
				framelist = [scaf,frame,start,end]
				mainlist.append(framelist)
				numlist = []
			scaf = line.split()[1]
		elif line[1:6] == "Frame": #indicates the start of a new blast result
			if len(numlist) != 0:
				start = (numlist[0])
				end = (numlist[-1])
				framelist = [scaf,frame,start,end]
				mainlist.append(framelist)
				numlist = []
			frame = line.split()[2]
		elif line[:5] == "Sbjct": #indicates a line with subject information
			start = line.split()[1]
			end = line.split()[3]
			numlist.append(start)
			numlist.append(end)
		elif line[:6] == "Lambda": #this indicates the end of the file.
			if len(numlist) != 0:
				start = numlist[0]
				end = numlist[-1]
				framelist = [scaf,frame,start,end]
				mainlist.append(framelist)
			break
	return mainlist

for exon in exonlist:
	fbid = exon[0]
	genename = exon[1].replace(' ','_')
	exID = exon[2]
	seq = exon[3]
	genemeta = open("%s/%s-meta.txt" %(metafolder,genename), "a")
	tempin = open("tempin.txt","w")
	tempin.write(seq)
	tempin.close()
	blastout = "blastoutput/%s_%s" %(fbid,exID,)
	blast = "%s -db %s -query tempin.txt -out %s" %(blasttype,genome,blastout)
	os.system(blast)
	blastout = open(blastout)	
	mainlist = blastreader(blastout)
	if len(mainlist) == 0:
		genemeta.write("%s results: None\n" %(exID))
		ex_in_scaf.write("%s,%s,%s\n" %(fbid,genename,exID))
	else:
		for i in range(len(mainlist)):
			scaffold = mainlist[i][0]
			dirx = mainlist[i][1]
			start = mainlist[i][2]
			end = mainlist[i][3]
			ex_in_scaf.write("%s,%s,%s,%s,%s,%s,%s,%s\n" %(fbid,genename,exID,i+1,scaffold,dirx,start,end))
		genemeta.write("%s results: %s\n" %(exID, len(mainlist)))
	genemeta.close()
	blastout.close()
	print "Completed blast for gene '%s', exon '%s' with %s results." %(genename,exID,len(mainlist))
ex_in_scaf.close()
tempin.close()
os.remove("tempin.txt")
