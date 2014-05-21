import csv, os, time

inputdb = "csv-input/segmentation_bicoid-output.csv"
genome = "/home/barbara/data/genomes/Ofasciatus/Ofas.scaffolds.fa"
blasttype = "tblastn"
output = inputdb[:-4] + "-blastres.csv"

tempin = open('tempin.txt','w') # will be used to store exon sequences temporarily in a file so they can be blasted
exonlist = csv.reader(open(inputdb))


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

ex_in_scaf = open(output,"w")
for exon in exonlist:
	fbid = exon[0]
	genename = exon[1]
	exID = exon[2]
	seq = exon[3]
	tempin = open("tempin.txt","w")
	tempin.write(seq)
	tempin.close()
	blast = "%s -db %s -query tempin.txt -out tempout.txt" %(blasttype,genome)
	os.system(blast)
	tempout = open("tempout.txt")	
	mainlist = blastreader(tempout)
	if len(mainlist) <= 0:
		print "no blastresults for exon", exID
		ex_in_scaf.write("%s,%s,%s\n" %(fbid,genename,exID))
	else:
		scaffold = mainlist[0][0]
		dirx = mainlist[0][1]
		start = mainlist[0][2]
		end = mainlist[0][3]
		ex_in_scaf.write("%s,%s,%s,%s,%s,%s,%s\n" %(fbid,genename,exID,scaffold,dirx,start,end))
	tempout.close()
ex_in_scaf.close()
tempin.close()
os.remove("tempin.txt")
os.remove("tempout.txt")
