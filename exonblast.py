import csv, os, time

inputdb = "csv-input/segmentation_short-output.csv"
genome = "/home/barbara/data/genomes/Ofasciatus/Ofas.scaffolds.fa"
blasttype = "tblastn"

time = int(time.time())
dir_out = "output-%s" %(time)
os.mkdir(dir_out)

tempin = open('tempin.txt','w') # will be used to store exon sequences temporarily in a file so they can be blasted
exonlist = csv.reader(open(inputdb))


'''
FROM READBLAST.PY:
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

fbid = ""
for exon in exonlist:
	genename = exon[1]
	exID = exon[2]
	seq = exon[3]
	#tempin.truncate(0)
	tempin = open("tempin.txt","w")
	tempin.write(seq)
	tempin.close()
	print seq
	blast = "%s -db %s -query tempin.txt -out tempout.txt" %(blasttype,genome)
	print blast
	os.system(blast)
	tempout = open("tempout.txt")
	if exon[0] == fbid:
		mainlist = blastreader(tempout)
		#ex_in_scaf.write("%s,%s,%s,%s,%s\n" %(exID,mainlist[0],mainlist[1],mainlist[2],mainlist[3]))
	else:
		#if fbid != "":
			#ex_in_scaf.close() # close old file (only if this is not the first loop)
		fbid = exon[0]
		#genomeloc = "%s/%s_genomeloc.csv" %(dir_out, fbid)
		#ex_in_scaf = open(genomeloc,"w")	
		mainlist = blastreader(tempout)
		#ex_in_scaf.write("%s,%s,%s,%s,%s\n" %(exID,mainlist[0],mainlist[1],mainlist[2],mainlist[3]))
	print genename, exID, mainlist
	tempout.close()

#ex_in_scaf.close()
tempin.close()
#os.remove("tempin.txt")
#os.remove("tempout.txt")


'''
	if exon[0] == fbid:
		# do the stuff you need to do when still working with the same gene
	else:
		# start a new gene
		fbid = exon[0]
'''
