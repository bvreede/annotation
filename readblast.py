'''
This script can be used to quickly extract the appropriate sequences from
a large genome sequence (fasta file), using the output from a local blast.
The script has been tested and developed for a genome sequence in scaffolds,
and works with NCBI blast 2.2.29+ in python 2.7.
Author: Barbara Vreede
Contact: b.vreede@gmail.com
Date: 29 April 2014
'''

import os, time, sys

### INPUT ARGUMENTS ###
if len(sys.argv) <= 2:
	sys.exit("USAGE: python readblast.py path/to/genome path/to/blast.output.file [n_nucleotides_to_add_before] [n_nucleotides_to_add_after]")

genome = sys.argv[1] #path to genome
blastoutput = sys.argv[2] #path to BLAST outputfile
if len(sys.argv) == 5:
	addbp_before=int(sys.argv[3])
	addbp_after=int(sys.argv[4])
elif len(sys.argv) == 4:
	addbp_before=int(sys.argv[3])
	addbp_after=int(sys.argv[3])
else:
	addbp_before = 100 #n of nucleotides to add before each sequence
	addbp_after = 100


### READ BLAST OUTPUT FILE ###
'''
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


### EXTRACT THE APPROPRIATE SCAFFOLDS ###
'''
Goes into genome fasta file and collects the scaffolds that
are in the scaffold list; saves them as separate fasta files.
'''
def scaffoldextract(readgenome,dir_out,scaflist):
	marker = 0 # is "1" when reading and saving scaffold sequence
	scafname = "notinlist" # is "notinlist" when scaffold name is not in list
	for line in readgenome:
		if line[0] == ">":
			scaf = line[1:].strip()
			if scaf in scaflist:
				marker = 1
				scafname = scaf
				outputfile = open("%s/%s-entire.fa" %(dir_out,scafname), "w")
				outputfile.write(line)
			else:
				marker = 0
				if scafname != "notinlist":
					outputfile.close()
					scafname = "notinlist"
		elif marker == 1:
			outputfile.write(line)#.strip())
			
### GO THROUGH THE SAVED SCAFFOLDS AND SELECT THE APPROPRIATE AREA ###
'''
Takes each hit and searches the appropriate scaffold file.
Finds the sequence that spans between start-end, and adds up
the extra span length that has been indicated by the user.
(Default is 100.)
Makes a new fasta file that contains the sequences of all hits
and describes in the header the start and end site, as well as
the frame, and how many nt were added.
'''
def hitextract(fragments,dir_out,b,a):
	for hit in fragments:
		scaffid = hit[0]
		frame = hit[1]
		span = [int(hit[2]),int(hit[3])]
		start = min(span)-b
		if start <= 1:
			start = 1
		end = max(span)+a
		scaffasta = open("%s/%s-entire.fa" %(dir_out,scaffid))
		line = scaffasta.readlines()
		linelen = len(line[1].strip()) # calculates the n of basepairs in each line
		linestart = start/50+1
		lineend = end/50+1
		if lineend >= len(line):
			lineend = len(line)
		start_inseq = min(span)-((linestart-1)*50)+1 # where in the sequence between linestart-lineend is the actual hit?
		end_inseq = max(span)-((linestart-1)*50)+1 # where in the sequence between linestart-lineend does the hit end?
		outfasta = open("%s/%s-blastresults.fa" %(dir_out,scaffid), "a")
		outfasta.write(">%s:%s-%s (hit: %s-%s in this sequence) frame: %s\n" %(scaffid,hit[2],hit[3],start_inseq,end_inseq,frame[0]))
		for i in range(linestart,lineend):
			outfasta.write(line[i])
		outfasta.write("\n")
		#results = open("%s/HIT-TABLE.txt" %(dir_out), "a")
		#results.write("%s\tstart: %s\tend: %s\tframe: %s\n" %(scaffid,hit[2],hit[3],frame))

blastres = open(blastoutput)
readgenome = open(genome)
time = int(time.time())
dir_out = "%s-%s" % (blastoutput, time)
os.mkdir(dir_out)
fragments = blastreader(blastres)
scaflist = []
#make a list of scaffolds to collect
for line in fragments:
	if line[0] in scaflist:
		continue
	else:
		scaflist.append(line[0])
scaffoldextract(readgenome,dir_out,scaflist)
hitextract(fragments,dir_out,addbp_before,addbp_after)
