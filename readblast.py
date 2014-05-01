'''
This script can be used to quickly extract the appropriate sequences from
a large genome sequence (fasta file), using the output from a local blast.
The script has been tested and developed for a genome sequence in scaffolds,
and works with NCBI blast 2.2.29+ in python 2.7.
Author: Barbara Vreede
Contact: b.vreede@gmail.com
Date: 29 April 2014
'''

import os, time

### INPUT VARIABLES TO BE GIVEN TO THE PROGRAMME ###
genome = "/home/barbara/data/genomes/Ofasciatus/Ofas.scaffolds.fa" #input # the fasta file with the genome (the same one that was used for the initial blast)
blastoutput = "testoutput" #input # blast output file
addbp = 100 #input # how many bp before and after the result sequence
#numseq = input # number of hits to look at

### OTHER PREDEFINED VARIABLES AND INFORMATION ###



### READ OUTPUT FILE ###
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
the frame, and how many bp were added.
'''
def hitextract(fragments,dir_out,n):
	for hit in fragments:
		scaffid = hit[0]
		frame = hit[1]
		span = [int(hit[2]),int(hit[3])]
		start = min(span)-n
		end = max(span)+n
		print start,end

'''


		scaffid = scaff[0]
		items = scaff[-1]
		values = scaff[1:-1]
		low = min(values) - distance
		if low <= 0:
			low = 0
		high = max(values) + distance
		readscaffold = open("%s/%s.txt" %(subdir_out, scaffid))
		for seq in readscaffold:
			if high >= len(seq):
				high = len(seq)
			resultfile = open("%s/%s_part.xdna" %(subdir_out, scaffid), "w")
			resultfile.write(seq[low:high])
			resultfile.close()
		decodefile.write("\nResults for %s\n" %scaffid)
		for i in range(items):
			start = scaff[i*2+1] - low
			end = scaff[i*2+2] - low
			score = i+1
			decodefile.write("Score %s:\t%s\t%s\n" %(score, start, end))
		readscaffold.close()
'''


blastres = open(blastoutput)
readgenome = open(genome)
time = int(time.time())
dir_out = "%s_seq-%s" % (blastoutput, time)
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
hitextract(fragments,dir_out,addbp)


