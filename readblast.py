'''
This script can be used to quickly extract the appropriate sequences from
a large genome sequence (fasta file), using the output from a local blast.
The script has been tested and developed for a genome sequence in scaffolds,
and works with NCBI blast 2.2.29+ in python 2.7.
Author: Barbara Vreede
Contact: b.vreede@gmail.com
Date: 29 April 2014
'''


### INPUT VARIABLES TO BE GIVEN TO THE PROGRAMME ###
#genome = input # the fasta file with the genome (the same one that was used for the initial blast)
blastres = open("testoutput") #input # blast output file
#addbp = input # how many bp before and after the result sequence
#numseq = input # number of hits to look at



### READ OUTPUT FILE ###
def blastreader(blast):
	'''
	reads the output file and returns info needed to extract from the
	genome fasta file:
	* Scaffold number
	* frame (not necessary to extract, but useful info for outputfile)
	* Start site
	* End site
	'''
	mainlist = [] #list of lists: scaffolds with their corresponding start-end sites
	numlist = [] #list of all start-end sites per scaffold
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
				start = numlist[0]
				end = numlist[-1]
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

fragments = blastreader(blastres)




