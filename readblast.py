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
genome = input # the fasta file with the genome [the same one that was used for the initial blast]
blastres = input # blast output file
addbp = input # how many bp before and after the result sequence
numseq = input # number of hits to look at



### READ OUTPUT FILE ###
def blastreader(blast):
'''
reads the output file and returns info needed to extract from the
genome fasta file:
* Scaffold number
* Start site
* End site
'''
	#prevline = ""
	#scorecount = 0
	#start = ""
	#startcount = 0
	mainlist = [] #list of lists: scaffolds with their corresponding start-end sites
	#scaflist = [] #list of scaffolds.
	numlist = [] #list of all start-end sites per scaffold


	for line in blast:
		if line[0] == ">":
			scaf = line.split()[1]
		elif line[1:6] == "Frame":
			if len(numlist) != 0:
				numlist.order
				framelist = [scaf,frame,start,end]
				mainlist.append(framelist)
			else:
				continue
			frame = line.split()[2]
		elif line[:5] == "Sbjct":
			if startcount == 0:
				start = line.split()[1]
				numlist.append(int(start))
			startcount += 1
			prevline = line
		elif line[1:6] == "Score":
			scorecount += 1
			startcount = 0
			if prevline == "":
				pass
			else:
				stop = prevline.split()[3]
				numlist.append(int(stop))
		elif line[:6] == "Lambda": #this indicates the end of the file.
			stop = prevline.split()[3]
			numlist.append(int(stop))
			numlist.append(scorecount)
			mainlist.append(numlist)
			break
	mainlist = mainlist[1:] #because the first item in the list is empty as it stores the scorecount 0
	return mainlist, scaflist


'''
		if line[2:10] == "Scaffold":
			scaffold = line.split()
			if scaffold[0] == ">":
				if prevline != "":
					stop = prevline.split()[3]
					numlist.append(int(stop))
					prevline = ""
				numlist.append(scorecount)
				mainlist.append(numlist)
				numlist = []
				numlist.append(scaffold[1])
				scorecount = 0
			else:
				scaflist.append(scaffold[0])#[8:])
'''


