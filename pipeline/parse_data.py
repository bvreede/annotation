'''
This script can be used in an automated pipeline to parse a blast output
file and collect the metadata for each hit in a csv file.
The script works with NCBI blast 2.2.29+ in python 2.7.
Author: Barbara Vreede
Contact: b.vreede@gmail.com
Date: 17 August 2014
'''

import os,sys

if len(sys.argv) <= 1:
	sys.exit("USAGE: python parse_data.py path/to/outputfolder")
infolder = sys.argv[1]
if infolder[-1:] == "/": #to prevent mishaps when folder name is specified with / at the end
	infolder = infolder[:-1]

'''
Writes the information currently in memory to the outputfile.
'''
def saveinfo(out,scaffold,evalue,frame,querylist,sbjlist):
	qstart,qend = min(querylist),max(querylist)
	sstart,send = min(sbjlist),max(sbjlist)
	wa = scaffold + "%3A"
	out.write("%s,%s,%s,%s,%s,%s,%s,%s%s..%s\n" %(scaffold,evalue,frame,qstart,qend,sstart,send,wa,sstart,send))
	querylist = []
	sbjlist = []
	return querylist,sbjlist

'''
Reads the blast output for a given filename, opens and saves the
appropriate outputfile.
Calls saveinfo for each line (=hit).
'''
def blastreader(filename):
	blast = open("%s/%s" %(infolder,filename))
	outname = filename[:-4] + "_smy.csv"
	out = open("%s/%s" %(infolder,outname), "w")
	#define parameters (empty or null)
	length,score = 0,0
	scaffold,evalue,frame = '','',''
	querylist,sbjlist = [],[]
	for line in blast:
		k = line.split()
		if len(k) == 0: #checks if line is not empty
			continue
		if line[:7] == "Length=" and length == 0: #it only takes length 1x (this is query length)
			length = line.strip()[7:]
			out.write("Querylength:,%s\n\n" %length)
			out.write("Scaffold,e-val,Frame,Query,,Subject,,Webapollo\n")
		elif k[0] == ">": # indicates the start of a new scaffold
			if len(querylist) >= 1: #this condition is not met at the first scaffold
				querylist,sbjlist = saveinfo(out,scaffold,evalue,frame,querylist,sbjlist)
			score = 0
			scaffold = k[1]
		elif k[0] == "Score" and scaffold != '':
			if score >= 1: #this condition is not met when this is the first hit on that scaffold
				querylist,sbjlist = saveinfo(out,scaffold,evalue,frame,querylist,sbjlist)
			evalue = k[7][:-1]
			score += 1
		elif k[0] == "Frame":
			frame = k[2]
		elif k[0] == "Query":
			querylist.append(int(k[1]))
			querylist.append(int(k[3]))
		elif k[0] == "Sbjct":
			sbjlist.append(int(k[1]))
			sbjlist.append(int(k[3]))
		elif k[0] == "Lambda":
			saveinfo(out,scaffold,evalue,frame,querylist,sbjlist)
	out.close()
	blast.close()

'''
Run the blastreader for each file in the folder:
'''
for filename in os.listdir(infolder):
	if filename[0] == ".": #So not to run on hidden files...
		print "Skipping hidden file:", filename
		continue
	elif filename[-4:] != ".txt":
		continue
	blastreader(filename)
