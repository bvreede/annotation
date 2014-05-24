import csv, sys, os


if len(sys.argv) <= 1:
	sys.exit("USAGE: python blastalign.py path/to/inputfile (result from exonblast.py)")

inputdb = sys.argv[1] # input file
blastlist = csv.reader(open(inputdb))
outputfolder = "alignments"

output = open("x.txt","w") # decoy outputfile so that close statement can be used

if os.path.exists(outputfolder):
	pass
else:
	os.mkdir(outputfolder)

fbid = ""
for line in blastlist:
	if line[0] != fbid:
		output.close()
		output = open("%s/%s-align.csv" %(outputfolder,fbid),"w")
		# proceed with analysing line
	else:
		# proceed with analysing line

output.close()
os.remove("x.txt") # remove decoy outputfile
