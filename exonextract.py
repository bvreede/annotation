'''
Very quick script that transforms parts of the meta-files
made by the annotation pipeline, specifically the part where
exons are collected and a list of exons per isoform is given,
to a workable table.
Takes as input: 
* a txt file with the exon list (copy-pasted from the meta file).
* a txt file with the isoforms and their exons (idem).
Spits out a csv file indicating which exon is in which isoform,
but in table format. Can be useful with many isoforms/exons!
'''

exons = open("exons.txt")
isoforms = open("isoforms.txt")
output = open("exontable.csv", "w")

exonlist=[] #collects all the exons in the exonlist file
for line in exons:
	exon = line.split()[0]
	exonlist.append(exon)
	output.write(",%s" %(exon))
exons.close()
output.write("\n")

isoexons = []
for line in isoforms:
	s = line.split()
	if len(s) > 0 and s[0] == "Isoform": #isoform id is collected
		output.write("%s," %(s[2]))
	elif len(s) == 1:
		isoexons.append(s[0].strip())	#exons for this isoform are collected
	elif len(s) == 0: #indicating the end of an isoform
		for e in exonlist:
			if e in isoexons:
				output.write("1,")
			else:
				output.write(",")
		output.write("\n")
		isoexons = []

if len(isoexons) > 0: #is only called if there is no empty final line to process the last exon.
	for e in exonlist:
		if e in isoexons:
			output.write("X,")
		else:
			output.write(",")
			
isoforms.close()
output.close()