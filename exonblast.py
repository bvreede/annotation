import csv, os

inputdb = "csv-input/segmentation_short-output.csv"
genome = "/home/barbara/data/genomes/Ofasciatus/Ofas.scaffolds.fa"
blasttype = "tblastn"

now = str(datetime.datetime.now())
date = now[:4] + now[5:7] + now[8:10] + now[11:13] + now[14:16] + now [17:19]
dir_out = "./blast_output_%s" % (date,)
os.mkdir(dir_out)

temp = open('temp.txt','w') # will be used to store exon sequences temporarily in a file
exonlist = csv.reader(open(inputdb))

fbid = ""
for exon in exonlist:
	temp.truncate(0)
	temp.write(exon[3])
	output = "%s/%s_%s.txt" %(dir_out,exon[1],exon[2])
	blast = "%s -db %s -query temp.txt -out %s" %(blasttype,genome,output)
	#os.system(blast)
	if exon[0] == fbid:
		# do the stuff you need to do when still working with the same gene
	else:
		# start a new gene
		fbid = exon[0]

temp.close()
os.remove("temp.txt")



