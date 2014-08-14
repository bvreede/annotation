import sys

if len(sys.argv) <= 1:
	sys.exit("USAGE: python readblast_mini.py path/to/blast.output.file")

blast = open(sys.argv[1])
outname = sys.argv[1] + "_smy"
out = open(outname, "w")

def saveinfo(scaffold,evalue,frame,querylist,sbjlist):
	qstart,qend = min(querylist),max(querylist)
	sstart,send = min(sbjlist),max(sbjlist)
	wa = scaffold + "%3A"
	for i in range(14-len(scaffold)):
		scaffold += ' '
	out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t\t%s%s..%s\n" %(scaffold,evalue,frame,qstart,qend,sstart,send,wa,sstart,send))
	querylist = []
	sbjlist = []
	return querylist,sbjlist


length,score = 0,0
scaffold,evalue,frame = '','',''
querylist,sbjlist = [],[]
for line in blast:
	k = line.split()
	if len(k) == 0: #checks if line is not empty
		continue
	if line[:7] == "Length=" and length == 0: #it only takes length 1x (this is query length)
		length = line.strip()[7:]
		out.write("Total length query: %s\n\n" %length)
		out.write("Scaffold\te-val\tFrame\tQuery\t\tSubject\t\t\tWebapollo\n")
	elif k[0] == ">": # indicates the start of a new scaffold
		if len(querylist) >= 1: #this condition is not met at the first scaffold
			querylist,sbjlist = saveinfo(scaffold,evalue,frame,querylist,sbjlist)
		score = 0
		scaffold = k[1]
	elif k[0] == "Score" and scaffold != '':
		if score >= 1: #this condition is not met when this is the first hit on that scaffold
			querylist,sbjlist = saveinfo(scaffold,evalue,frame,querylist,sbjlist)
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
		saveinfo(scaffold,evalue,frame,querylist,sbjlist)
