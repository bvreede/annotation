import csv, urllib2, re, sys, os

if len(sys.argv) <= 1:
	sys.exit("USAGE: python getprot.py path/to/inputfile \n(inputfile is list of flybase IDs and genenames in csv)")

inputdb = sys.argv[1] # input file
inputname = inputdb.split("/")[-1]
outputfolder = inputname.split(".csv")[0]

if os.path.exists(outputfolder):
	pass
else:
	os.mkdir(outputfolder)

genelist = csv.reader(open(inputdb))


def locfinder(fbid):
	'''
	Takes flybase ID of a gene and returns its sequence
	location, which includes the chromosome as well as the
	orientation of the gene.
	'''
	url = "http://flybase.org/reports/%s.html" %(fbid)
	xml = urllib2.urlopen(url)
	seqloc = re.compile('[1-3|XRL]{,2}:[0-9|,]{,15}\.\.[0-9|,]{,15} \[[+-]\]')
	chromloc = "err" # default: prevents error in case no sequence location can be found [yes this happens!].
	for line in xml:
		if seqloc.search(line) != None:			
			chromloc = seqloc.search(line).group()
	return chromloc

def seqparse(fbid,dirx,chrom):
	'''
	Takes a gene's flybase ID, chromosome and orientation and goes
	to its translation fasta page. Returns a dictionary of different
	protein isoforms shown there, with the corresponding sequence as value.
	'''
	transl_url = "http://flybase.org/cgi-bin/getseq.html?source=dmel&id=%s&chr=%s&dump=PrecompiledFasta&targetset=translation" %(fbid,chrom)
	transl_xml = urllib2.urlopen(transl_url)
	seq,protid = "","" #empty before starting
	seqdict = {}
	for line in transl_xml:
		if line[0] == ">":
			if seq != "":
				seqdict[protid] = seq
				seq = ""
			protid = line.split()[0][1:]
		else:
			seq += line.strip()
	seqdict[protid] = seq # save the collected sequence with its ID in a dictionary
	return seqdict


for gene in genelist:
	genename,fbid = gene #each 'gene' contains the full input line, with both the fbid and the name
	genename = genename.replace(" ","_")
	if fbid == "":
		print "Error in document: gene %s has no FBid." %(genename)
		continue
	chromloc = locfinder(fbid)
	if chromloc != "err":
		dirx = chromloc[-2]
		chrom = chromloc.split(":")[0]
		seqdict = seqparse(fbid,dirx,chrom)
	else:
		continue
		print "Error finding sequence location of gene %s. Sequence finding aborted." %(genename)
	for k in seqdict:
		output = open("%s/%s_%s" %(outputfolder,genename,k), "w")
		output.write(seqdict[k])
		output.write("\n")
		output.close()
	print "Saved %s proteins for gene %s." %(len(seqdict),genename)
