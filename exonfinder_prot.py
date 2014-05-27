import csv, urllib2, re, sys, os
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


if len(sys.argv) <= 1:
	sys.exit("USAGE: python exonfinder_prot.py path/to/inputfile (list of flybase IDs and genenames in csv)")

inputdb = sys.argv[1] # input file
inputname = inputdb.split("/")[-1]

csvfolder = "csv"
if os.path.exists(csvfolder):
	pass
else:
	os.mkdir(csvfolder)

outputfolder = "results"
if os.path.exists(outputfolder):
	pass
else:
	os.mkdir(outputfolder)

genelist = csv.reader(open(inputdb))
output = open("%s/%s-exons.csv" %(csvfolder,inputname[:-4]),"w")

def metaextract(header,genemeta):
	'''
	Takes the header of each isoform and extracts meta-information
	(i.e. which isoform uses which exons)
	Saves that information in the metafile.
	'''
	protID = header.split()[0][1:]
	tag = re.compile('loc=[1-4|XRL]{,2}:[a-z|0-9|\(|\.|,]*').search(header)
	if tag == None:
		exons = '(exons)'
	else:	
		tag1 = tag.group() #parses out the exon locations used for this isoform
		tag2 = re.split('loc=[1-4|XRL]{,2}[a-z|\(|:]*', tag1) #splits off the useless info
		tag3 = "".join(tag2)
		exons = tag3.replace(',','\n')
	tag4 = re.compile('name=[A-Z|a-z|\-|0-9]*').search(header)
	if tag4 == None:
		name = ''
	else:
		tag5 = tag4.group() # parses out isoform name
		name = tag5[5:]
	genemeta.write("ProtID = %s\nIsoform = %s\n%s\n\n" %(protID,name,exons))

def prot_exons(isolist):
	'''
	Takes the list of exon locations per isoform ('nnnnnn..nnnnnn')
	and finds the corresponding start and end sites in the protein sequence
	of that isoform.
	Returns those locations as integers in a list that is twice as long
	as the initial list: per entry in 'isolist' there is a 'start' and an 'end'
	entry in the 'isolist_p'.
	'''
	isonum = len(isolist)
	isolist_p = []
	isolist_g = []
	texd_i=0
	for i in range(isonum):
		ex_i = isolist[i].split('..')
		if len(ex_i) <2: # if a non-exon is in the list it will not be used.
			continue
		lexd_i = int(ex_i[1])-int(ex_i[0])+1 #the length of the exon in basepairs
		sexp_i = int(texd_i/3+0.1) #generates aminoacid start location: ensures 1.0=1; 1.3=1; 1.6=1; 1.99=2
		isolist_p.append(sexp_i) #appends start location of exon to list
		texd_i += lexd_i
		eexp_i = int(texd_i/3-0.1+1) #generates aminoacid end location: ensures 1.0=0; 1.3=1; 1.6=1; 1.99=1; adds 1 to INCLUDE that amino acid (it is the last amino acid, not the end site!)
		isolist_p.append(eexp_i) #appends end location of exon to list
	return isolist_p
		
def prot_dict(isolist,isolist_p,isoseq,dirx):
	'''
	Takes the list of exon locations per isoform ('nnnnnn..nnnnnn') and uses it as keys;
	uses the isolist_p list (made in the prot_exons module) with start and end sites to
	find the corresponding protein sequence of that isoform.
	Returns a dictionary with exon locations as keys, and the sequence as values.
	'''
	isodict = {}
	if dirx == '-':
		isoseq = isoseq[::-1]
	for i in range(len(isolist)):
		exon = isolist[i] #identifier of the exon is the 'nnnnn..nnnn' genomic location
		start = isolist_p[i*2] #start site is found in the protein index list isolist_p; there are twice as many items in this list (start, end, start, end) so identifier has to be multiplied
		end=isolist_p[i*2+1]
		sequence = isoseq[start:end] #gets the right string for this exon from the total sequence
		if dirx == '-':
			sequence = sequence[::-1]
		isodict[exon] = sequence #saves in dictionary: the exon-specific sequence with the exon identifier as key
	return isodict

"""
def isolator(isolist):
	'''
	takes the total list of exons from all isoforms of a gene (produced by prot_exons)
	and finds the duplicates and overlaps. Returns a set of keys
	that can be removed from the exon dictionary.
	'''
	isoset=set(isolist) #puts all items in a set. This deletes duplicates.
	isolist2 = []
	for i in isoset:
		spliti = i.split('..') #splits start and end location
		if len(spliti) == 2:
			ilist = [int(spliti[0]),int(spliti[1])] #makes a new sublist with start and end locations
			isolist2.append(ilist)
	isolist2.sort() #sorts on start location (first item in the sublists)
	if len(isolist2) > 0:
		check = isolist2[0] #the sublist (start,end) against which each test entry will be checked
	pop = []
	for n in range(1,len(isolist2)):	#go through each entry after the first (the first is 'check', the others are 'test')
		if isolist2[n][0] == check[0]: 		#if the start sequences are the same...
			if isolist2[n][1] >= check[1]:		#...and the end of the test is after the end of the check
				pop.append(check)			#then discard the check
				check = isolist2[n]			#and make the test entry a new check
			else:					#...and the end of the check is after the end of the test
				pop.append(isolist2[n])			#discard the test entry
		else:					#if the start sequences are not equal (which means test starts after check)...
			if isolist2[n][1] <= check[1]:		#...and the end of the test is before or equal to the end of the check
				pop.append(isolist2[n])			#then discard the test entry
			else:					#...and the end of the test is after the end of the check
				check = isolist2[n]			#make the test entry a new check
	return pop
"""

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

def proteinscan(fbid,chrom,dirx,genename):
	'''
	Takes a gene's flybase ID, chromosome and orientation and goes
	to its translation fasta page.
	Calls prot_exons and prot_dict, and isolator to remove duplicates.
	Returns dictionary of individual exons and their sequence.
	'''
	transl_url = "http://flybase.org/cgi-bin/getseq.html?source=dmel&id=%s&chr=%s&dump=PrecompiledFasta&targetset=translation" %(fbid,chrom)
	transl_xml = urllib2.urlopen(transl_url)
	genemeta = open("%s/%s-meta.txt" %(outputfolder,genename),"w")
	genemeta.write("++++ISOFORMS:++++\n\n")
	isoseq = ""
	isolist = []
	isolist_collect = []
	genedict = {}
	# finding the headers and the sequences
	for line in transl_xml:
		if line[0]== ">": # this line has the fasta header
			metaextract(line,genemeta)
			isolist_p = prot_exons(isolist) #get corresponding protein start-end sites from chromosomal start-end sites
			if len(isolist_p) != 0:
				isodict = prot_dict(isolist,isolist_p,isoseq,dirx) #put all in dictionary
				genedict.update(isodict)
			tag = re.compile('loc=[1-4|XRL]{,2}:[a-z|0-9|\(|\.|,]*').search(line)
			if tag != None:
				tag1 = tag.group() #parses out the exon locations used for this isoform
			else:
				genemeta.write("error finding exon location in line:\n%s\n\n" %(line))
				tag1 = 'err'
			tag2 = re.split('loc=[1-4|XRL]{,2}[a-z|\(|:]*', tag1) #splits off the useless info
			tag3 = "".join(tag2)
			isolist = tag3.split(',') #creates list of exons per isoform
			isolist_collect.extend(isolist)
			isoseq = ""
		#collect protein sequences per isoform
		else: #lines with sequence data
			isoseq += line.strip()
	#what follows now is a repeat so that the last exon is included
	isolist_p = prot_exons(isolist) #get corresponding protein start-end sites from chromosomal start-end sites
	if len(isolist_p) != 0:
		isodict = prot_dict(isolist,isolist_p,isoseq,dirx) #put all in dictionary
		genedict.update(isodict)
	#end of repeat code
	"""
	pop = isolator(isolist_collect)
	for k in pop:
		i = [str(k[0]),str(k[1])] #make string to enable joining
		popkey = "..".join(i)
		del genedict[popkey]
	"""
	genemeta.close()
	return genedict

for gene in genelist:
	fbid = gene[0]
	genename = gene[1].replace(' ','_')
	chromloc=locfinder(fbid)
	chrom = chromloc.split(':')[0]
	if chrom == "err":
		print "Error finding flybase info gene '%s' at http://flybase.org/reports/%s.html.\nContinuing..." %(genename, fbid)
		continue
	dirx = chromloc[-2]
	genedict = proteinscan(fbid,chrom,dirx,genename)
	exons = 0
	for exon in genedict:
		output.write("%s,%s,%s,%s\n" %(fbid,genename,exon,genedict[exon]))
		exons += 1
	print "For gene '%s': %s exons saved." %(genename,exons)

output.close()
