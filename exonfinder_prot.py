import csv, urllib2, re
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

inputdb = "csv-input/segmentation_short.csv"

genelist = csv.reader(open(inputdb))
output = open("%s-output.csv" %(inputdb[:-4]),"w")


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
		lexd_i = int(ex_i[1])-int(ex_i[0])+1 #the length of the exon in basepairs
		sexp_i = int(texd_i/3+0.1) #generates aminoacid start location: ensures 1.0=1; 1.3=1; 1.6=1; 1.99=2
		isolist_p.append(sexp_i) #appends start location of exon to list
		texd_i += lexd_i
		eexp_i = int(texd_i/3-0.1+1) #generates aminoacid end location: ensures 1.0=0; 1.3=1; 1.6=1; 1.99=1; adds 1 to INCLUDE that amino acid (it is the last amino acid, not the end site!)
		isolist_p.append(eexp_i) #appends end location of exon to list
	return isolist_p
		
def prot_dict(isolist,isolist_p,isoseq):
	'''
	Takes the list of exon locations per isoform ('nnnnnn..nnnnnn') and uses it as keys;
	uses the isolist_p list (made in the prot_exons module) with start and end sites to
	find the corresponding protein sequence of that isoform.
	Returns a dictionary with exon locations as keys, and the sequence as values.
	'''
	isodict = {}
	for i in range(len(isolist)):
		exon = isolist[i] #identifier of the exon is the 'nnnnn..nnnn' genomic location
		start = isolist_p[i*2] #start site is found in the protein index list isolist_p; there are twice as many items in this list (start, end, start, end) so identifier has to be multiplied
		end=isolist_p[i*2+1]
		sequence = isoseq[start:end] #gets the right string for this exon from the total sequence 
		isodict[exon] = sequence #saves in dictionary: the exon-specific sequence with the exon identifier as key
	return isodict

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
		ilist = [int(spliti[0]),int(spliti[1])] #makes a new sublist with start and end locations
		isolist2.append(ilist)
	isolist2.sort() #sorts on start location (first item in the sublists)
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

def chrom(fbid):
	'''
	Goes to flybase page of gene and returns the sequence
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

def proteinscan(fbid,chrom,genename,revflag):
	transl_url = "http://flybase.org/cgi-bin/getseq.html?source=dmel&id=%s&chr=%s&dump=PrecompiledFasta&targetset=translation" %(fbid,chrom)
	transl_xml = urllib2.urlopen(transl_url)
	isoseq = ""
	isolist = []
	isolist_collect = []
	isoname = ""
	genedict = {}
	for line in transl_xml:
		if line[0]== ">": # this line has the fasta header
			isolist_p = prot_exons(isolist) #get corresponding protein start-end sites from chromosomal start-end sites
			isodict = prot_dict(isolist,isolist_p,isoseq) #put all in dictionary
			genedict.update(isodict)
			tag = re.compile('loc=[1-3|XRL]{,2}:[a-z|0-9|\(|\.|,]*').search(line).group() #parses out the exon locations used for this isoform
			tag2 = re.split('loc=[1-3|XRL]{,2}[a-z|\(|:]*', tag) #splits off the useless info
			tag3 = "".join(tag2)
			isolist = tag3.split(',') #creates list of exons per isoform
			isolist_collect.extend(isolist)
			isoname = re.compile('name=[a-z|\-|A-Z]*;').search(line).group()[5:-1] #parses out name
			isoseq = ""
		#collect protein sequences per isoform
		if line[0]!= ">": #lines with sequence data
			isoseq += line.strip()
	#what follows now is a repeat so that the last exon is included
	isolist_p = prot_exons(isolist)
	isodict = prot_dict(isolist,isolist_p,isoseq) #put all in dictionary
	genedict.update(isodict)
	#end of repeat code
	pop = isolator(isolist_collect)
	for k in pop:
		i = [str(k[0]),str(k[1])] #make string to enable joining
		popkey = "..".join(i)
		del genedict[popkey]
	print genedict


for gene in genelist:
	fbid = gene[0]
	genename = gene[1]
	chromloc=chrom(fbid)
	print genename, chromloc
	





