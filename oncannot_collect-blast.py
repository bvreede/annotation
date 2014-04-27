import csv, urllib2, re
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

inputdb = "segmentation_short.csv"

genelist = csv.reader(open(inputdb))
output = open("%s-output.csv" %(inputdb[:-4]),"w")

def prot_exons(isolist):
	isonum = len(isolist)
	isolist_p = []
	texp_i=0
	for i in range(isonum):
		ex_i = isolist[i].split('..')
		lexp_i = (int(ex_i[1])-int(ex_i[0]))/3 #the length of the exon in aminoacids
		#isolist_p.append(texp_i) #appends start location of exon to list
		texp_i += lexp_i
		isolist_p.append(texp_i) #appends end location of exon to list
	return isolist_p
		
	

def reverser(CDSlist):
	print "reversing sequence..."
	CDSlist_rev = CDSlist[::-1] #reverses the entries in the CDSlist
	CDSlist = [] #empties CDSlist so reverse complements can append
	for exon in CDSlist_rev:
		seqobj = Seq(exon,IUPAC.unambiguous_dna) # puts sequence in a sequence object
		seqinv = str(seqobj.reverse_complement())
		CDSlist.append(seqinv)
	return CDSlist

def exonfinder(fbid,genename): # collects individual exons from a flybase gene entry
	# open flybase page with corresponding ID: genome region
	print "reading FlyBase entry for %s..." %genename	
	region_url = "http://flybase.org/cgi-bin/getseq.html?source=dmel&id=%s&dump=DecoratedFasta" %fbid
	region_xml = urllib2.urlopen(region_url)
	exons = [] #will collect lines that include an exon
	pre_exon = ""
	exoff = 0 #collecting exon: yes (1) or no (0)
	chrom = "err"

	# parse xml document: find CDS by colour
	for line in region_xml:
		if line.find('<title>')!=-1:#collect chromosome data
			if re.compile('<title>[1-3|XRL]{,2}:[0-9]{3}').search(line):
				chrom = re.compile('<title>[1-3|XRL]{,2}:[0-9]{3}').search(line).group()[7:-4]
		if exoff == 1 and line.find('#E0FFFF">')==-1 and line.find('#C0FEFE">')==-1: #identifies a 'middle' line
			pre_exon += line.strip() # add line to exon
		if exoff == 1 and (line.find('#E0FFFF">')!=-1 or line.find('#C0FEFE">')!=-1): #identifies an 'end' line
			pre_exon += line.strip() # add line to exon
			exons.append(pre_exon) # add exon to exon list
			pre_exon = "" # empty exon list
			exoff = 0
		if line.find('#0000F')!=-1 and exoff == 0: #identifies a 'start' line
			pre_exon += line.strip() # add line to exon
			exoff = 1

	# parse exon list: remove html code and introns/excess sequence
	if len(exons)!=0: #prevents crash in case the gene gives an error
		CDSlist = [] # will collect the cleaned up exons
		for each in exons: 
			exon = re.split('</span><span style="Background-color: #[CE]0F[FE]F[FE]; Color: #0000F.">', each) #splits at html info that starts exon, as well as any html info in the middle of the exon
			exon2 = "".join(exon[1:]) # removes intron before exon and turns list into string
			exon3 = re.split('</span><span style="Background-color: #[CE]0F[FE]F[FE]">', exon2) # splits at html info that ends the exon
			CDSlist.append(exon3[0]) # only appends exon to list, excluding any introns that follow	
		totalgene = "".join(CDSlist)
		if totalgene[:3] == "ATG": #checks whether gene has been collected sense or reverse complement
			pass
		else:
			CDSlist = reverser(CDSlist) # or use simply reverser(CDSlist)?
			totalgene = "".join(CDSlist)
			if totalgene[:3] == "ATG": #checks whether gene has been correctly reversed
				pass
			else:
				print "error in reversing sequence for ", genename
		output.write("%s,%s,whole,%s\n" %(fbid,genename,totalgene)) #entry with whole genome info
		for i in range(len(CDSlist)):
			output.write("%s,%s,exon%s,%s\n" %(fbid,genename,i+1,CDSlist[i]))
		CDSlist = [] #empties CDSlist for the next entry
	else:
		output.write("%s,%s,error\n" %(fbid,genename)) #entry when gene gives an error
	return chrom

def proteinscan(fbid,chrom,genename):
	transl_url = "http://flybase.org/cgi-bin/getseq.html?source=dmel&id=%s&chr=%s&dump=PrecompiledFasta&targetset=translation" %(fbid,chrom)
	transl_xml = urllib2.urlopen(transl_url)
	isoseq = ""
	isolist = []
	isoname = ""
	for line in transl_xml:
		#get info on isoforms: which exons are used?
		if line[0]== ">": # this line has the fasta header
			#save previous isoseq somewhere here
			isolist_p = prot_exons(isolist)
			print isolist_p
			#including previous tag information
			tag = re.compile('loc=[1-3|XRL]{,2}:[a-z|0-9|\(|\.|,]*').search(line).group() #parses out the exon locations used for this isoform
			tag2 = re.split('loc=[1-3|XRL]{,2}[a-z|\(|:]*', tag) #splits off the useless info
			tag3 = "".join(tag2)
			isolist = tag3.split(',') #creates list of exons per isoform
			isoname = re.compile('name=[a-z|\-|A-Z]*;').search(line).group()[5:-1] #parses out name
			isoseq = ""
		#collect protein sequences per isoform
		if line[0]!= (">" and "<"): #lines with sequence data
			isoseq += line.strip()
			

for gene in genelist:
	fbid = gene[0]
	genename = gene[1]
	chrom = exonfinder(fbid,genename)
	if chrom != "err":
		proteinscan(fbid,chrom,genename)



### perform blast on all protein sequences ###

