import csv, urllib2, re
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

inputdb = "segmentation_bazooka.csv"

genelist = csv.reader(open(inputdb))
output = open("%s-output.csv" %(inputdb[:-4]),"w")

def reverser(CDSlist):
	print "reversing sequence..."
	CDSlist_rev = CDSlist[::-1] #reverses the entries in the CDSlist
	CDSlist = [] #empties CDSlist so reverse complements can append
	for exon in CDSlist_rev:
		seqobj = Seq(exon,IUPAC.unambiguous_dna)
		seqinv = str(seqobj.reverse_complement())
		CDSlist.append(seqinv)
	return CDSlist

def exonfinder(fbid, genename): # collects individual exons from a flybase gene entry
	# open flybase page with corresponding ID: genome region
	print "reading FlyBase entry for %s..." %genename	
	region_url = "http://flybase.org/cgi-bin/getseq.html?source=dmel&id=%s&dump=DecoratedFasta" %fbid
	region_xml = urllib2.urlopen(region_url)
	exons = [] #will collect lines that include an exon
	pre_exon = ""
	exoff = 0 #collecting exon: yes (1) or no (0)

	# parse xml document: find CDS by colour
	countline = 0
	for line in region_xml:
		countline += 1
		if exoff == 1 and line.find('#E0FFFF">')==-1 and line.find('#C0FEFE">')==-1: #identifies a 'middle' line
			pre_exon += line.strip() # add line to exon
			print "added midline ", countline
		if exoff == 1 and (line.find('#E0FFFF">')!=-1 or line.find('#C0FEFE">')!=-1): #identifies an 'end' line
			pre_exon += line.strip() # add line to exon
			print "added endline ", countline
			exons.append(pre_exon) # add exon to exon list
			pre_exon = "" # empty exon list
			exoff = 0
		if line.find('#0000F')!=-1 and exoff == 0: #identifies a 'start' line
			pre_exon += line.strip() # add line to exon
			print "added startline ", countline
			exoff = 1
	#print exons
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
				print "success"
			else:
				print "error in reversing sequence for " + genename
				CDSlist = reverser(CDSlist)
				#print CDSlist
				print len(exons)
				print len(CDSlist)
		output.write("%s,%s,whole,%s\n" %(fbid,genename,totalgene)) #entry with whole genome info
		for i in range(len(CDSlist)):
			output.write("%s,%s,exon%s,%s\n" %(fbid,genename,i+1,CDSlist[i]))
		#print CDSlist
		CDSlist = [] #empties CDSlist for the next entry
	else:
		output.write("%s,%s,error\n" %(fbid,genename)) #entry when gene gives an error



for gene in genelist:
	fbid = gene[0]
	genename = gene[1]
	exonfinder(fbid,genename)
	


	transl_url = "http://flybase.org/cgi-bin/getseq.html?source=dmel&id=%s&chr=3L&dump=PrecompiledFasta&targetset=translation" %fbid
	transl_xml = urllib2.urlopen(transl_url)


### perform blast on all protein sequences ###

