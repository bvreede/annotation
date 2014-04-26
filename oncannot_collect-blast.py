import csv, urllib2, re
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

inputdb = "segmentation_short.csv"

genelist = csv.reader(open(inputdb))
output = open("%s-output.csv" %(inputdb[:-4]),"w")

def reverser(CDSlist):
	# TODO should put gene entry in reverse complement.
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
	for line in region_xml:
		linecheck = 0 #to ensure that lines are not added in duplicate
		if line.find('#E0FFFF">')!=-1:
			if exoff == 1:
				if linecheck == 0:
					pre_exon += line.strip()
					linecheck = 1
				exons.append(pre_exon)
				pre_exon = ""
				exoff = 0
		if line.find('#0000F')!=-1:
			if linecheck == 0:
				pre_exon += line.strip()
				linecheck = 1
			linecheck = 1
			exoff = 1
		if exoff == 1:
			if linecheck == 0:
				pre_exon += line.strip()
				linecheck = 1
	if len(exons)!=0: #prevents crash in case the gene gives an error
		CDSlist = [] #will collect the cleaned up exons
		for each in exons: 
			exon = re.split('</span><span style="Background-color: #E0FFFF; Color: #0000F.">', each) #removes html info
			exon2 = "".join(exon[1:]).split('</span><span style="Background-color: #E0FFFF">') # removes intron-before-exon, joins the list, splits off the intron-after-exon
			CDSlist.append(exon2[0])		
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
		output.write("%s,%s,whole,%s\n" %(fbid,genename,totalgene)) #entry with whole genome info
		for i in range(len(CDSlist)):
			output.write("%s,%s,exon%s,%s\n" %(fbid,genename,i+1,CDSlist[i]))
		return CDSlist[0][:3] #returns first codon; if not ATG then gene entry is reverse complement
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

