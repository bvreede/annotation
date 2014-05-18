import csv, urllib2, re
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

# currently incomplete code that consists of 
# modules written for oncannot_collect-blast.py
# but removed as easier methods were found.
# contains a translation module, a sequence reverser
# and a module that finds exons in flybase HTML
# of the DNA sequence of a genomic locus.

def translater(seq):
	'''
	takes one DNA sequence and translates all three forward frames
	(but only complete codons)
	returns frame1, frame2, frame3 as translated sequence objects.
	'''
	if len(seq) == 0:
		print "error: empty sequence entered in translater"
		return ["err","err","err"]
	elif len(seq) % 3 == 0:
		seq1 = seq
		seq2 = seq[1:-2]
		seq3 = seq[2:-1]
	elif len(seq) % 3 == 1:
		seq1 = seq[:-1]
		seq2 = seq[1:]
		seq3 = seq[2:-2]
	else:
		seq1 = seq[:-2]
		seq2 = seq[1:-1]
		seq3 = seq[2:]
	trans = {}
	for seq, i in zip((seq1, seq2, seq3), range(3)):
		seqt = Seq(seq, IUPAC.unambiguous_dna)
		trans[i]=seqt.translate()
	return trans[0], trans[1], trans[2]


def reverser(CDSlist):
	'''
	takes list of sequences and reverse-complements them one by one
	returns same list, but backwards, and reverse complement.
	eg.: ['ATTG', 'CCC', 'AAT']  becomes ['ATT', 'GGG', 'CAAT']
	'''
	print "reversing sequence..."
	CDSlist_rev = CDSlist[::-1] #reverses the entries in the CDSlist
	CDSlist = [] #empties CDSlist so reverse complements can append
	for exon in CDSlist_rev:
		seqobj = Seq(exon,IUPAC.unambiguous_dna) # puts sequence in a sequence object
		seqinv = str(seqobj.reverse_complement())
		CDSlist.append(seqinv)
	return CDSlist

def exonfinder(fbid,genename): 
	'''
	collects individual exons from a flybase gene entry
	'''
	# open flybase page with corresponding ID: genome region
	print "reading FlyBase entry for %s..." %genename	
	region_url = "http://flybase.org/cgi-bin/getseq.html?source=dmel&id=%s&dump=DecoratedFasta" %fbid
	region_xml = urllib2.urlopen(region_url)
	exons = [] #will collect lines that include an exon
	pre_exon = ""
	exoff = 0 #collecting exon: yes (1) or no (0)
	chrom = "err"
	revflag = 0 #will be turned on if the gene is reverse complement

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
			revflag = 1
			totalgene = "".join(CDSlist)
			if totalgene[:3] == "ATG": #checks whether gene has been correctly reversed
				pass
			else:
				print "error in reversing sequence for ", genename
		output.write("%s,%s,whole,%s\n" %(fbid,genename,totalgene)) #entry with whole genome info
		for i in range(len(CDSlist)):
			translation=translater(CDSlist[i])
			output.write("%s,%s,exon%s,%s,%s,%s,%s\n" %(fbid,genename,i+1,CDSlist[i],translation[0],translation[1],translation[2]))
		CDSlist = [] #empties CDSlist for the next entry
	else:
		output.write("%s,%s,error\n" %(fbid,genename)) #entry when gene gives an error
	return chrom, revflag
