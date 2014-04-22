import csv, urllib2, re

inputdb = "segmentation_bicoid.csv"

genelist = csv.reader(open(inputdb))
output = open("%s-output.csv" %(inputdb[:-4]),"w")

### collect gene and protein sequences and save ###
for gene in genelist:
	fbid = gene[0]
	genename = gene[1]
	# open flybase page with corresponding ID: genome region
	print "reading FlyBase entry for %s..." %genename	
	region_url = "http://flybase.org/cgi-bin/getseq.html?source=dmel&id=%s&dump=DecoratedFasta" %fbid
	print "URL: %s" %region_url
	region_xml = urllib2.urlopen(region_url)
	exons = []
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
		CDSlist = []
		for each in exons:
			exon = re.split('</span><span style="Background-color: #E0FFFF; Color: #0000F.">', each) #removes html info that signifies exon
			exon2 = "".join(exon[1:]).split('</span><span style="Background-color: #E0FFFF">') # removes intron-before-exon, joins the list, splits off the intron-after-exon
			CDSlist.append(exon2[0])		
		totalgene = "".join(CDSlist)
		print totalgene
		output.write("%s,%s,whole,%s\n" %(fbid,genename,totalgene)) #entry with whole genome info
		for i in range(len(CDSlist)):
			print i
		CDSlist = []
	else:
		output.write("%s,%s,error\n" %(fbid,genename)) #entry when gene gives an error
	


	transl_url = "http://flybase.org/cgi-bin/getseq.html?source=dmel&id=%s&chr=3L&dump=PrecompiledFasta&targetset=translation" %fbid
	transl_xml = urllib2.urlopen(transl_url)


### perform blast on all protein sequences ###

