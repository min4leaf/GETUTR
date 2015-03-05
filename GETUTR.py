from rpy import *
from types import *
import numpy
import random
import annotationImporter_min as Anno
import sys

import filemanager as FM
import RNAseq_smoothing as smooth
import cleavageDetection as cDetec


outputSmoothed = 1
outputCleavage = 1



##### Input Arguments and PreProcessing #####
method = 10
fileCore = 'default'
inputBED = 'default.bed'
fileCore = 'default'

ptr = 1
while ptr<len(sys.argv) and len(sys.argv)>2:
	if sys.argv[ptr]=='-i':
		ptr += 1
		inputFile = sys.argv[ptr]
	elif sys.argv[ptr]=='-o':
		ptr += 1
		fileCore = sys.argv[ptr]
	elif sys.argv[ptr]=='-m':
		ptr += 1
		method = int(sys.argv[ptr])
	else:
		ptr += 1	
	ptr += 1

if method==10:
	methods = "PAVA"
elif method==0:
	methods = "Max.fit"
elif method==1:
	methods = "Min.fit"
else:
	methods = "PAVA"

if inputFile[-3:]=='bam' or inputFile[-3:]=='BAM':
	if fileCore=='default':
		fileCore = inputFile[:-3]
	inputBED = FM.bamtobed(inputFile)
	#inputBED = FM.readbam(inputFile)
elif inputFile[-3:]=='bed' or inputFile[-3:]=='BED':
	if fileCore=='default':
		fileCore = inputFile[:-3]
	inputBED = FM.readbed(inputFile)
else:
	if fileCore=='default':
		fileCore = inputFile
	inputBED = FM.readbed(inputFile)
#sys.exit(0)

##### TEST #####
#for key in inputBED.keys():
#	for i in inputBED[key]:
#		print i[0]
#sys.exit(0)
#for key in inputBED.keys():
#	for key2 in inputBED[key]:
#		for i in inputBED[key][key2]:
#			print key2, i[0]
#sys.exit(0)
estimatedUTR1 = smooth.smoothing(inputBED, methods)

def writeUTR(UTR, fileCore):
	#track type=wiggle_0 name="HeLa_124_ratio_+_max" description="HeLa_124_ratio_+_max"
	#chrY    150800  150801  1.33439188984	
	fout = open(fileCore, 'w')
	geneAnno = Anno.getRefFlatUniq('hg19')
	for geneid in UTR.keys():
		gene = geneAnno[geneid]
		chrom = gene.chr()
		sense = gene.sense()
		usages = UTR[geneid]
		if len(usages)>0:
			fout.write('track type=wiggle_0 name=\"' + geneid + '\" description=\"' + geneid + '\"\n')
			for usage in usages:
				maxv = max(usage[1:3])
				minv = min(usage[1:3])
				usage[1] = minv
				usage[2] = maxv
				usage = map(str,usage)
				fout.write('\t'.join(usage) + '\n')
				#fout.write(chrom + '\t' + '\t'.join(usage) + '\n')
				##print chrom, '\t'.join(usage)
	fout.close()

if outputSmoothed == 1:
	writeUTR(estimatedUTR1, fileCore+"."+methods+".smoothed.bed")
print 'done'
#sys.exit(0)

#Traceback (most recent call last):
#  File "GETUTR.py", line 94, in <module>
#    writeUTR(estimatedUTR1, fileCore+"."+methods+".smoothed.bed")
#  File "GETUTR.py", line 76, in writeUTR
#    for geneid in UTR.keys():
#AttributeError: 'tuple' object has no attribute 'keys'

###################################



PCS = cDetec.cleavDetect(estimatedUTR1)
estimatedUTR2, PCS = cDetec.makeSeq(estimatedUTR1, PCS)
print 'done'
#sys.exit(0)

if outputCleavage == 1:
	fout = open(fileCore+"."+methods+"_PCS", "w")
	for geneid in PCS.keys():
		for item in PCS[geneid]:
			fout.write(str(item[0])+"\t"+str(item[1])+"\t"+str(item[2])+"\t"+str(item[3])+"\n")
	fout.close()
	
	fout = open(fileCore+"."+methods+"_UTR", "w")
	for geneid in estimatedUTR2.keys():
		for item in estimatedUTR2[geneid]:
			fout.write(str(item[0])+"\t"+str(item[1])+"\t"+str(item[2])+"\t"+str(item[3])+"\t"+str(item[4])+"\n")
	fout.close()

