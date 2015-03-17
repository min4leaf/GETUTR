import os
import sys
from subprocess import Popen, PIPE, STDOUT

inputFile = sys.argv[1]

def readbam(inputFile):
	print "SAM"
	maxLeng = 0
	bed = dict()
	commandBase = 'samtools view ' + inputFile
	output=''
	if True:
		p=Popen(commandBase, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
		output = p.stdout.readlines()
	for lines in output:
		oneline = lines.strip('\n').split('\t')
		chr = oneline[2]
		if oneline[1]=='0':	## 0:+, 16:-
			sense = '+'
		else:
			sense = '-'
		start = int(oneline[3])
		end = start
		if oneline[5][-1]=='M':
			end += int(oneline[5][:-1])
		#else:
		#	print lines
		#	sys.exit(0)
		tag = oneline[0]
		qual = oneline[4]
		key = int(start/10000)
                if not bed.has_key(chr):
                        #bed[chr] = []
                        bed[chr] = dict()
                if not bed[chr].has_key(key):
                        bed[chr][key] = []
                bed[chr][key].append([start, end, sense, qual, tag])
		if maxLeng<end-start:
			maxLeng = end-start
	return bed, maxLeng

def bamtobed(inputFile):
	print "BAMBED"
	maxLeng = 0
	#bed = {'chr1':[], 'chr2':[], 'chr3':[], 'chr4':[], 'chr5':[], 'chr6':[], 'chr7':[], 'chr8':[], 'chr9':[], 'chr10':[], 'chr11':[], 'chr12':[], 'chr13':[], 'chr14':[], 'chr15':[], 'chr16':[], 'chr17':[], 'chr18':[], 'chr19':[], 'chr20':[], 'chr21':[], 'chr22':[], 'chrX':[], 'chrY':[], 'chrM':[], 'etc':[]}
	bed = dict()
	commandBase = 'bedtools bamtobed -i ' + inputFile
	output=''
	if True:
		p=Popen(commandBase, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
		output = p.stdout.readlines()
	for lines in output:
		oneline = lines.strip('\n').split('\t')
		chr = oneline[0]
		sense = oneline[5]
		start = int(oneline[1])+1
		end = int(oneline[2])+1
		tag = oneline[3]
		qual = int(oneline[4])
		key = int(start/10000)
		if not bed.has_key(chr):
			#bed[chr] = []
			bed[chr] = dict()
		if not bed[chr].has_key(key):
			bed[chr][key] = []
		bed[chr][key].append([start, end, sense, qual, tag])
		if maxLeng<end-start:
			maxLeng = end-start
		#print 'bed:'+str(start)+' '+str(end)
#	for chr in bed.keys():
#		bed[chr].sort()
	return bed, maxLeng


def readbed(inpufFile):
	print "BED"
	maxLeng = 0
	#bed = {'chr1':[], 'chr2':[], 'chr3':[], 'chr4':[], 'chr5':[], 'chr6':[], 'chr7':[], 'chr8':[], 'chr9':[], 'chr10':[], 'chr11':[], 'chr12':[], 'chr13':[], 'chr14':[], 'chr15':[], 'chr16':[], 'chr17':[], 'chr18':[], 'chr19':[], 'chr20':[], 'chr21':[], 'chr22':[], 'chrX':[], 'chrY':[], 'chrM':[], 'etc':[]}
	bed = dict()
	fin = fopen(inputFile, "r")
	for lines in fin:
		oneline = line.split('\t')
		chr = oneline[0]
		sense = oneline[5]
		start = int(oneline[1])+1
		end = int(oneline[2])+1
		tag = oneline[3]
		qual = int(oneline[4])
		key = int(start/10000)
		if not bed.has_key(chr):
			bed[chr] = dict()
		if not bed[chr].has_key(key):
			bed[chr][key] = []
		bed[chr][key].append([start, end, sense, qual, tag])
		if maxLeng<end-start:
			maxLeng = end-start
#	for chr in bed.keys():
#		bed[chr].sort()
	return bed, maxLeng
