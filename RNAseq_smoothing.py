from rpy import *
import sys
import math
import random
import regression as RGR
import utilityModule_min as util
import annotationImporter_min as Anno


def get_geneIntervals(geneAnno):
	geneChrD=dict(); geneAnno2=dict()
	for gene in geneAnno.values():
		if not(geneChrD.has_key(gene.chr())): geneChrD[gene.chr()]=[]
		geneChrD[gene.chr()].append([gene.txLocus().start(), gene.txLocus().end(), gene])

	geneInfo=dict()
	chromosomes = geneChrD.keys()
	def cmp0(a1, a2): return cmp(a1[0], a2[0])
	for chrom in chromosomes:
		genes = geneChrD.get(chrom,[])
		genes.sort(cmp0)
		geneNum= len(genes)
		ngenes=[]
		for gene in genes: geneInfo[gene[2].geneID()]= [gene[2].txLocus().chr(), gene[2].sense()] 
		if geneNum>0:
			predistance=2000
			preEnd=[0]
			for g in genes[:-1]:
				j=0
				while (genes[j][0] < g[1] or genes[j][2].sense()!=g[2].sense()):
					j+=1
					if j>geneNum-3: break
				if j>geneNum-1: break
				ng = genes[j]
				distance = ng[0]-g[1]+1
				if distance<0: distance=2000
				if g[2].sense()=='+': g.append(distance)
				else:
					preEnd.sort()
					for i in xrange(len(preEnd)):
						if preEnd[i]>g[0]: break
					predistance = g[0]-preEnd[i-1]+1
					g.append(predistance)

				ngenes.append(g)
				if g[2].sense()=='-':  preEnd.append(g[1])
				geneAnno2[g[2].geneID()]=[g[2], g[-1]]
	return geneAnno2

#def senseDetection(list):
	# input::
	# list list = [string chr, int position, int value]
	# output::
	# string '+'/'-'
#	sense = ''
#	xVal = []
#	yVal = []
#	for i in xrange(len(list)):
#		xVal.append(list[i][1])
#		yVal.append(list[i][2])
#	set_default_mode(NO_CONVERSION)
#	linear_model = r.lm(r("y ~ x"), data = r.data_frame(x=xVal, y=yVal))
#	set_default_mode(BASIC_CONVERSION)
#	gradient = linear_model.as_py()['coefficients']['x']
#	if gradient>0:
#		sense = '-'
#	else:
#		sense = '+'
#	return sense

def simple1(utrSignalRaw, utrPosRaw, chr, sense):
	# input::
	# list utrSignalRaw = [float value]
	# list utrPosRaw =  [int value]
	# sense = '+'/'-'
	# output::
	# list utrValue = [string chr, int start, int end, float value]
	utrValue = []
	if len(utrSignalRaw)>0 and len(utrSignalRaw)==len(utrPosRaw):
		utrSignal = []
		utrPos = []
		utrSignal.append(utrSignalRaw[0])
		utrPos.append(utrPosRaw[0])
		count = 1.0
		for i in xrange(1,len(utrSignalRaw)-1):
			if utrPosRaw[i]==utrPosRaw[i-1]:
				count += 1.0
				utrSignal[-1] += utrSignalRaw[i]
			else:
				utrSignal[-1] /= count
				count = 1.0
				utrSignal.append(utrSignalRaw[i])
				utrPos.append(utrPosRaw[i])
		if utrPosRaw[len(utrPosRaw)-1]==utrPosRaw[len(utrPosRaw)-2]:
			count += 1.0
			utrSignal[-1] += utrSignalRaw[len(utrPosRaw)-1]
			utrSignal[-1] /= count
		else:
			utrSignal[-1] /= count
			utrSignal.append(utrSignalRaw[len(utrPosRaw)-1])
			utrPos.append(utrPosRaw[len(utrPosRaw)-1])

		start = utrPos[0]
		value = utrSignal[0]
		for i in xrange(len(utrSignal)-1):
			if utrSignal[i]!=value:
				end = utrPos[i]
				utrValue.append([chr, start, end, value, sense])
				start = end
				value = utrSignal[i]
		if utrSignal[len(utrSignal)-1]!=value:
			end = utrPos[len(utrSignal)-1]
			utrValue.append([chr, start, end, value, sense])
			start = end
		utrValue.append([chr, start, utrPos[len(utrSignal)-1]+1, utrSignal[len(utrSignal)-1], sense])
	return utrValue
	
def simple2(utrValue):
	# input::
	# list utrValue = [string chr, int start, int end, float value, sense]
	# output::
	# list utrSValue = [string chr, int start, int end, float value, sense]
	utrSValue = []
	if len(utrValue)>0:
		chr = utrValue[0][0]
		start = utrValue[0][1]
		value = utrValue[0][3]
		sense = utrValue[0][4]
		for i in xrange(len(utrValue)-1):
			if utrValue[i][3]!=value:
				end = utrValue[i][1]
				utrSValue.append([chr, start, end, value, sense])
				start = end
				value = utrValue[i][3]
		if utrValue[len(utrValue)-1][3]!=value:
			end = utrValue[len(utrValue)-1][1]
			utrSValue.append([chr, start, end, value, sense])
			start = end
		utrSValue.append([chr, start, utrValue[len(utrValue)-1][2], utrValue[len(utrValue)-1][3], sense])
	return utrSValue

def preprocessing():
	geneAnno1st = Anno.getRefFlatUniq('hg19')
	geneAnno2nd = get_geneIntervals(geneAnno1st)
	return geneAnno2nd

def len_dict(dic):
	length = 0
	for key in dic.keys():
		length += len(dic[key])
	return length

def getPos(RNA, locus):
	#
	#
	ret = []
	chr = locus.chr()
	sense = locus.sense()
	start = locus.start()
	end = locus.end()
	startKey = int(start/10000)
	endKey = int(end/10000)
	if RNA.has_key(chr):
		for i in xrange(startKey, endKey+1):
			if RNA[chr].has_key(i):
				for rna in RNA[chr][i]:
					if rna[0]<=end and start<=rna[1]:
						#print rna[0], rna[1], rna[2], rna[3], rna[4]
						if sense==rna[2]:
							ret.append(rna)
	return ret

def add(list, sub):
	for i in xrange(len(sub)):
		list.append(sub[i])
	return list

def smoothing(inputBED, method='PAVA'):
	# input::
	# list inputBED = [start, end, sense, qual, tag]
	# string method 
	# output::
	# list ret_list = [string chr, string sense, int start_pos, int end_pos, double value]

	usageRatioD = dict()
	tUTRLoci = dict()

	if len_dict(inputBED)==0:
		return usageRatioD
	
	geneAnno = preprocessing()
	
	
	for gene, interval in geneAnno.values():
		geneid = gene.geneID()
		print geneid
		#chr = gene.chr()
		if geneid[:2]!='NR':
			interval = min(interval, 3000)
			beforeUTR = gene.fpExons(gene.sense()) + gene.cdExons()
			beforeUTRLen = sum(map(lambda x: x.len(), beforeUTR))
			txExon = gene.txExons()
			if gene.sense()=='+':
				txExon[-1] = util.Locus(gene.chr(), txExon[-1].start(), txExon[-1].end()+interval, '+')
			else:
				txExon[0] = util.Locus(gene.chr(), txExon[0].start()-interval, txExon[0].end(), '-')
			print geneid
			usageRatio = []
			#if geneid=='NM_001038633':
			#	sys.exit(0)
			#print txExon[0].start(), txExon[-1].end()
			usageRatio = sAlgorithm(txExon, inputBED, beforeUTRLen+1, method)
			#	sys.exit(0)
			usageRatioD[geneid] = usageRatio

	return usageRatioD

def sAlgorithm(txExons, inputBED, utrstart=0, method='PAVA', density_kernel='gaussian'):
	res=[]; resTemp=[]; resPosA=[]; mpos=dict(); mpos2=dict()
	chr = txExons[0].chr()
	sense = txExons[0].sense()
	if sense=='+':
		rp=0
		for tx in txExons:
			print tx
			#res += SM.getBam(bamfile, locus, 'sp')#; print locus,
			Pos = getPos(inputBED, tx)
			resTemp = add(res, Pos)
			res = resTemp
			for i in range(tx.start(), tx.end()+1):
				mpos[i]=rp
				mpos2[rp]=i
				rp+=1
	elif sense=='-':
		rp = sum(map(lambda x: x.len(), txExons))-1
		for tx in txExons: 
			print tx
			#res += SM.getBam(bamfile, locus, 'sp')#; print locus,
			Pos = getPos(inputBED, tx)
			resTemp = add(res, Pos)
			res = resTemp
			for i in range(tx.start(), tx.end()+1):
				mpos[i]=rp
				mpos2[rp]=i
				rp-=1

	if sense=='+':
		utrStart = txExons[0].start()
	else:
		utrStart = txExons[-1].end()
	#print "length"+str(len(res))
	for re in res:
		relPos = -1
		if sense=='+':
			relPos = mpos.get(int(re[0]), -1)
			#print str(re[0])+"::"+str(relPos)
		else:
			relPos = mpos.get(int(re[0])+35, -1)
			#print str(re[0]+35)+" "+str(relPos)
		#print re[0]
		if relPos>=0:
			resPosA.append(int(relPos))
	#utrLocus = util.Locus(txExons[0].chr(), txExons[0].start(), txExons[-1].end(), txExons[0].sense())
	usageRatio=[]
   	utrLen = sum(map(lambda x: x.len(), txExons))
	rpk = 1000.0*len(resPosA)/float(utrLen)
	# print min(mpos.values()), max(mpos.values())
	newUsages=[]
	#print rpk, len(resPosA)
	if rpk>=1 and len(resPosA)>1:
		smoothedValue = r.density(resPosA, kernel=density_kernel, adjust=1.0)
		if smoothedValue.has_key('x'):
			#try:
			#print txExons[0].chr()
			usageRatio, usageRaw = makeSmoothed(smoothedValue, utrstart, utrLen, chr, sense, method)# method='maxStart' or 'neverLower'
			for usage in usageRatio:
				if mpos2.has_key(usage[1]) and mpos2.has_key(usage[2]):
					usage[1] = mpos2[usage[1]]
					usage[2] = mpos2[usage[2]]
					newUsages.append(usage)

	return newUsages



def makeSmoothed(smoothedValue, utrstart, utrLen, chr, sense, method):
	arrayStart = len(filter(lambda x: x<utrstart, smoothedValue['x']))
	utrSignal = map(float, smoothedValue['y'][arrayStart:])
	utrPos = map(int, smoothedValue['x'][arrayStart:])
	#print 'utrSignal'
	#for i in xrange(len(utrSignal)):
	#	print utrPos[i], utrSignal[i]
	utrRaw = simple1(utrSignal, utrPos, chr, sense)
	#print 'utrRaw'
	#for i in xrange(len(utrRaw)):
	#	print utrRaw[i]
	utrValue = []
	if len(utrRaw)>0:
		maxValue = 0
		maxPtr = 0
		for i in xrange(len(utrRaw)):
			if utrRaw[i][3]>maxValue:
				maxValue = utrRaw[i][3]
				maxPtr = i
			if sense=='-':
				temp = utrRaw[i][1]
				utrRaw[i][1] = utrRaw[i][2]
				utrRaw[i][2] = temp

		if method=="PAVA":
			utrRaw.reverse()
			weight = []
			value = []
			for i in xrange(len(utrRaw)):
				weight.append(abs(utrRaw[i][1]-utrRaw[i][2]))
				value.append(utrRaw[i][3])
			value = RGR.isoReg(value, weight)
			maxValue = max(value)
			for i in xrange(len(utrRaw)):
				utrRaw[i][3] = value[i]
				utrValue.append([chr, utrRaw[i][1], utrRaw[i][2], float(utrRaw[i][3]/maxValue), utrRaw[i][4]])
			if sense=='+':
				utrRaw.reverse()
				utrValue.reverse()
		else:
			ptr = 0
			while ptr<maxPtr:	
#				print utrRaw[ptr][1], utrRaw[ptr][2], utrRaw[ptr][3], utrRaw[ptr][4]
				utrRaw[ptr][3] = maxValue
				ptr += 1
			
			if method=="Max.fit":	#neverLower
				utrRaw.reverse()
				criterion = utrRaw[0][3]
				for i in xrange(len(utrRaw)):
					criterion = max(criterion, utrRaw[i][3])
					utrRaw[i][3] = criterion
					utrValue.append([chr, utrRaw[i][1], utrRaw[i][2], float(utrRaw[i][3]/maxValue), utrRaw[i][4]])
				if sense=='+':
					utrRaw.reverse()
					utrValue.reverse()
			elif method=="Min.fit":	#neverUpper
				criterion = utrRaw[0][3]
				for i in xrange(len(utrRaw)):
					criterion = min(criterion, utrRaw[i][3])
					utrRaw[i][3] = criterion
					utrValue.append([chr, utrRaw[i][1], utrRaw[i][2], float(utrRaw[i][3]/maxValue), utrRaw[i][4]])
				if sense=='+':
					utrRaw.reverse()
					utrValue.reverse()
				utrRaw.reverse()
				utrValue.reverse()
		if sense=='+':
			utrRaw[0][1] = utrstart-1
			utrValue[0][1] = utrstart
		else:
			utrRaw[-1][2] = utrstart+1
			utrValue[-1][2] = utrstart+1


	return simple2(utrValue), simple2(utrRaw)
		
		
