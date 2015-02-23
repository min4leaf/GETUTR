import os
from subprocess import Popen, PIPE, STDOUT

def samtoolReheader(headerInSam, headerOutBam):
    commandBase = 'samtools view -H ' + headerInSam + '>' + headerInSam + '.header'
    print commandBase
    p=Popen(commandBase, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
    p.wait()
    commandBase = 'samtools reheader ' + headerInSam + '.header' + ' ' + headerOutBam
    print commandBase
    p=Popen(commandBase + ' >' + headerOutBam + '2', shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
    p.wait()
    p = Popen('rm ' + headerInSam +'.header', shell=True); p.wait()
    p = Popen('rm ' + headerOutBam, shell=True); p.wait()
    p = Popen('mv ' + headerOutBam +'2 ' + headerOutBam, shell=True); p.wait()

def sortBed(bedfile, outputfile):
    commandBase = 'bedtools sort -i ' + bedfile + ' >' + outputfile
    print commandBase
    p=Popen(commandBase, shell=True)
    p.wait()
  
def samtoolMerge(bamfilelist, outputfilename):
    commandBase = 'samtools merge -nrf '+ outputfilename + ' ' + ' '.join(bamfilelist)
    print commandBase
    p=Popen(commandBase, shell=True)
    p.wait()
    samtoolSort(outputfilename)
    samtoolIndex(outputfilename)
    samtoolMappingInfo(outputfilename)

def samToBam(samfilename):
    bamfilename = '.'.join(samfilename.split('.')[:-1]) + '.bam'
    commandBase = 'samtools view -Sb ' + samfilename + ' >' + bamfilename
    print commandBase
    p=Popen(commandBase, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
    p.wait()
    os.remove(samfilename)
    
def samtoolSort(outputfilename):
    sortedFile = '.'.join(outputfilename.split('.')[:-1]) + 'Sorted'
    commandBase = 'samtools sort ' + outputfilename + ' ' + sortedFile
    print commandBase
    p=Popen(commandBase, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
    p.wait()
    os.remove(outputfilename)
    os.renames(sortedFile+'.bam', outputfilename)

def samtoolRetrieve(bamfilename, coord):
    commandBase = 'samtools view ' + bamfilename + ' ' + coord
    p=Popen(commandBase, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
    output = p.stdout.read()
    return output

def samtoolIndex(bamfilename):
    commandBase = 'samtools index ' + bamfilename
    print commandBase
    p=Popen(commandBase, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
    p.wait()

def samtoolMappingInfo(bamfilename):
    commandBase = 'samtools idxstats ' + bamfilename #the bamfile should be indexed
    p=Popen(commandBase, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
    output = p.stdout.read()
    return output

def getBamTotalRead(bamfilename):
    output = samtoolMappingInfo(bamfilename)
    output = output.split('\n')
    #print output
    #print output[-2].split('\t')
    
    read = float(output[-2].split('\t')[-1])
    return read

def getTotalHitFromBam(bamfilename):
    output = samtoolMappingInfo(bamfilename)
    output = output.split('\n')
    return sum(map(float, map(lambda x: x.split('\t')[2], output[:-1])))
  
def samtoolPrepare0(samfilename):
    #samToBam(samfilename)
    bamfilename = '.'.join(samfilename.split('.')[:-1]) + '.bam'
    delSeqInSam(samfilename, bamfilename)
    samtoolSort(bamfilename)
    samtoolIndex(bamfilename)

def samtoolPrepare(bamfilename):
    samtoolSort(bamfilename)
    samtoolIndex(bamfilename)

def samtoolPrepare2(samfilename):
    samToBam(samfilename)
    bamfilename = '.'.join(samfilename.split('.')[:-1]) + '.bam'
    samtoolSort(bamfilename)
    samtoolIndex(bamfilename)

def delSeqInSam(samfilename, bamfilename):
    def change(line):
        if line[0]!='@':
          line = line.split('\t'); line[9]='*'; line[10]='*'
          return '\t'.join(line)
        else: return line
    fin = open(samfilename); lines = fin.readlines(); fin.close()
    os.remove(samfilename)
    #lines = filter(lambda x: x[0]!='@', lines)
    lines = map(lambda x: change(x), lines)
    tempsam =samfilename +'.temp'
    fout = open(tempsam, 'w')
    fout.write(''.join(lines)); fout.close()
    
    commandBase = 'samtools view -Sb ' + tempsam + ' >' + bamfilename
    p=Popen(commandBase, shell=True)
    p.wait()
    os.remove(tempsam)
    #samtoolPrepare(bamfilename)

def delSeqInBam(bamfilename):
    tempsam = '/'.join(bamfilename.split('/')[:-1]) + '/temp.sam'
    p=Popen('samtools view -h ' + bamfilename + ' >' + tempsam, shell=True)
    p.wait()
    os.remove(bamfilename)
    def change(line): 
	if line[0]!='@':
	  line = line.split('\t'); line[9]='*'; line[10]='*'
	  return '\t'.join(line)
	else: return line
    fin = open(tempsam); lines = fin.readlines(); fin.close()
    #lines = filter(lambda x: x[0]!='@', lines)
    lines = map(lambda x: change(x), lines)
    os.remove(tempsam)
    fout = open(tempsam, 'w')
    fout.write(''.join(lines)); fout.close()

    commandBase = 'samtools view -Sb ' + tempsam + ' >' + bamfilename
    p=Popen(commandBase, shell=True)
    p.wait()
    os.remove(tempsam)
    samtoolPrepare(bamfilename)

def bedToSam(bedfile, genome):
    bamfilename = '.'.join(bedfile.split('.')[:-1]) +'.bam'
    commandBase=''
    if genome=='ce6': commandBase = 'bedtools bedtobam -i ' + bedfile +' -g /lab/solexa_bartel1/jwnam/RNASeq/Celegans/Fire_lab/ce6.genome >' + bamfilename
    elif genome=='hg19': commandBase = 'bedtools bedtobam -i ' + bedfile +' -g /lab/bartel3_ata/jwnam/h.sapiens/hg19/hg19.genome >' + bamfilename
    elif genome=='dm3': commandBase = 'bedtools bedtobam -i ' + bedfile +' -g /lab/bartel3_ata/jwnam/d.melanogaster/dm3.genome >' + bamfilename
    elif genome=='mm9': commandBase = 'bedtools bedtobam -i ' + bedfile +' -g /lab/bartel6_ata/nam/mouse/genome/mm9.genome >' + bamfilename
    elif genome=='danRer7': commandBase = 'bedtools bedtobam -i ' + bedfile +' -g /lab/bartel3_ata/jwnam/zebrafish/danRer7/danRer7.genome >' + bamfilename
    p=Popen(commandBase, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
    p.wait()
    print commandBase
    samtoolSort(bamfilename)
    samtoolIndex(bamfilename)   

def gffToSam(gffFile):
    bedFilename = '.'.join(gffFile.split('.')[:-1]) + '.bed'
    commandBase = 'python gtfToBed.py ' + gffFile
    print commandBase
    p=Popen(commandBase, shell=True)
    p.wait()
    bedToSam(bedFilename)

def getBambyChrom(bamfile, chrom, strandType='np'):
    senseCode={'+':'0', '-':'16'}
    commandBase = 'samtools view ' + bamfile + ' ' + chrom
    print commandBase
    p=Popen(commandBase, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
    output = p.stdout.read()
    output = output.split('\n')
    output = filter(lambda x: x.split('\t')[5].find('N')<0, output[:-1])
    return output

def getBam(bamfile, locus, strandType='np'):
    senseCode={'+':'0', '-':'16'}
    commandBase = 'samtools view ' + bamfile + ' ' + locus.coord()
    #print commandBase
    output=''
    #if locus.coord()!='chrY:10035874-10036679' and locus.coord()!='chr17:16342640-16343567' and locus.coord()!='chr1:145103907-145104017' and locus.coord()!='chr1:145109523-145109684': #and locus.coord()!='chr1:173835665-173835934' and locus.coord()!='chr9:139006426-139008870':
    if True:
    	#p=Popen(commandBase, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
    	p=Popen(commandBase, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
    	#outputE = p.stderr.read(); print p
    	output = p.stdout.read()
    	output = output.split('\n')
	#print len(output)
    	sense = locus.sense()
    	if strandType=='np': output = filter(lambda x: x.split('\t')[5].find('N')<0, output[:-1])
    	else: output = filter(lambda x: x.split('\t')[1]==senseCode[sense], output[:-1])
	#print commandBase
	#print output
    return output

def bamToBedGraph(bamfilename, chromInfo, strand='.'):
    def update(x):
		x = x.split('\t')
		x[-1] = '-'+x[-1]
		return '\t'.join(x)
    if strand!='.':
        bedGraphfilename = '.'.join(bamfilename.split('.')[:-1]) + '.' + strand + '.bed'
	print 'genomeCoverageBed -bg -strand ' + strand + ' -ibam ' + bamfilename + ' -g ' + chromInfo + ' > ' + bedGraphfilename
    	a = Popen('genomeCoverageBed -bg -strand ' + strand + ' -ibam ' + bamfilename +' -g ' + chromInfo + ' > ' + bedGraphfilename, shell=True) # compute coverage
    	a.wait()
	if strand=='-':
		filein = open(bedGraphfilename)
		lines = filein.readlines();filein.close()
		lines = map(lambda x: update(x), lines)
        	bedGraphfilenameT = bedGraphfilename + '.tmp'
		fileout = open(bedGraphfilenameT, 'w')
		fileout.write(''.join(lines)); fileout.close()
		b = Popen('rm ' + bedGraphfilename, shell=True); b.wait()
		b = Popen('mv ' + bedGraphfilenameT + ' ' + bedGraphfilename, shell=True); b.wait()
    else:
        bedGraphfilename = '.'.join(bamfilename.split('.')[:-1]) + '.bed'
	print 'genomeCoverageBed -bg -ibam ' + bamfilename + ' -g ' + chromInfo + ' > ' + bedGraphfilename
    	a = Popen('genomeCoverageBed -bg -ibam ' + bamfilename +' -g ' + chromInfo + ' > ' + bedGraphfilename, shell=True) # compute coverage
    	a.wait()

def bamToCoverage(bamfilename): 
	command = 'samtools depth ' + bamfilename
	p = Popen(command, shell=True, stdout=PIPE, stdin=PIPE, stderr=PIPE)
 	return p.stdout
	
#import sys
#gffToSam(sys.argv[1]) 
#bedfile = sys.argv[1]    
#bedToSam(bedfile, sys.argv[2])
#bamfilename = '.'.join(bedfile.split('.')[:-1]) +'.bam'
#bamToBedGraph(bamfilename, '/lab/solexa_bartel1/jwnam/RNASeq/Celegans/Fire_lab/ce6.genome' , '+')
#bamToBedGraph(bamfilename, '/lab/solexa_bartel1/jwnam/RNASeq/Celegans/Fire_lab/ce6.genome' , '-')

'''
bedfiles=[0]*9
bedfiles[8] = '/home/jwnam/theOther/Project/Celegans/3UTR/rawData/mixed_b_pA.bed'
bedfiles[0] = '/home/jwnam/theOther/Project/Celegans/3UTR/rawData/mixed_a_pA.bed'
bedfiles[1] = '/home/jwnam/theOther/Project/Celegans/3UTR/rawData/L1_pA.bed'
bedfiles[2] = '/home/jwnam/theOther/Project/Celegans/3UTR/rawData/L2_pA.bed'
bedfiles[3] = '/home/jwnam/theOther/Project/Celegans/3UTR/rawData/L3_pA.bed'
bedfiles[4] = '/home/jwnam/theOther/Project/Celegans/3UTR/rawData/L4_pA.bed'
bedfiles[5] = '/home/jwnam/theOther/Project/Celegans/3UTR/rawData/adult_pA.bed'
bedfiles[6] = '/home/jwnam/theOther/Project/Celegans/3UTR/rawData/glp4_pA.bed'
bedfiles[7] = '/home/jwnam/theOther/Project/Celegans/3UTR/rawData/dauer_pA.bed'
#for bedfile in bedfiles: bedToSam(bedfile)
bamfiles=[]
for bedfile in bedfiles: bamfile = '.'.join(bedfile.split('.')[:-1]) + '.bam'; bamfiles.append(bamfile)
samtoolMerge(bamfiles, 'all_3PSeq.bam')
'''
