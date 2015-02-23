def cleavDetect(UTR):
	# input::
	# dict UTR = [string geneid:[string chr, int start_pos, int end_pos, float value, string sense]]
	# output::
	# dict cleavage = [string geneid:[string chr, string sense, int pos, float value]]
	cleavage = dict()
	for geneid in UTR.keys():
		grad = []
		if len(UTR[geneid])>0:
			cleavage[geneid] = []
			chr = UTR[geneid][0][0]
			sense = UTR[geneid][0][4]
			if sense=='+':
				for i in xrange(len(UTR[geneid])-1):
					diff = UTR[geneid][i][3] - UTR[geneid][i+1][3]
					if diff!=0:
						grad.append([i, diff])
				grad.append([len(UTR[geneid])-1, UTR[geneid][len(UTR[geneid])-1][3]])
			else:
				grad.append([0, -UTR[geneid][0][3]])
				for i in xrange(len(UTR[geneid])-2):
					diff = UTR[geneid][i+1][3] - UTR[geneid][i+2][3]
					if diff!=0:
						grad.append([i+1, diff])
			if len(grad)==1:
				ptr = grad[0][0]
				if sense=='+':
					cleavage[geneid].append([chr, sense, UTR[geneid][ptr][2], grad[0][1]])
				else:
					cleavage[geneid].append([chr, sense, UTR[geneid][ptr][1], -grad[0][1]])
			elif len(grad)==2:
				if sense=='+':
					if grad[0][1]>grad[1][1]:
						ptr = grad[0][0]
						cleavage[geneid].append([chr, sense, UTR[geneid][ptr][2], grad[0][1]])
					else:
						ptr = grad[1][0]
						cleavage[geneid].append([chr, sense, UTR[geneid][ptr][2], grad[1][1]])
				else:
					if grad[0][1]<grad[1][1]:
						ptr = grad[0][0]
						cleavage[geneid].append([chr, sense, UTR[geneid][ptr][1], -grad[0][1]])
					else:
						ptr = grad[1][0]
						cleavage[geneid].append([chr, sense, UTR[geneid][ptr][1], -grad[1][1]])
			elif len(grad)>2:
				if sense=='+':
					if grad[0][1]>grad[1][1]:
						ptr = grad[0][0]
						cleavage[geneid].append([chr, sense, UTR[geneid][ptr][2], grad[0][1]])
					for i in xrange(1, len(grad)-1):
						if grad[i][1]>grad[i-1][1] and grad[i][1]>grad[i+1][1]:
							ptr = grad[i][0]
							cleavage[geneid].append([chr, sense, UTR[geneid][ptr][2], grad[i][1]])
					if grad[len(grad)-1][1]>grad[len(grad)-2][1]:
						ptr = grad[len(grad)-1][0]
						cleavage[geneid].append([chr, sense, UTR[geneid][ptr][2], grad[len(grad)-1][1]])
				else:
					if grad[0][1]<grad[1][1]:
						ptr = grad[0][0]
						cleavage[geneid].append([chr, sense, UTR[geneid][ptr][1], -grad[0][1]])
					for i in xrange(1, len(grad)-1):
						if grad[i][1]<grad[i-1][1] and grad[i][1]<grad[i+1][1]:
							ptr = grad[i][0]
							cleavage[geneid].append([chr, sense, UTR[geneid][ptr][1], -grad[i][1]])
					if grad[len(grad)-1][1]<grad[len(grad)-2][1]:
						ptr = grad[len(grad)-1][0]
						cleavage[geneid].append([chr, sense, UTR[geneid][ptr][1], -grad[len(grad)-1][1]])
	return cleavage

def makeSeq(UTR, cleavage):
	# input::
	# dict UTR = [string geneid:[string chr, int start_pos, int end_pos, float value, string sense]]
	# dict cleavage = [string geneid:[string chr, string sense, int pos, float value]]
	# output::
	# dict nUTR = [string geneid:[string chr, int start_pos, int end_pos, float value, string sense]]
	# dict cleavage = [string geneid:[string chr, string sense, int pos, float value]]
	nUTR = dict()
	for geneid in UTR.keys():
		if cleavage.has_key(geneid):
			nUTR[geneid] = []
			if len(cleavage[geneid])>0 and len(UTR[geneid])>0:
				chr = cleavage[geneid][0][0]
				sense = cleavage[geneid][0][1]
				ptr = 0
				maxLength = -1
				maxPtr = -1
				for i in xrange(len(cleavage[geneid])):
					maxLength = -1
					if sense=='+':
						ptr = 0
						while UTR[geneid][ptr][2]<cleavage[geneid][i][2]:
							if abs(UTR[geneid][ptr][2] - UTR[geneid][ptr][1])>maxLength:
								maxLength = abs(UTR[geneid][ptr][2] - UTR[geneid][ptr][1])
								maxPtr = ptr
							ptr += 1
						cleavage[geneid][i][3] = UTR[geneid][maxPtr][3]
					else:
						ptr = len(UTR[geneid])-1
						while cleavage[geneid][-i-1][2]<UTR[geneid][ptr][1]:
							if abs(UTR[geneid][ptr][2] - UTR[geneid][ptr][1])>maxLength:
								maxLength = abs(UTR[geneid][ptr][2] - UTR[geneid][ptr][1])
								maxPtr = ptr
							ptr -= 1
						cleavage[geneid][-i-1][3] = UTR[geneid][maxPtr][3]
				if sense=='+':
					if cleavage[geneid][0][3]!=1:
						for i in xrange(1,len(cleavage[geneid])+1):
							cleavage[geneid][-i][3] /= cleavage[geneid][0][3]
					for i in xrange(len(cleavage[geneid])-1):
						cleavage[geneid][i][3] -= cleavage[geneid][i+1][3]
				else:
					if cleavage[geneid][-1][3]!=1:
						for i in xrange(0,len(cleavage[geneid])):
							cleavage[geneid][i][3] /= cleavage[geneid][-1][3]
					for i in xrange(1, len(cleavage[geneid])):
						cleavage[geneid][-i][3] -= cleavage[geneid][-i-1][3]
				
				if sense=='+':
					start = UTR[geneid][0][1]
					end = cleavage[geneid][0][2]
					value = 1.0
					for i in xrange(1, len(cleavage[geneid])):
						nUTR[geneid].append([chr, start, end, value, sense])
						start = cleavage[geneid][i-1][2]
						end = cleavage[geneid][i][2]
						value -= cleavage[geneid][i-1][3]
					nUTR[geneid].append([chr, start, end, value, sense])
					#nUTR.append([chr, sense, end, list[-1][3], 0.0])
				else:
					value = 0.0
					for i in xrange(0, len(cleavage[geneid])-1):
						start = cleavage[geneid][i][2]
						end = cleavage[geneid][i+1][2]
						value += cleavage[geneid][i][3]
						nUTR[geneid].append([chr, start, end, value, sense])
					nUTR[geneid].append([chr, sense, cleavage[geneid][-1][2], UTR[geneid][-1][2], 1.0])
	return nUTR, cleavage
