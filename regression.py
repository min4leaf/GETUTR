import sys

def midCut(list):
        ret = []
        for i in list:
                ret.append(i)
        cut_mode = 0
        cut_value = list[0]
        max = ret[0]
        min = ret[0]
        inc = 1
        for i in xrange(1,len(ret)):
                if cut_mode==1:
                        if ret[i]<cut_value:
                                ret[i] = cut_value
                        else:
                                inc = 1
                                max = ret[i]
                                cut_mode = 0
                else:
                        if ret[i-1]>ret[i]:
                                if i==len(ret)-1:
                                        cut_value = (max+ret[i])/2.0
                                        ret[i] = cut_value
                                        ptr = i-1
                                        while (ptr>-1) and (ret[ptr]<cut_value):
                                                ret[ptr] = cut_value
                                                ptr -=1
                                        while (ptr>-1) and (ret[ptr]>cut_value):
                                                ret[ptr] = cut_value
                                                ptr -=1
                                else:
                                        inc = 0 #increase
                                        min = ret[i]
                        else:
                                if inc==0:
                                        if max<ret[i]:
                                                max = ret[i]
                                        cut_value = (min+max)/2.0
                                        if ret[i]<cut_value:
                                                ret[i] = cut_value
                                                cut_mode = 1
                                        ptr = i-1
                                        while (ptr>-1) and (ret[ptr]<cut_value):
                                                ret[ptr] = cut_value
                                                ptr -= 1
                                        while (ptr>-1) and (ret[ptr]>cut_value):
                                                ret[ptr] = cut_value
                                                ptr -= 1
                                inc = 1
                                max = ret[i]
        return ret

def isoReg(list, w):
	y = dict()
	weight = dict()
	s = dict()

	y[0] = list[0]
	weight[0] = w[0]
	ptr = 0
	s[-1] = 0
	s[0] = 1
	for i in xrange(1,len(list)):
		ptr +=1
		y[ptr] = list[i]
		weight[ptr] = w[i]
		while (ptr>0) and (y[ptr]<y[ptr-1]):
			y[ptr-1] = (weight[ptr]*y[ptr] + weight[ptr-1]*y[ptr-1])/(weight[ptr]+weight[ptr-1])
                        weight[ptr-1] = weight[ptr] + weight[ptr-1]
                        ptr -= 1
		s[ptr] = i
	retDic = dict()
	for i in xrange(0,ptr+1):
		for j in xrange(s[i-1],s[i]+1):
			retDic[j] = y[i]
	ret = []
	for key in retDic.keys():
		ret.append(retDic[key])	
	return ret
	'''
	ptr = 0
        y = dict()
        weight = dict()
        y[0] = list[0]
        weight[0] = w[0]
        for i in xrange(1,len(list)):
                ptr += 1
                y[ptr] = list[i]
                weight[ptr] = w[i]
                while (ptr>0) and (y[ptr-1] > y[ptr]):
                        y[ptr-1] = (weight[ptr]*y[ptr] + weight[ptr-1]*y[ptr-1])/(weight[ptr]+weight[ptr-1])
                        weight[ptr-1] = weight[ptr] + weight[ptr-1]
                        ptr -= 1
        ret_spread = []
	ret = []
        for i in xrange(len(y)):
                for j in xrange(0,weight[i]):
                        ret_spread.append(y[i])
                if len(ret_spread)==len(list):
                        break;
	ptr = -1
	if len(ret_spread)<w[len(w)-1]:
		print "!!!!!!"
		print len(w)
		print len(list)
		print len(ret_spread)
		print w
		print list
		print ret_spread
	for i in xrange(len(w)):
		ptr += w[i]
		ret.append(ret_spread[ptr])
        return ret_spread
	'''

