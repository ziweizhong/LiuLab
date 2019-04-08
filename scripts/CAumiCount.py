import numpy as np
import ast
import pickle

def getFreq(caDict):
	c = []
	f = []

	for i in caDict:
		if i != '' or i != 'count':
			for j in caDict[i]:
				if j != 'count':
					for k in caDict[i][j]:
						if k != 'count':
							try:
								# print("i:%s j:%s k:%s"%(i,j,k))
								# print(caDict[i][j][k]['count'])
								ind = c.index(caDict[i][j][k]['count'])
								f[ind] += 1
							except:
								# print(caDict[i][j][k]['count'])
								c.append(caDict[i][j][k]['count'])
								f.append(1)

	return c, f


count = [1]
freq = [0]

for j in range(1,9):
	for k in range(1,4):
		run = 1
		if j == 1 and k == 1:
			run = 0
		elif j == 6 and k == 1:
			with open('out.CA_9_1.fastq.pickle', 'rb') as handle:
				data = pickle.load(handle)
		else:
			with open('out.CA_%d_%d.fastq.pickle'%(j,k), 'rb') as handle:
				data = pickle.load(handle)



		if run:
			a,b = getFreq(data)
			for i in range(0,len(a)):
				try:
					ind = count.index(a[i])
					freq[ind] += b[i]
				except:
					count.append(a[i])
					freq.append(b[i])

print(count)
print(freq)


