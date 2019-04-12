import numpy as np
import ast
import pickle
from scipy import stats
import sys

a = [2,7,13,19]
b = [3,8,14,20]
c = [4,9,15,21]
d = [5,10,16,22]
e = [6,11,17,23]
f = [6,12,18,24]
aC = [17521,18403,22831,18130]
bC = [20911,14551,24096,16992]
cC = [16315,18292,24545,15282]
dC = [12257,17125,20647,14510]
eC = [22194,22814,15626,17654]
fC = [22194,21371,13984,17074]

mutPosV1 = [37, 71, 85, 124, 182, 209, 226]
mutPosV2 = [15, 55, 71, 72, 74, 92, 139, 190, 209, 235]
mutPosV3 = [7, 15, 36, 61, 124, 143, 207, 234]
mutPosV4 = [15, 37, 167, 185, 209]
mutPosV5 = [60, 71, 81, 107, 125, 127, 133, 235]

eSets = [a,b,c,d,e,f]


#########################################
#	Groups results based on which set 	#
#	they belong to. For example, set a, #
#	which includes wells 2, 7, 13 and 19#
#	are all timepoints of the same 		#
#	library. Thus, data from them are	#
#	grouped.							#
#########################################
def collateResults(indMatrix, sel = ''):
	dataA = []

	for i in indMatrix:
		index = [((i-1)%8)+1, int((i-1)/8+1)]
		# print(index)
		if index[0] == 6 and index[1] == 1:
			with open('out.CA_9_1.fastq.%spickle'%(sel), 'rb') as handle:
				print('out.CA_9_1.fastq.%spickle'%(sel))
				data = pickle.load(handle)
			dataA.append(data)
		else:
			with open('out.CA_%d_%d.fastq.%spickle'%(index[0],index[1],sel), 'rb') as handle:
				print('out.CA_%d_%d.fastq.%spickle'%(index[0],index[1],sel))
				data = pickle.load(handle)
			dataA.append(data)

	results = {}

	for i in range(0,len(dataA)):
		for muts in dataA[i]:
			if muts != '':
				if muts in results:
					results[muts][i] = dataA[i][muts]['count']
				else:
					results[muts] = {i:dataA[i][muts]['count']}
	return results






#########################################################################
#	Changes a string of mutations into an array with true				#
#	or false based on the presence or absence of the mutation.			#
#		For example, for mutPosV1 = [37, 71, 85, 124, 182, 209, 226]	#
#		If the string is '37-209', the following array will be returned	#
#		[1,0,0,0,0,1,0]
#########################################################################
def getMuts(muts,mutPos):
	mutPresence = np.zeros(len(mutPos))
	if muts == 'WT':
		return mutPresence
	else:
		sepPos = [i for i, a in enumerate(str(muts)) if a == '-']
		if len(sepPos) == 0:
			try:
				mutInd = mutPos.index(int(muts))
				mutPresence[mutPos.index(int(muts))] = 1
			except:
				pass
		else:
			mutsMtx = []
			mutsMtx.append(int(muts[0:sepPos[0]]))
			for i in range(0,len(sepPos)-1):
				mutsMtx.append(int(muts[(sepPos[i]+1):sepPos[i+1]]))
			mutsMtx.append(int(muts[(sepPos[len(sepPos)-1]+1):len(muts)]))
			
			for i in mutsMtx:
				try:
					mutInd = mutPos.index(i)
					mutPresence[mutInd] = 1
				except:
					continue
					

	return mutPresence

# print(results)







#########################################################################
#	Start of actual script												#
#########################################################################
try:
	pickleSel = sys.argv[1].upper() + '.'
except:	
	pickleSel = ''

print(pickleSel)



for s in eSets:
	fitnessDict = {}
	counts = [1,1,1,1]
	gen = [0.0,7.64,6.64,5.64]

	if s == a:
		counts = aC
		mutPos = mutPosV3
		name = '2-7'
	elif s == b:
		counts = bC
		mutPos = mutPosV1
		name = '3-8'
	elif s == c:
		counts = cC
		mutPos = mutPosV2
		name = '4-9'
	elif s == d:
		counts = dC
		mutPos = mutPosV4
		name = '5-10'
	elif s == e:
		counts = eC
		gen = [0.0,6.64,5.64,5.64]
		mutPos = mutPosV5
		name = '6-11'
	elif s == f:
		counts = fC
		gen = [0.0,6.64,5.64,5.64]
		mutPos = mutPosV5
		name = '6-12'

	with open ('0Output-%s-enrich.txt'%(str(s)),'w') as outFit:
		outFit.write("Name,Fitness,Count t0,Count t1,Count t2,Count t3,Total Count,R-squared,Number of Mutations")
		for i in range(0,len(mutPos)):
			outFit.write(","+str(mutPos[i]))
		outFit.write("\n")

		results = collateResults(s,pickleSel)
		print(s)
		
		for muts in results:
			mutMtx = getMuts(muts,mutPos)
			
			freq = [0.0,0.0,0.0,0.0]
			counts2 = [0,0,0,0]
			totalCount = 0

			for i in results[muts]:
				freq[i] = results[muts][i]/counts[i]
				totalCount += results[muts][i]
				counts2[i] = results[muts][i]


			fit = np.polyfit(gen,freq,deg=1)
			slope, intercept, r_value, p_value, std_err = stats.linregress(gen, freq)

			numMuts = len([i for i, a in enumerate(str(muts)) if a == '-']) + 1

			outFit.write("%s,%f,%d,%d,%d,%d,%d,%f,%d"%(str(muts),slope*100,counts2[0], counts2[3], counts2[2], counts2[1],totalCount,r_value**2,numMuts))
			for i in range(0,len(mutMtx)):
				outFit.write(","+str(mutMtx[i]))
			outFit.write("\n")

			fitnessDict[muts] = slope*100

	with open('fitness-'+name+'.pickle', 'wb') as handle:
		pickle.dump(fitnessDict, handle)

