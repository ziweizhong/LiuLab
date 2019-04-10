import numpy as np
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



def getMuts(muts):
	# if s == a:
	# 	mutPos = mutPosV3
	# elif s == b:
	# 	mutPos = mutPosV1
	# elif s == c:
	# 	mutPos = mutPosV2
	# elif s == d:
	# 	mutPos = mutPosV4
	# elif s == e or s == f:
	# 	mutPos = mutPosV5
	# else:
	# 	mutPos = []

	mutPos = mutPosV5

	mutPresence = np.zeros(len(mutPos))
	if muts == 'WT':
		return mutPresence
	else:
		sepPos = [i for i, a in enumerate(str(muts)) if a == '-']
		# print (mutPos)
		# print(int(muts))
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




asdf = getMuts(sys.argv[1])
print(mutPosV5)
print(asdf)








