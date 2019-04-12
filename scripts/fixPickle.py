import pickle
import sys
import pprint

selectedMuts = [60, 71, 81, 107, 125, 127, 133, 235]

with open(sys.argv[1], 'rb') as handle:
	print('%s loaded!'%sys.argv[1])
	dataOrg = pickle.load(handle)

newDict = {}

for muts in dataOrg:
# muts = sys.argv[1]
	sepPos = [i for i, a in enumerate(str(muts)) if a == '-']
	mutsMtx = []
	if muts == '':
		newMuts = 'WT'
	elif muts == 'WT':
		newMuts = muts
	elif len(sepPos) == 0:
		mutsMtx.append(int(muts))
	else:
		mutsMtx.append(int(muts[0:sepPos[0]]))
		for i in range(0,len(sepPos)-1):
			mutsMtx.append(int(muts[(sepPos[i]+1):sepPos[i+1]]))
		mutsMtx.append(int(muts[(sepPos[len(sepPos)-1]+1):len(muts)]))

	if len(mutsMtx) > 0:
		newMuts = ''
		for i in mutsMtx:
			if i in selectedMuts:
				newMuts += (str(i) + '-')

		newMuts = newMuts[0:len(newMuts)-1]
	# print(muts)
	# print(mutsMtx)
	# print(newMuts)
	if newMuts in newDict:
		pp = pprint.PrettyPrinter(indent=4)
		# print("Before Merge")
		# pp.pprint(newDict[newMuts])

		# print("To Merge")
		# pp.pprint(dataOrg[muts])
		newCount = newDict[newMuts]['count'] + dataOrg[muts]['count']
		newDict[newMuts] = {**newDict[newMuts],**dataOrg[muts]}
		newDict[newMuts]['count'] = newCount
		# print("After Merge")
		# pp.pprint(newDict[newMuts])

	else:
		newDict[newMuts] = dataOrg[muts]



with open(sys.argv[2], 'wb') as handle:
		pickle.dump(newDict, handle)
