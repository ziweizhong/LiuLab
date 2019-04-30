import pickle
import sys
import pprint


with open(sys.argv[1], 'rb') as handle:
	print('%s loaded!'%sys.argv[1])
	dataOne = pickle.load(handle)

with open(sys.argv[2], 'rb') as handle:
	print('%s loaded!'%sys.argv[2])
	dataTwo = pickle.load(handle)

r = 0
for muts in dataTwo:
	if r < 5:
		r += 1
		print (dataOne[muts]['count'])
		print (dataTwo[muts]['count'])
		if muts in dataOne:
			
			newCount = dataOne[muts]['count'] + dataTwo[muts]['count']
			dataOne[muts] = {**dataOne[muts],**dataTwo[muts]}
			dataOne[muts]['count'] = newCount


		else:
			dataOne[muts] = dataTwo[muts]
		print (dataOne[muts]['count'])


with open(sys.argv[3], 'wb') as handle:
		pickle.dump(dataOne, handle)
