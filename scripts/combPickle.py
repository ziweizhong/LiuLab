import pickle
import sys
import pprint


with open(sys.argv[1], 'rb') as handle:
	print('%s loaded!'%sys.argv[1])
	dataOne = pickle.load(handle)

with open(sys.argv[2], 'rb') as handle:
	print('%s loaded!'%sys.argv[2])
	dataTwo = pickle.load(handle)


for muts in dataTwo:
	
	if muts in dataOne:
		
		newCount = dataOne[muts]['count'] + dataTwo[muts]['count']
		dataOne[muts] = {**dataOne[muts],**dataTwo[muts]}
		dataOne[muts]['count'] = newCount


	else:
		dataOne[muts] = dataTwo[muts]



with open(sys.argv[3], 'wb') as handle:
		pickle.dump(dataOne, handle)
