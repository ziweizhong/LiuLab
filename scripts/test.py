import pickle
import sys

selectedMuts = [2,4,6,8,10]

# try:
totalCounts = 0
blankCount = 0
with open(sys.argv[1], 'rb') as handle:
	print('%s loaded!'%sys.argv[1])
	dataOrg = pickle.load(handle)

for muts in dataOrg:
# muts = sys.argv[1]
	if muts == '':
		blankCount += dataOrg[muts]['count']
	totalCounts += dataOrg[muts]['count']

try:
	print("WT Counts: %d"%dataOrg['WT']['count'])
except:
	print("WT Counts: NONE!")
print("Blank Counts: %d"%blankCount)
print("Total Counts: %d"%totalCounts)
# except:
# 	for i in range(0,10):
# 		for j in range(0,4):
# 			fn = "out.CA_%d_%d.fastq.pickle"%(i,j)
# 			totalCounts = 0
# 			blankCount = 0
# 			with open(fn, 'rb') as handle:
# 				print('%s loaded!'%fn)
# 				dataOrg = pickle.load(handle)

# 			for muts in dataOrg:
# 			# muts = sys.argv[1]
# 				if muts == '':
# 					blankCount += dataOrg[muts]['count']
# 				totalCounts += dataOrg[muts]['count']

# 			print("Blank Counts: %d"%blankCount)
# 			print("Total Counts: %d"%totalCounts)