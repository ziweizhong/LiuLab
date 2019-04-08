import numpy as np
import ast
import pickle

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

eSets = [a,b,c,d,e,f]

def collateResults(indMatrix):
	dataA = []
	for i in indMatrix:
		index = [((i-1)%8)+1, int((i-1)/8+1)]
		# print (index)
		if index[0] == 6 and index[1] == 1:
			with open('out.CA_9_1.fastq.pickle', 'rb') as handle:
	 			data = pickle.load(handle)
			dataA.append(data)
		else:
			with open('out.CA_%d_%d.fastq.pickle'%(index[0],index[1]), 'rb') as handle:
				data = pickle.load(handle)
			dataA.append(data)

	# print(type(dataA))
	# print(len(dataA))
	# print(type(dataA[1]))



	#########################
	#        SET A          #
	#########################


	results = {}

	for i in range(0,len(dataA)):
		# print(dataA[i])
		for muts in dataA[i]:
			if muts != '':
				if muts in results:
					results[muts][i] = dataA[i][muts]['count']
				else:
					results[muts] = {i:dataA[i][muts]['count']}
				# print(type(results))
				# print(type(dataA))
				# print(muts)
				# print(dataA[i][muts]['count'])
				# results[muts] = {'count':dataA[i][muts]['count']}
	return results

# print(results)
with open ('enrichData.txt','w') as out:
	with open ('enrichDataValues.txt','w') as outVal:
		for s in eSets:
			results = collateResults(s)
			print(s)
			out.write(str(s))
			out.write(',--,--,---,---,')
			out.write('\n')
			counts = [1,1,1,1]

			gen = [0.0,7.64,6.64,5.64]
			if s == a:
				counts = aC
			elif s == b:
				counts = bC
			elif s == c:
				counts = cC
			elif s == d:
				counts = dC
			elif s == e:
				counts = eC
				gen = [0.0,6.64,5.64,5.64]
			elif s == f:
				counts = fC
				gen = [0.0,6.64,5.64,5.64]

			for muts in results:
				freq = [0.0,0.0,0.0,0.0]
				for i in results[muts]:		
					out.write("%s,%d,%d,%f\n"%(muts,i,results[muts][i],results[muts][i]/counts[i]))
					freq[i] = results[muts][i]/counts[i]

				fit = np.polyfit(gen,freq,deg=1)
				out.write("%s fitness:\t%f\n"%(str(muts),fit[0]))

				outVal.write(str(muts))
				outVal.write('\n')
				outVal.write(str(gen)[1:len(str(gen))-1])
				outVal.write('\n')
				outVal.write(str(freq)[1:len(str(freq))-1])
				outVal.write('\n')
				outVal.write('int:\t%f\nslope:\t%f\n'%(fit[1],fit[0]))


for s in eSets:
	with open ('%s-enrich.txt'%(str(s)),'w') as out:
		with open ('%s-enrich-RAW.txt'%(str(s)),'w') as outVal:
			with open ('%s-enrich-FITTED.txt'%(str(s)),'w') as outFit:
				results = collateResults(s)
				print(s)
				out.write(str(s))
				out.write(',--,--,---,---,')
				out.write('\n')
				counts = [1,1,1,1]

				gen = [0.0,7.64,6.64,5.64]
				if s == a:
					counts = aC
				elif s == b:
					counts = bC
				elif s == c:
					counts = cC
				elif s == d:
					counts = dC
				elif s == e:
					counts = eC
					gen = [0.0,6.64,5.64,5.64]
				elif s == f:
					counts = fC
					gen = [0.0,6.64,5.64,5.64]

				for muts in results:
					freq = [0.0,0.0,0.0,0.0]
					totalCount = 0
					for i in results[muts]:		
						out.write("%s,%d,%d,%f\n"%(muts,i,results[muts][i],results[muts][i]/counts[i]))
						freq[i] = results[muts][i]/counts[i]
						totalCount += results[muts][i]

					fit = np.polyfit(gen,freq,deg=1)
					out.write("%s fitness:\t%f\n"%(str(muts),fit[0]))
					outFit.write("%s fitness:\t%f\t%d\n"%(str(muts),fit[0],totalCount))
					outVal.write(str(muts))
					outVal.write('\n')
					outVal.write(str(gen)[1:len(str(gen))-1])
					outVal.write('\n')
					outVal.write(str(freq)[1:len(str(freq))-1])
					outVal.write('\n')
					outVal.write('int:\t%f\nslope:\t%f\n'%(fit[1],fit[0]))




# for k in range(1,4):
# 	row = []

# 	for j in range(1,9):
# 		run = 1

# 		print("%d %d"%(k,j))

# 		if j == 1 and k == 1:
# 			run = 0
# 			row.append({})
# 		elif j == 6 and k == 1:
# 			with open('out.CA_9_1.fastq.pickle', 'rb') as handle:
# 				data = pickle.load(handle)
# 			row.append(data)
# 		else:
# 			with open('out.CA_%d_%d.fastq.pickle'%(j,k), 'rb') as handle:
# 				data = pickle.load(handle)
# 			row.append(data)





