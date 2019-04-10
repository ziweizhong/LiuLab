import pickle


for i in range(0,10):
	for j in range(0,4):
		print('out.CA_%d_%d.fastq.pickle'%(i,j))

		try:
			with open('out.CA_%d_%d.fastq.pickle'%(i,j), 'rb') as handle:
				dataWithWT = pickle.load(handle)
			wtData = dataWithWT["WT"]
		except:
			break

		with open('out.CA_%d_%d.fastq.LENIENT.pickle'%(i,j), 'rb') as handle:
			tempData = pickle.load(handle)
		tempData["WT"] = wtData

		with open('out.CA_%d_%d.fastq.LENIENT.pickle'%(i,j), 'wb') as handle:
			pickle.dump(tempData, handle)

		with open('out.CA_%d_%d.fastq.ALL.pickle'%(i,j), 'rb') as handle:
			tempData = pickle.load(handle)
		tempData["WT"] = wtData

		with open('out.CA_%d_%d.fastq.ALL.pickle'%(i,j), 'wb') as handle:
			pickle.dump(tempData, handle)