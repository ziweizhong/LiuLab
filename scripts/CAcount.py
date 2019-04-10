from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import Align
import numpy as np
import multiprocessing as mp
import datetime
import sys
import pprint
import pickle

print(datetime.datetime.now())



def mtAlign(data):
	limitReads = 8000000
	aligner = Align.PairwiseAligner()

	BARCODE_SEQUENCE = Seq("aaatcataagaaattcgcNNNNNNNNNNNNNNTTAGCGAGCATATCTCTTC".upper(),generic_dna)
	UMI_SEQUENCE = Seq("GACGTGTGCTCTTCCGATCnnnnnnnnnnnnnnnnnnnncataaatcataagaaattcgc".upper(),generic_dna)

	HAref = Seq("TGTAAAATGCATGTTGTCCCGGCGATAGATCTCTTCAGAGGAAAGGTAGCGAGGATGATAAAAGGAAGAAAAGAGAACACCATATTTTACGAAAAAGATCCCGTAGAACTGGTGGAAAAACTCATCGAAGAGGGATTCACACTGATTCACGTGGTGGATCTCTCGAATGCGATAGAAAACAGCGGCGAGAATCTTCCAGTTCTCGAGAAACTCTCTGAATTTGCCGAGCACATACAGATCGGAGGCGGGATCAGATCGCTCGATTACGCGGAAAAACTCCGAAAGCTGGGATACAGAAGACAGATCGTGAGCTCAAAGGTTCTGGAAGATCCTTCTTTCCTGAAATCCCTGAGAGAAATCGATGTGGAGCCCGTGTTCAGTCTGGACACTCGAGGTGGAAGAGTAGCGTTCAAAGGGTGGCTGGCGGAAGAGGAGATCGACCCTGTTTCTCTTCTGAAGAGACTGAAAGAATACGGCCTTGAAGAGATCGTACACACGGAGATCGAAAAAGATGGCACTCTTCAGGAGCACGATTTTTCTCTCACCAAAAAGATAGCGATCGAAGCTGAAGTGAAAGTCCTCGCAGCGGGTGGTATCTCTTCGGAGAACTCTTTGAAAACAGCGCAGAAGGTTCACACAGAAACGAACGGGCTTCTCAAAGGTGTGATCGTGGGAAGGGCGTTTCTGGAGGGAATTCTCACAGTTGAGGTGATGAAGAGATATGCTCGCTAA", generic_dna)

	HArefFINAL = Seq("ATGCATGTTGTCCCGGCGATAGATCTCTTCAGAGGAAAGGTAGCGAGGATGATAAAAGGAAGAAAAGAGAACACCATATTTTACGAAAAAGATCCCGTAGAACTGGTGGAAAAACTCATCGAAGAGGGATTCACACTGATTCACGTGGTGGATCTCTCGAATGCGATAGAAAACAGCGGCGAGAATCTTCCAGTTCTCGAGAAACTCTCTGAATTTGCCGAGCACATACAGATCGGAGGCGGGATCAGATCGCTCGATTACGCGGAAAAACTCCGAAAGCTGGGATACAGAAGACAGATCGTGAGCTCAAAGGTTCTGGAAGATCCTTCTTTCCTGAAATCCCTGAGAGAAATCGATGTGGAGCCCGTGTTCAGTCTGGACACTCGAGGTGGAAGAGTAGCGTTCAAAGGGTGGCTGGCGGAAGAGGAGATCGACCCTGTTTCTCTTCTGAAGAGACTGAAAGAATACGGCCTTGAAGAGATCGTACACACGGAGATCGAAAAAGATGGCACTCTTCAGGAGCACGATTTTTCTCTCACCAAAAAGATAGCGATCGAAGCTGAAGTGAAAGTCCTCGCAGCGGGTGGTATCTCTTCGGAGAACTCTTTGAAAACAGCGCAGAAGGTTCACACAGAAACGAACGGGCTTCTCAAAGGTGTGATCGTGGGAAGGGCGTTTCTGGAGGGAATTCTCACAGTTGAGGTGATGAAGAGATATGCTCGCTAA", generic_dna)

	#with open ('out.'+data,'a') as outFile:
	#	outFile.write("%s\n"%(datetime.datetime.now()))

	caF = Seq("tcaggctgctttccaaatgc".upper(), generic_dna)
	caR = Seq("GGACGTGTGCTCTTCCGATC".upper(), generic_dna)

	file = open(data,'r',buffering=1000)
	i = 0
	r = -1

	sampleIndex = [(i+1) for i,a in enumerate(data) if a == '_']
	if sampleIndex[0] == 9 and sampleIndex[1] == 1:
		sampleNum = 6
	else:
		sampleNum = (int(data[sampleIndex[1]]) - 1) * 8 + int(data[sampleIndex[0]])

	# print(data)
	# print(sampleIndex[0])
	# print(sampleIndex[1])
	# print(sampleNum)

	a = [2,7,13,19]
	b = [3,8,14,20]
	c = [4,9,15,21]
	d = [5,10,16,22]
	e = [6,11,12,17,18,23,24]

	resultsDict = {}
	resultsDictALL = {}
	resultsDictLENIENT = {}
	ntDict = {}
	aaDict = {}

	mutPosV1 = [37, 71, 85, 124, 182, 209, 226]
	mutPosV2 = [15, 55, 71, 72, 74, 92, 139, 190, 209, 235]
	mutPosV3 = [7, 15, 36, 61, 124, 143, 207, 234]
	mutPosV4 = [15, 37, 167, 185, 209]
	mutPosV5 = [60, 71, 81, 107, 125, 127, 133, 235]

	if sampleNum in a:
		mutPos = mutPosV3
	elif sampleNum in b:
		mutPos = mutPosV1
	elif sampleNum in c:
		mutPos = mutPosV2
	elif sampleNum in d:
		mutPos = mutPosV4
	elif sampleNum in e:
		mutPos = mutPosV5
	else:
		mutPos = list(set(mutPosV1) | set(mutPosV2) | set(mutPosV3) | set(mutPosV4) | set(mutPosV5))





	# print(mutPos)

	mutPosLenient = []
	# print(mutPosLenient)
	for i in mutPos:
		mutPosLenient.append(i)
		if not( (i+1) in mutPos ):
			mutPosLenient.append(i+1)
		if not( (i-1) in mutPos ):
			mutPosLenient.append(i-1)
	# print(mutPosLenient)

	for line in file:
		r += 1
		
		# print(r)

		# if r%10 == 0:
		# 	print("%s: %d"%(data,r))
		if r%4 == 0:
			nameIn = line
		elif r%4 == 1:
			seqIn = line
		elif r%4 == 2:
			continue
		elif r%4 == 3 and r < limitReads:
			# print("4")
			try:
				qualIn = line

				HAorg = Seq(seqIn.strip(), generic_dna)
				HARCorg = HAorg.reverse_complement()
				HA = HAorg[0:40]
				HARC = HARCorg[0:40]

				alnF = pairwise2.align.localms(caF, HA, 1,-2, -2, -1)
				alnR = pairwise2.align.localms(caR, HARC, 1,-2, -2, -1)
				fal1,fal2, fscore, fbegin, fend = alnF[0]
				ral1,ral2, rscore, rbegin, rend = alnR[0]

				doAnalysis = 0

				if (fscore > 8 and rscore > 8) or (fscore > 10) or (rscore > 10):
					doAnalysis = 1
					HA = HAorg
					HARC = HARCorg
				else:
					fscoreOLD = fscore
					rscoreOLD = rscore

					alnF = pairwise2.align.localms(caR, HA, 1,-2, -2, -1)
					alnR = pairwise2.align.localms(caF, HARC, 1,-2, -2, -1)
					fal1,fal2, fscore, fbegin, fend = alnF[0]
					ral1,ral2, rscore, rbegin, rend = alnR[0]
					
					if (fscore > 8 and rscore > 8) or (fscore > 10) or (rscore > 10):
						doAnalysis = 1
						HA = HARCorg
						HARC = HAorg
					else:
						doAnalysis = 0
						print(data)
						print(nameIn)
						print('DID NOT DO ANALYSIS!!!!')
						print("fscoreOLD:\t%s\trscoreOLD:\t%s"%(fscoreOLD,rscoreOLD))
						print("fscore:\t%s\trscore:\t%s"%(fscore,rscore))


				
				# HA = Seq(seqIn.strip(), generic_dna)
				# HARC = HA.reverse_complement()

				# aln = pairwise2.align.localms(HA, HAref, 1,-1, -4, -1)
				# al1,al3, score, begin, end = aln[0]
				# # print(format_alignment(*aln[0]))
				# doAnalysis = 1

				# if score < 400:
				# 	HAtemp = HA
				# 	HA = HARC
				# 	HARC = HAtemp
					
				# if score < 400:
				# 	doAnalysis = 0

				if doAnalysis:

					HARCshort = HARC[0:120]
					#################################################################
					# Obtains the UMI for the sequence				                #
					#   This is only necessary for C/A plasmids, as p1 plasmids		#
					#    do not have UMI 											#
					#################################################################
					alnUMI = pairwise2.align.localms(UMI_SEQUENCE, HARCshort, 1, 0, -4, -4)

					# print(format_alignment(*alnUMI[0]))

					ualn_temp = format_alignment(*alnUMI[0])
					ind = [i for i, a in enumerate(ualn_temp) if a == '\n']
					ual2 = ualn_temp[ind[0]+1:ind[1]]
					ual1,ual3, uscore, ubegin, uend = alnUMI[0]

					ual1 = ual1[ubegin:uend]
					ual2 = ual2[ubegin:uend]
					ual3 = ual3[ubegin:uend]


					ind = [i for i, a in enumerate(ual1) if a == 'N']
					umi = ual3[ind[0]:(ind[len(ind)-1]+1)]
					umi = umi.replace("-","")

					#################################################################
					# Obtains the barcode for the specific mutations plasmid				                #
					#   This is only necessary for C/A plasmids, as p1 plasmids		#
					#    do not have barcodes 											#
					#################################################################
					alnBC = pairwise2.align.localms(BARCODE_SEQUENCE, HARCshort, 1, 0, -4, -4)

					# print(format_alignment(*alnBC[0]))

					bcaln_temp = format_alignment(*alnBC[0])
					ind = [i for i, a in enumerate(bcaln_temp) if a == '\n']
					bcal2 = ualn_temp[ind[0]+1:ind[1]]
					bcal1,bcal3, bcscore, bcbegin, bcend = alnBC[0]

					bcal1 = bcal1[bcbegin:bcend]
					bcal2 = bcal2[bcbegin:bcend]
					bcal3 = bcal3[bcbegin:bcend]


					ind = [i for i, a in enumerate(bcal1) if a == 'N']
					barcode = bcal3[ind[0]:(ind[len(ind)-1]+1)]
					barcode = barcode.replace("-","")














					# with open ("temptest.txt",'a') as tempOut:
					# 	tempOut.write('---------\n')
					# 	tempOut.write(str(HA))
					# 	tempOut.write('\n')
					aln = pairwise2.align.localms(HA, HAref, 1,-1, -4, -1)
					al1,al3, score, begin, end = aln[0]
					aln2 = format_alignment(*aln[0])

					# for i in (1:len(aln2))
					ind = [i for i, a in enumerate(aln2) if a == '\n']
					# print (ind)
					al2 = aln2[ind[0]+1:ind[1]]

					# with open ("temptest.txt",'a') as tempOut:
					# 	tempOut.write(al1)
					# 	tempOut.write('\n')
					# 	tempOut.write(al2)
					# 	tempOut.write('\n')
					# 	tempOut.write(al3)
					# 	tempOut.write('\n')


					al1 = al1[begin:end]
					al2 = al2[begin:end]
					al3 = al3[begin:end]
					


					try:
						if str(al2) in ntDict:
							HA2 = ntDict[str(al2)]
							print(data)
							print("USED DICTIONARY!!!")
						else:
							al2OLD = al2

							HA2 = ''
							skip = 0

							for i in range(0,len(al2)):
								if skip:
									skip = 0
									continue
								if al2[i] == '|':
									HA2 += al1[i]
								elif al2[i] == '.':
									HA2 += al1[i]
								elif al2[i] == ' ' and al1[i] == '-':
									HA2 += al3[i]
								elif al2[i] == ' ' and al3[i] == '-' and i < (len(al2) - 3):
									if al1[i+1] == al1[i+2]:
										HA2 += al1[i]
										skip = 1
										continue
									elif al1[i-2] == al1[i-1]:
										HA2 = HA2[0:len(HA2)-1]
										HA2 += al1[i]
									elif al1[i] == al1[i+1]:
										continue
					except:
						print(data)
						print("DICTIONARY ERROR!")



						# aln = pairwise2.align.localms(HA2, HAref, 1,-1, -3, -1)

						# aln2 = format_alignment(*aln[0])


						# # for i in (1:len(aln2))
						# ind = [i for i, a in enumerate(aln2) if a == '\n']
						# al2 = aln2[ind[0]+1:ind[1]]

						# al1,al3, score, begin, end = aln[0]


						# al1 = al1[begin:end]
						# al2 = al2[begin:end]
						# al3 = al3[begin:end]

						ntDict[str(al2OLD)] = HA2



					HA2seq = Seq(HA2,generic_dna)
					HAref_AA = HArefFINAL.translate()
					HA2_AA = HA2seq.translate()

					try:
						if str(HA2_AA) in aaDict:
							al1 = aaDict[str(HA2_AA)]['al1']
							al2 = aaDict[str(HA2_AA)]['al2']
							al3 = aaDict[str(HA2_AA)]['al3']
							print(data)
							print("USED DICTIONARY!!!")
						else:
							aln = pairwise2.align.localms(HA2_AA, HAref_AA, 1,-1, -10, -10)
							# print (format_alignment(*aln[0]))
							aln2 = format_alignment(*aln[0])

							ind = [i for i, a in enumerate(aln2) if a == '\n']
							al2 = aln2[ind[0]+1:ind[1]]
							al1,al3, score, begin, end = aln[0]

							al1 = al1[al3.find("M"):al3.find("*")]
							al2 = al2[al3.find("M"):al3.find("*")]
							al3 = al3[al3.find("M"):al3.find("*")]

							aaDict[str(HA2_AA)] = {'al1':al1,'al2':al2,'al3':al3} 
					except:
						print(data)
						print("DICTIONARY ERROR!")

					# with open ("temptest.txt",'a') as tempOut:
					# 	tempOut.write (al1)
					# 	tempOut.write('\n')
					# 	tempOut.write (al2)
					# 	tempOut.write('\n')
					# 	tempOut.write (al3)
					# 	tempOut.write('\n')
					# with open (data+"2.txt",'a') as diag:
					# 	diag.write(al1+"\n")
					# 	diag.write(al2+"\n")
					# 	diag.write(al3+"\n")


					# Generates the dictionary key. This is only for mutations that
					# are supposed to appear
					ind = [i for i, a in enumerate(al2) if a == '.']
					dictKey = ""
					if len(ind) == 0:
						dictKey = "WT"
						dictKeyALL = "WT"
						dictKeyLENIENT = "WT"
					else:
						for i in ind:
							# print('%s%d%s' %(al3[i],i+1,al1[i]))
							if (i+1) in mutPos:
								dictKey = dictKey + str(i+1) + "-"
						# print(dictKey[0:(len(dictKey)-1)])
						dictKey = dictKey[0:(len(dictKey)-1)]


						# Generates the dictionary key for ALL mutations
						dictKeyALL = ""
						for i in ind:
							dictKeyALL = dictKeyALL + str(i+1) + "-"
						# print(dictKey[0:(len(dictKey)-1)])
						dictKeyALL = dictKeyALL[0:(len(dictKeyALL)-1)]

						# Generates the dictionary for LENIENT mutations
						dictKeyLENIENT = ""
						for i in ind:
							# print('%s%d%s' %(al3[i],i+1,al1[i]))
							if (i+1) in mutPosLenient:
								dictKeyLENIENT = dictKeyLENIENT + str(i+1) + "-"
						# print(dictKey[0:(len(dictKey)-1)])
						dictKeyLENIENT = dictKeyLENIENT[0:(len(dictKeyLENIENT)-1)]
					


					# print("Barcode(%d):\t%s"%(len(barcode),barcode))
					# print("UMI(%d):\t%s"%(len(umi),umi))


					
					# print (resultsDict)
					if dictKey in resultsDict:

						resultsDict[dictKey]['count'] += 1
						if barcode in resultsDict[dictKey]:

							resultsDict[dictKey][barcode]['count'] += 1
							if umi in resultsDict[dictKey][barcode]:

								resultsDict[dictKey][barcode][umi]['count'] += 1
							else:
								resultsDict[dictKey][barcode][umi] = {'count':1}
						else:
							resultsDict[dictKey][barcode] = {'count':1}
							resultsDict[dictKey][barcode][umi] = {'count':1}
					else:
							resultsDict[dictKey] = {'count':1}
							resultsDict[dictKey][barcode] = {'count':1}
							resultsDict[dictKey][barcode][umi] = {'count':1}


					if dictKeyALL in resultsDictALL:

						resultsDictALL[dictKeyALL]['count'] += 1
						if barcode in resultsDictALL[dictKeyALL]:

							resultsDictALL[dictKeyALL][barcode]['count'] += 1
							if umi in resultsDictALL[dictKeyALL][barcode]:

								resultsDictALL[dictKeyALL][barcode][umi]['count'] += 1
							else:
								resultsDictALL[dictKeyALL][barcode][umi] = {'count':1}
						else:
							resultsDictALL[dictKeyALL][barcode] = {'count':1}
							resultsDictALL[dictKeyALL][barcode][umi] = {'count':1}
					else:
							resultsDictALL[dictKeyALL] = {'count':1}
							resultsDictALL[dictKeyALL][barcode] = {'count':1}
							resultsDictALL[dictKeyALL][barcode][umi] = {'count':1}


					if dictKeyLENIENT in resultsDictLENIENT:

						resultsDictLENIENT[dictKeyLENIENT]['count'] += 1
						if barcode in resultsDictLENIENT[dictKeyLENIENT]:

							resultsDictLENIENT[dictKeyLENIENT][barcode]['count'] += 1
							if umi in resultsDictLENIENT[dictKeyLENIENT][barcode]:

								resultsDictLENIENT[dictKeyLENIENT][barcode][umi]['count'] += 1
							else:
								resultsDictLENIENT[dictKeyLENIENT][barcode][umi] = {'count':1}
						else:
							resultsDictLENIENT[dictKeyLENIENT][barcode] = {'count':1}
							resultsDictLENIENT[dictKeyLENIENT][barcode][umi] = {'count':1}
					else:
							resultsDictLENIENT[dictKeyLENIENT] = {'count':1}
							resultsDictLENIENT[dictKeyLENIENT][barcode] = {'count':1}
							resultsDictLENIENT[dictKeyLENIENT][barcode][umi] = {'count':1}
			except:
				continue



	# print(resultsDict)
	# pp = pprint.PrettyPrinter(indent=4)
	# pp.pprint(resultsDict)
	ppp = pprint.PrettyPrinter(indent=4,stream=open('out.'+data,'a'))
	ppp.pprint(resultsDict)

	with open('out.'+data+'.pickle', 'wb') as handle:
		pickle.dump(resultsDict, handle)

	with open('out.'+data+'.ALL.pickle', 'wb') as handle:
		pickle.dump(resultsDictALL, handle)


	with open('out.'+data+'.LENIENT.pickle', 'wb') as handle:
		pickle.dump(resultsDictLENIENT, handle)

	# with open('out.'+data+'.pickle', 'rb') as handle:
	# 	b = pickle.load(handle)

	print(data)
	print(datetime.datetime.now())


				# with open ("data_"+data+"2.txt",'a') as outFile2:
				# 	outFile2.write("\n")
				# 	for i in ind:
				# 		outFile2.write('%s%d%s\n' %(al3[i],i+1,al1[i]))
				# 		oneD[i] += 1

	# with open (("data_"+data+".txt"),'w') as outFile:
	# 	for i in range(0,len(oneD)):
	# 		outFile.write(str(oneD[i])+"\n")
			# if oneD[i] > 0:
				# print ("%d: %d"%(i+1,oneD[i]))













# sample_barcode_pairs = []
# with open(sample_barcode_file) as file:
# 	for line in file:
# 		line = line.split('\t')
# 		sample = line[0]
# 		fwd_barcode = line[1]
# 		rev_barcode = line[2]
# 		sample_barcode_pairs.append([sample, fwd_barcode, rev_barcode])

jobs = [] # create list of jobs to be multithreaded
# jobs.append("p1_4_4.fastq")

if sys.argv[1] == "all":
	for j in range(0,10):
		for k in range(0,7):
#		for l in range(0,19):
			try:
				with open ("CA_%d_%d.fastq"%(j,k)):
					jobs.append("CA_%d_%d.fastq"%(j,k))
			except:
				print("CA_%d_%d.fastq does not exist!"%(j,k))
elif sys.argv[1] == "set2":
	jobs.append("CA_1_2.fastq")
	jobs.append("CA_2_2.fastq")
	jobs.append("CA_4_2.fastq")
	jobs.append("CA_5_2.fastq")
	jobs.append("CA_7_2.fastq")
	jobs.append("CA_8_2.fastq")

else:
	jobs.append(sys.argv[1])
# for sample_barcode_pair in sample_barcode_pairs:
# 	job = {}
# 	job["sample"] = sample_barcode_pair[0]
# 	job["fwd_barcode"] = sample_barcode_pair[1]
# 	job["rev_barcode"] = sample_barcode_pair[2]
# 	jobs.append(job)

pool = mp.Pool(processes=mp.cpu_count())
print("CPU count: %f"%(mp.cpu_count()))
print(jobs)
multithread = pool.map(mtAlign,jobs)

print(datetime.datetime.now())
