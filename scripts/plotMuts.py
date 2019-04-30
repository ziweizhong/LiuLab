from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import Align
import numpy as np
# import multiprocessing as mp
import datetime
import matplotlib
import matplotlib.pyplot as plt
import sys
import random

print(datetime.datetime.now())


# def mtAlign(data):
namePartB = [',H,',',H,',',W,',',W,',',UL,',',UL,']



if sys.argv[1] == 'timepoints':
	for l in range(0,6):
		for m in range(0,8):
			f, subplt = plt.subplots(2, 3, sharey=True, sharex=True)
			namePartA = str(int(m/4+1))
			namePartC = str(int(m%4+1+(l%2)*4))
			title = namePartA+namePartB[l]+namePartC


			for j in range(0,3):
				for k in range(0,2):

					count = []

					fIndex = l*3+j
					rIndex = m*2+k


					# file = open("data_p1_"+str(fIndex+1)+"_"+str(jIndex+1)+".fastq.txt",'r',buffering=1000)

					try:
						file = open("data_p1_"+str(fIndex+1)+"_"+str(rIndex+1)+".fastq.txt",'r',buffering=1000)
						run = 1
					except:
						print("file open error!")
						run = 0

					if run:
						for line in file:
							count.append(int(line.strip()))

						


						pos = range(1,len(count)+1)
						count = np.array(count)
						pos = np.array(pos)

						countFreq = count/max(count)


						with open("p1_binary.txt","a") as binary:
							with open("p1_percent.txt","a") as percent:
								binary.write("%s,%d,"%(title,k*3+j+1))
								percent.write("%s,%d,"%(title,k*3+j+1))
								for o in range(0,len(count)):
									binary.write("%d,"%(int((count[o]/max(count) >= 0.25) and (max(count) >= 4))))
									try:
										percent.write("%.2f,"%(count[o]/max(count)))
									except:
										percent.write("0,")
								binary.write("\n")
								percent.write("\n")

						# fig,ax = plt.subplots()
						subplt[k][j].plot(pos,countFreq)
						highCount = max(count)


						subplt[k][j].set_title(highCount)
						for i, a in enumerate(countFreq):
							if a > (0.25):
								# print("greater than 0.25")
								subplt[k][j].annotate(pos[i],(pos[i], countFreq[i]+(-1*random.random()/10)),fontsize=5)

						# fig.savefig("pngFreq_%d_%d.png" %(j,k))

			
			# f.suptitle(title)
			# f.savefig("%s.png"%(title),dpi=800)
			# f.close()

if sys.argv[1] == 'plateaus':
	f, subplt = plt.subplots(2, 2, sharey=True, sharex=True)
	f1, subplt1 = plt.subplots(2, 2, sharey=True, sharex=True)

	for j in range(0,8):

		count = []

		try:
			file = open("data_p1_"+str(j+1)+"_19"+".fastq.txt",'r',buffering=1000)
			run = 1
		except:
			print("file open error!")
			run = 0

		if run:
			for line in file:
				count.append(int(line.strip()))

			pos = range(1,len(count)+1)
			count = np.array(count)
			pos = np.array(pos)

			countFreq = count/max(count)
			highCount = max(count)
			# fig,ax = plt.subplots()

			if j in [0, 1, 4, 5]:
				subplt[j%4][int(j/4)].plot(pos,countFreq)
				subplt[j%4][int(j/4)].set_title("%d %s"%(j,highCount))
				for i, a in enumerate(countFreq):
					if a > (0.25):
						# print("greater than 0.25")
						subplt[j%4][int(j/4)].annotate(pos[i],(pos[i], countFreq[i]+(-1*random.random()/10)),fontsize=5)
			else:
				subplt1[(j-2)%4][int(j/4)].plot(pos,countFreq)
				subplt1[(j-2)%4][int(j/4)].set_title("%d %s"%(j,highCount))
				for i, a in enumerate(countFreq):
					if a > (0.25):
						# print("greater than 0.25")
						subplt1[(j-2)%4][int(j/4)].annotate(pos[i],(pos[i], countFreq[i]+(-1*random.random()/10)),fontsize=5)

			with open("p1PLAT_binary.txt","a") as binary:
				with open("p1PLAT_percent.txt","a") as percent:
					binary.write("PLAT-%d,"%(j+1))
					percent.write("PLAT-%d,"%(j+1))
					for o in range(0,len(count)):
						binary.write("%d,"%(int((count[o]/max(count) >= 0.25) and (max(count) >= 4))))
						try:
							percent.write("%.2f,"%(count[o]/max(count)))
						except:
							percent.write("0,")
					binary.write("\n")
					percent.write("\n")

			# fig.savefig("pngFreq_%d_%d.png" %(j,k))

	# namePartA = str(int(m/4+1))
	# namePartC = str(int(m%4+1+(l%2)*4))
	# title = namePartA+namePartB[l]+namePartC

	# f.suptitle('Plateaus Set 1,2,5,6')
	# f.savefig("Plat_1256.png",dpi=800)
	# # f.close()

	# f1.suptitle('Plateaus Set 3,4,7,8')
	# f1.savefig("Plat_3478.png",dpi=800)
	# f1.close()











# sample_barcode_pairs = []
# with open(sample_barcode_file) as file:
# 	for line in file:
# 		line = line.split('\t')
# 		sample = line[0]
# 		fwd_barcode = line[1]
# 		rev_barcode = line[2]
# 		sample_barcode_pairs.append([sample, fwd_barcode, rev_barcode])

# jobs = [] # create list of jobs to be multithreaded
# # jobs.append("p1_4_4.fastq")
# if sys.argv[1] == "all":
# 	for j in range(1,20):
# 		for k in range(1,20):
# 			try:
# 				with open ("p1_%d_%d.fastq"%(j,k)):
# 					jobs.append("p1_%d_%d"%(j,k))
# 			except:
# 				print("p1_%d_%d,fastq does not exist!"%(j,k))
# else:
# 	jobs.append(sys.argv[1])


# for sample_barcode_pair in sample_barcode_pairs:
# 	job = {}
# 	job["sample"] = sample_barcode_pair[0]
# 	job["fwd_barcode"] = sample_barcode_pair[1]
# 	job["rev_barcode"] = sample_barcode_pair[2]
# 	jobs.append(job)

# pool = mp.Pool(processes=mp.cpu_count())
# print("CPU count: %f"%(mp.cpu_count()))
# print(jobs)
# multithread = pool.map(mtAlign,jobs)

print(datetime.datetime.now())