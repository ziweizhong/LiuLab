from graphviz import Digraph
import pickle
import pprint
import numpy as np
from itertools import permutations

def getName(muts):
	sepPos = [i for i, a in enumerate(str(muts)) if a == '-']
	if len(sepPos) == 0:
		name = 'WT'
	else:
		mutsMtx = []
		# mutsMtx.append(int(muts[0:sepPos[0]]))
		for i in range(0,len(sepPos)-1):
			mutsMtx.append(int(muts[(sepPos[i]+1):sepPos[i+1]]))
		mutsMtx.append(int(muts[(sepPos[len(sepPos)-1]+1):len(muts)]))
		
		mutsMtx.sort()

		name = ''
		for i in mutsMtx:
			name += (str(i) + "-")

		name = name[0:len(name)-1]

	return name



def cleanName(muts):
	sepPos = [i for i, a in enumerate(str(muts)) if a == '-']
	if len(sepPos) == 0:
		name = muts
	else:
		mutsMtx = []
		mutsMtx.append(int(muts[0:sepPos[0]]))
		for i in range(0,len(sepPos)-1):
			mutsMtx.append(int(muts[(sepPos[i]+1):sepPos[i+1]]))
		mutsMtx.append(int(muts[(sepPos[len(sepPos)-1]+1):len(muts)]))
		
		mutsMtx.sort()

		name = ''
		for i in mutsMtx:
			name += (str(i) + "-")

		name = name[0:len(name)-1]

	return name

def getFitness(muts,setName):
	name = getName(muts)

	with open('fitness-'+setName+'.pickle', 'rb') as handle:
		fitnessDict = pickle.load(handle)
		try:
			fitness = fitnessDict[name]
		except:
			fitness = fitnessDict['WT']

	return fitness

def getFitness2(muts,setName):
	name = muts

	with open('fitness-'+setName+'.pickle', 'rb') as handle:
		fitnessDict = pickle.load(handle)
		try:
			fitness = fitnessDict[name]
		except:
			fitness = fitnessDict['WT']

	return fitness


# def getColor(value,*,minV = -.0,maxV = .286,colorMin = 'ed6e93',colorMax = '40b43c'):
def getColor(value,*,minV = -.0,maxV = .3,colorMin = 'ece7f2',colorMax = '2b8cbe'):


	rMin = int(colorMin[0:2],16)
	gMin = int(colorMin[2:4],16)
	bMin = int(colorMin[4:6],16)
	rMax = int(colorMax[0:2],16)
	gMax = int(colorMax[2:4],16)
	bMax = int(colorMax[4:6],16)

	if value > maxV:
		value = maxV
	if value < minV:
		value = minV

	# print(value)

	percent = (value - minV)/(maxV - minV)

	# print(percent)

	rSet = int(percent * (rMax - rMin) + rMin)
	gSet = int(percent * (gMax - gMin) + gMin)
	bSet = int(percent * (bMax - bMin) + bMin)

	# print('%d %d %d'%(rSet,gSet,bSet))

	colorHex = '#' + hex(rSet)[2:4] + hex(gSet)[2:4] + hex(bSet)[2:4]

	return colorHex


def addNodesALL(currentNode,graph,remain,setName):
	remaining = remain.copy()
	# print("remaining: %s"%remaining)
	# print("len(remaning): %d"%(len(remaining)))
	if len(remaining) > 0:
		# print("len(remaning) 2: %d"%(len(remaining)))
		for i in range(0,len(remaining)):
			# print("i: %d"%i)
			# print("currentNode: %s\tremaining: %s"%(currentNode,remaining))
			nextNode = currentNode + "-" + str(remaining[i])

			nextFitness = getFitness(nextNode,setName)
			currentFitness = getFitness(currentNode,setName)
			# print(nextNode)
			# print(nextFitness)

			

			
			arrowColor = 'black'
			if nextFitness < currentFitness:
				arrowColor = 'red'

			arrowWidth = int((nextFitness - currentFitness)*20)
			if arrowWidth > 15:
				arrowWidth = 15
			elif arrowWidth < 1:
				arrowWidth  = 1
			arrowWidth = str(arrowWidth)

			graph.node(nextNode,(nextNode[3:len(nextNode)]+'\n%.4f'%(nextFitness)),fillcolor=getColor(nextFitness),style='filled')
			graph.edge(currentNode,nextNode,color=arrowColor,penwidth=arrowWidth)

			nextSet = remaining.copy()
			nextSet.remove(remaining[i])
			# print("nextNode %s\tnextSet %s\tlen(nextSet) %d"%(nextNode,nextSet,len(nextSet)))
			if len(nextSet) > 0:
				graph = addNodesALL(nextNode,graph,nextSet,setName)
			# print("currentNode: %s\tremaining: %s"%(currentNode,remaining))

	return graph

def addNodes(currentNode,graph,remain,setName):
	remaining = remain.copy()
	# print("remaining: %s"%remaining)
	# print("len(remaning): %d"%(len(remaining)))
	if len(remaining) > 0:
		# print("len(remaning) 2: %d"%(len(remaining)))
		isEnd = np.zeros(len(remaining))

		for i in range(0,len(remaining)):
			# print("i: %d"%i)
			# print("currentNode: %s\tremaining: %s"%(currentNode,remaining))
			nextNode = currentNode + "-" + str(remaining[i])

			nextFitness = getFitness(nextNode,setName)
			currentFitness = getFitness(currentNode,setName)
			# print(nextNode)
			# print(nextFitness)

			

			if (nextFitness > currentFitness):
				isEnd[i] = 1
				arrowColor = 'black'
				arrowWidth = int((nextFitness - currentFitness)*20)
				if arrowWidth > 15:
					arrowWidth = 15
				elif arrowWidth < 1:
					arrowWidth  = 1
				arrowWidth = str(arrowWidth)
				arrowWidth = '1'

				simpleNextName = getName(nextNode)
				simpleCurrentName = getName(currentNode)
				graph.node(nextNode,(nextNode[3:len(nextNode)]+'\n%.4f'%(nextFitness)),fillcolor=getColor(nextFitness),style='filled')
				# graph.node(simpleNextName)
				graph.edge(currentNode,nextNode,color=arrowColor,penwidth=arrowWidth)
				# graph.edge(simpleCurrentName,simpleNextName)

				nextSet = remaining.copy()
				nextSet.remove(remaining[i])
				# print("nextNode %s\tnextSet %s\tlen(nextSet) %d"%(nextNode,nextSet,len(nextSet)))
				# if len(nextSet) > 0:
				graph = addNodes(nextNode,graph,nextSet,setName)
				# print("currentNode: %s\tremaining: %s"%(currentNode,remaining))
		if isEnd.max() < 1:
			#CHECK THE SYNTAX
			# print("AAA%s"%currentNode)
			currentNodeName = getName(currentNode)
			if currentNodeName in countDict:
				countDict[currentNodeName] += 1
			else:
				countDict[currentNodeName] = 1
	else:
		# print("BBB%s"%currentNode)
		currentNodeName = getName(currentNode)
		if currentNodeName in countDict:
			countDict[currentNodeName] += 1
		else:
			countDict[currentNodeName] = 1

	return graph

def addNodes2(currentNode,graph,remain,setName,addedDict):
	remaining = remain.copy()
	# print(currentNode)
	# print(remaining)

	if len(remaining) > 0:
		# print(remaining)
		edgeWeights = np.zeros(len(remaining))
		for i in range(0,len(remaining)):
			# print(i)
			# print(currentNode)
			# print(remaining)


			nextNode = currentNode + "-" + str(remaining[i])

			nextFitness = getFitness(nextNode,setName)
			currentFitness = getFitness(currentNode,setName)
				

			if (nextFitness > currentFitness):
				simpleNextName = getName(nextNode)
				simpleCurrentName = getName(currentNode)
				if not(getName(nextNode) in addedDict):
					addedDict[getName(nextNode)] = {'ChanceValue':0}

					graph.node(simpleNextName,(simpleNextName+'\n%.4f'%(nextFitness)),fillcolor=getColor(nextFitness),style='filled')

					nextSet = remaining.copy()
					nextSet.remove(remaining[i])
					# print(nextNode)
					# print(nextSet)
					# print(addedDict)
					addNodes2(nextNode,graph,nextSet,setName,addedDict)
					# print('After addNodes2')
					# print(nextNode)
					# print(nextSet)

				edgeWeights[i] = nextFitness - currentFitness
				edgeWeights[i] = nextFitness
				

			# print(edgeWeights)
			# print(i)

		for i in range(0,len(remaining)):
			if edgeWeights[i] != 0:
				nextNode = getName(currentNode + "-" + str(remaining[i]))

				arrowWidth = str(5 * edgeWeights[i] / np.sum(edgeWeights))
				arrowColor = 'black'
				graph.edge(getName(currentNode),nextNode,color=arrowColor,penwidth=arrowWidth)

				addedDict[getName(currentNode)][nextNode] = edgeWeights[i] / np.sum(edgeWeights)

				# addedDict[nextNode]['ChanceValue'] += addedDict[getName(currentNode)]['ChanceValue'] * edgeWeights[i] / np.sum(edgeWeights)


#Start from arbitrary point, not just WT
def addNodes3(currentNode,graph,remain,setName):
	remaining = remain.copy()
	# print("remaining: %s"%remaining)
	# print("len(remaning): %d"%(len(remaining)))
	if len(remaining) > 0:
		# print("len(remaning) 2: %d"%(len(remaining)))
		isEnd = np.zeros(len(remaining))

		for i in range(0,len(remaining)):
			# print("i: %d"%i)
			# print("currentNode: %s\tremaining: %s"%(currentNode,remaining))
			nextNode = currentNode + "-" + str(remaining[i])

			print("currentNode: %s\nNextNode: %s\ngetName:%s"%(currentNode,nextNode,cleanName(currentNode)))

			nextFitness = getFitness2(cleanName(nextNode),setName)
			currentFitness = getFitness2(cleanName(currentNode),setName)
			# print(nextNode)
			# print(nextFitness)

			print("currentFitness: %s\nNextFitness: %s\n"%(currentFitness,nextFitness))

			if (nextFitness > currentFitness):
				isEnd[i] = 1
				arrowColor = 'black'
				arrowWidth = int((nextFitness - currentFitness)*20)
				if arrowWidth > 15:
					arrowWidth = 15
				elif arrowWidth < 1:
					arrowWidth  = 1
				arrowWidth = str(arrowWidth)
				arrowWidth = '1'


				print(nextNode)

				graph.node(nextNode,nextNode+'\n%.4f'%(nextFitness),fillcolor=getColor(nextFitness),style='filled')
				# graph.node(simpleNextName)
				graph.edge(currentNode,nextNode,color=arrowColor,penwidth=arrowWidth)
				# graph.edge(simpleCurrentName,simpleNextName)

				nextSet = remaining.copy()
				nextSet.remove(remaining[i])

				graph = addNodes3(nextNode,graph,nextSet,setName)

		if isEnd.max() < 1:
			#CHECK THE SYNTAX
			# print("AAA%s"%currentNode)
			currentNodeName = getName(currentNode)
			if currentNodeName in countDict:
				countDict[currentNodeName] += 1
			else:
				countDict[currentNodeName] = 1
	else:
		# print("BBB%s"%currentNode)
		currentNodeName = getName(currentNode)
		if currentNodeName in countDict:
			countDict[currentNodeName] += 1
		else:
			countDict[currentNodeName] = 1

	return graph



a = [37, 71, 85, 124, 182, 209, 226]
b = [15, 55, 71, 72, 74, 92, 139, 190, 209, 235]
c = [7, 15, 36, 61, 124, 143, 207, 234]
d = [15, 37, 167, 185, 209]
e = [60, 71, 81, 107, 125, 127, 133, 235]
f = [60, 71, 81, 107, 125, 127, 133, 235]

mutSet = []
mutSet.append(a)
mutSet.append(b)
mutSet.append(c)
mutSet.append(d)
mutSet.append(e)
mutSet.append(f)
setA = '3-8'
setB = '4-9'
setC = '2-7'
setD = '5-10'
setE = '6-11'
setF = '6-12'
setMtx = [setA, setB, setC, setD, setE, setF]

countDict = {}
addedDict = {}
m = b
currSet = setB




# for i in range(0,6):
# 	m = mutSet[i]
# 	currSet = setMtx[i]
# 	graph = Digraph(comment='Graph Title')

# 	wtFitness = getFitness('WT',currSet)


# 	# color will be : #f0e3ff
# 	graph.node('WT','WT\n%.4f'%(wtFitness),fillcolor=getColor(wtFitness),style='filled')

# 	graph = addNodes('WT',graph,m,currSet)

# 	# print(graph.source)
# 	graph.render(currSet)

# 	pp = pprint.PrettyPrinter(indent=4)
# 	pp.pprint(countDict)

currSet = setE
graph = Digraph(comment='Graph Title')
wtFitness = getFitness2('71',currSet)
addedDict['71'] = {'ChanceValue':1}
# color will be : #f0e3ff
graph.node('71','71\n%.4f'%(wtFitness),fillcolor=getColor(wtFitness),style='filled')
addNodes3('71',graph,[60, 81, 107, 125, 127, 133, 235],setE)

graph.render('6-11_71')





# graph2 = Digraph(comment='Graph Title')
# graph2.node('WT','WT\n%.4f'%(addedDict['WT']['ChanceValue']),fillcolor=getColor(wtFitness),style='filled')

# for k in range(0,5):
# 	if k == 0:
# 		for item in addedDict['WT']:
# 			if item in addedDict:
# 				# graph2.node(item,item+'\n%.3f'%(getFitness2(item,setD)),fillcolor=getColor(getFitness2(item,setD)),style='filled')
# 				addedDict[item]['ChanceValue'] += addedDict['WT']['ChanceValue'] * addedDict['WT'][item]
# 				arrowWidth = addedDict['WT']['ChanceValue'] * addedDict['WT'][item]*15
# 				AW2 = "%.3f"%(arrowWidth/15)
# 				if arrowWidth < 0.25:
# 					arrowWidth = str(0.25)
# 				else:
# 					arrowWidth = str(arrowWidth)
# 				graph2.edge('WT',item,color='black',penwidth=arrowWidth,label=AW2)
# 	else:
# 		perm = permutations(d,k)
# 		for i in perm:
# 			m2 = ''
# 			for j in i:
# 				m2 += (str(j) + '-')
# 			m2 = m2[0:len(m2)-1]

# 			if m2 in addedDict:
# 				for item in addedDict[m2]:
# 					if item in addedDict:
# 						# graph2.node(item,item+'\n%.3f'%(getFitness2(item,setD)),fillcolor=getColor(getFitness2(item,setD)),style='filled')
# 						addedDict[item]['ChanceValue'] += addedDict[m2]['ChanceValue'] * addedDict[m2][item]
# 						arrowWidth = addedDict[m2]['ChanceValue'] * addedDict[m2][item]*15
# 						AW2 = "%.3f"%(arrowWidth/15)
# 						if arrowWidth < 0.25:
# 							arrowWidth = str(0.25)
# 						else:
# 							arrowWidth = str(arrowWidth)
# 						graph2.edge(m2,item,color='black',penwidth=arrowWidth,label=AW2)

# for item in addedDict:
# 	graph2.node(item,item+'\n%.3f'%(addedDict[item]['ChanceValue']),fillcolor=getColor(getFitness2(item,setD)),style='filled')

# graph2.render('5-10_2')
# graph2.format = 'svg'
# graph2.render('5-10_2')
# print(graph2.source)

# pp = pprint.PrettyPrinter(indent=10)
# pp.pprint(addedDict)