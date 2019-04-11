from graphviz import Digraph
import pickle

def getFitness(muts,setName):
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

	with open('fitness-'+setName+'.pickle', 'rb') as handle:
		fitnessDict = pickle.load(handle)
		fitness = fitnessDict[name]

	return fitness

def getColor(value,*,minV = -.12,maxV = .7,colorMin = 'ed6e93',colorMax = '40b43c'):


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
		for i in range(0,len(remaining)):
			# print("i: %d"%i)
			# print("currentNode: %s\tremaining: %s"%(currentNode,remaining))
			nextNode = currentNode + "-" + str(remaining[i])

			nextFitness = getFitness(nextNode,setName)
			currentFitness = getFitness(currentNode,setName)
			# print(nextNode)
			# print(nextFitness)

			

			if (nextFitness > currentFitness):
				arrowColor = 'black'
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
					graph = addNodes(nextNode,graph,nextSet,setName)
				# print("currentNode: %s\tremaining: %s"%(currentNode,remaining))

	return graph


a = [37, 71, 85, 124, 182, 209, 226]
b = [15, 55, 71, 72, 74, 92, 139, 190, 209, 235]
c = [7, 15, 36, 61, 124, 143, 207, 234]
d = [15, 37, 167, 185, 209]
ef = [60, 71, 81, 107, 125, 127, 133, 235]
setA = '3-8'
setB = '4-9'
setC = '2-7'
setD = '5-10'
setE = '6-11'
setF = '6-12'


m = a
currSet = setA



graph = Digraph(comment='Graph Title')

wtFitness = getFitness('WT',currSet)


# color will be : #f0e3ff
graph.node('WT','WT\n%.4f'%(wtFitness),fillcolor=getColor(wtFitness),style='filled')

graph = addNodesALL('WT',graph,m,currSet)

# print(graph.source)
graph.render('graph2ALL')
