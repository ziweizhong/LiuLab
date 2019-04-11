from graphviz import Digraph
def getName(muts):
	sepPos = [i for i, a in enumerate(str(muts)) if a == '-']
	if len(sepPos) == 0:
		return ''
	else:
		mutsMtx = []
		# mutsMtx.append(int(muts[0:sepPos[0]]))
		for i in range(0,len(sepPos)-1):
			mutsMtx.append(int(muts[(sepPos[i]+1):sepPos[i+1]]))
		mutsMtx.append(int(muts[(sepPos[len(sepPos)-1]+1):len(muts)]))
		
		for i in mutsMtx:
			try:
				mutInd = mutPos.index(i)
				mutPresence[mutInd] = 1
			except:
				continue

		name = ''
		for i in mutsMtx:
			name += (str(i) + ".")

		name = name[0:len(name)-1]

		return name


def addNodes(currentNode,graph,remain):
	remaining = remain.copy()
	# print("remaining: %s"%remaining)
	# print("len(remaning): %d"%(len(remaining)))
	if len(remaining) > 0:
		# print("len(remaning) 2: %d"%(len(remaining)))
		for i in range(0,len(remaining)):
			# print("i: %d"%i)
			# print("currentNode: %s\tremaining: %s"%(currentNode,remaining))
			name = currentNode + "-" + str(remaining[i])

			name2 = getName(name)
			print(name)
			print(name2)
			print(type(name2))

			fitness = 2
			graph.node(name,(name+'\n'+str(fitness)))
			arrowColor = 'red'
			graph.edge(currentNode,name,color=arrowColor)
			nextSet = remaining.copy()
			nextSet.remove(remaining[i])
			# print("name %s\tnextSet %s\tlen(nextSet) %d"%(name,nextSet,len(nextSet)))
			if len(nextSet) > 0:
				graph = addNodes(name,graph,nextSet)
			# print("currentNode: %s\tremaining: %s"%(currentNode,remaining))

	return graph


a = [1,2,3]

graph = Digraph(comment='Graph Title')

graph.node('WT','WT')

graph = addNodes('WT',graph,a)

# print(graph.source)
graph.render('graph')
