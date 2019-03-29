def isStructureResidue(r):
	if(r == 'H'):
		return(1)
	if(r == 'E'):
		return(2)
	return(0)

def assignAminoIdx(fragment, amino):
	fragIdx = []
	for i in range(len(fragment)):
		fragIdx.append(amino.index(fragment[i]))
	return(fragIdx)

def getNeighboringAminos(fragment, window, amino):
	sideLen = int(window/2) # length of each side
	neighbors = []
	neighbors += assignAminoIdx(fragment[0:sideLen], amino)
	neighbors += assignAminoIdx(fragment[(sideLen+1):window], amino)
	
	return(neighbors)

def getAvgNeighboringSecStructScore(fragment, window, featureScore, amino):
	sideLen = int(window/2) # length of each side
	totalScore = 0.0
	neighbors = []
	neighbors += assignAminoIdx(fragment[0:sideLen], amino)
	neighbors += assignAminoIdx(fragment[(sideLen+1):window], amino)
	
	for n in neighbors:
		totalScore += featureScore[n]
	
	return(totalScore/(window - 1))