def isAlphaHelix(h):
	if(h == 'H'):
		return(1)
	return(0)

def isBetaSheet(e):
	if(e == 'E'):
		return(1)
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

def getAvgNeighboringHelixScore(fragment, window, helixScore, amino):
	sideLen = int(window/2) # length of each side
	totalHelixScore = 0.0
	neighbors = []
	neighbors += assignAminoIdx(fragment[0:sideLen], amino)
	neighbors += assignAminoIdx(fragment[(sideLen+1):window], amino)
	
	for n in neighbors:
		totalHelixScore += helixScore[n]
	
	return(totalHelixScore/(window - 1))