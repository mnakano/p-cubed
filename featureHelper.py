'''
featureHelper.py
Contains functions used to calculate and evaluate feature-related values for the training dataset.
'''

'''
Determines if a resiude is either a helix- or sheet-forming amino acid.
'''
def isStructureResidue(r):
	if(r == 'H'):
		return(1)
	if(r == 'E'):
		return(2)
	return(0)

'''
Takes in a fragment of amino acid sequence, and assigns a numerical value 
corresponding to each of the amino acids in the sequence fragment.
'''
def assignAminoIdx(fragment, amino):
	fragIdx = []
	for i in range(len(fragment)):
		fragIdx.append(amino.index(fragment[i]))
	return(fragIdx)

'''
Takes in a fragment of amino acid sequence, window length and a list of 20 amino acids.
Returns a list of numerical values correspnding to amino acids that are surrounding the amino acid that is being evaluated.
'''
def getNeighboringAminos(fragment, window, amino):
	sideLen = int(window/2) # length of each side
	neighbors = []
	neighbors += assignAminoIdx(fragment[0:sideLen], amino)
	neighbors += assignAminoIdx(fragment[(sideLen+1):window], amino)
	
	return(neighbors)

'''
Takes in a fragment of amino acid sequence, window length, list of feature score for 20 amino acids, and a list of 20 amino acids.
Calculates average feature score of the neighboring amino acids that are surrounding the amio acid that is being evaluated.
'''
def getAvgNeighboringSecStructScore(fragment, window, featureScore, amino):
	sideLen = int(window/2) # length of each side
	totalScore = 0.0
	neighbors = []
	neighbors += assignAminoIdx(fragment[0:sideLen], amino)
	neighbors += assignAminoIdx(fragment[(sideLen+1):window], amino)
	
	for n in neighbors:
		totalScore += featureScore[n]
	
	return(totalScore/(window - 1))