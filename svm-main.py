from sklearn import svm # scikit learn module for machine learning
import re # Python's regex module
import os # Python's directory navigation module
import pandas

##### functions #####
def isAlphaHelix(h):
	if(h == 'H'):
		return(1)
	return(0)

def isBetaSheet(e):
	if(e == 'E'):
		return(1)
	return(0)

def readFile(filename):
	with open(filename, encoding='utf-8') as f:
		lines = f.readlines()
	return(''.join(lines))
	
def readCSV(filename):
	df = pandas.read_csv(filename)
	return(df[df['has_nonstd_aa'] == False][['seq', 'sst3']])
	
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

def getAvgNeighboringCPHelix(fragment, window, cpHelix, amino):
	sideLen = int(window/2) # length of each side
	totalCPHelix = 0.0
	neighbors = []
	neighbors += assignAminoIdx(fragment[0:sideLen], amino)
	neighbors += assignAminoIdx(fragment[(sideLen+1):window], amino)
	
	for n in neighbors:
		totalCPHelix += cpHelix[n]
	
	return (totalCPHelix)/(window - 1)

def buildTrainingSet(dir, structType):
	amino = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
	cpHelix = []
	aminoCount = [0] * 20 # Ni
	aminoHelixCount = [0] * 20 # nij
	totalNumResidues = 0 # Nt
	totalNumHelixResidues = 0 # Nj
	
	filenames = os.listdir(dir)
	sequences = []
	ssSequences = []
	samples = []
	labels = []
	
	for filename in filenames:
		
		content = readFile(dir + filename)
		extracted = extractSequenceSet(content)
		seq = extracted[0]
		ssSeq = extracted[1]
		sequences.append(seq)
		ssSequences.append(ssSeq)
		
		seqList = list(seq)
		ssSeqList = list(ssSeq)
		
		totalNumResidues += len(seqList)
		totalNumHelixResidues += ssSeqList.count('H')
		
		for aa in amino:
			aminoCount[amino.index(aa)] += seqList.count(aa)
			
		for i, h in enumerate(ssSeqList):
			if(h == 'H'):
				aminoHelixCount[amino.index(seqList[i])] += 1
	
	print(totalNumResidues)
	print(totalNumHelixResidues)
	print(amino)
	print(aminoCount)
	print(aminoHelixCount)	
	
	freqHelix = totalNumHelixResidues / totalNumResidues # fj
	
	for aa in amino:
		freqHelixAmino = aminoHelixCount[amino.index(aa)] / aminoCount[amino.index(aa)] # fij
		cpHelix.append(freqHelixAmino / freqHelix) # Pij
	cpHelix.append(0.0) # Pij for '-'
	print(cpHelix)
	
	for i, seq in enumerate(sequences):

		learningSet = buildDataSet(seq, ssSequences[i], cpHelix, 'H')
		samples += learningSet[0]
		labels += learningSet[1]
	
	return([samples, labels, cpHelix])	

def extractSequenceSet(content):
	seqSet = []
	match = re.search("OrigSeq:(.*)\n", content)
	span = match.span()	
	seq = content[span[0]+8:span[1]-1]
	seqSet.append(seq)
	
	match = re.search("cons:(.*)\n", content)
	span = match.span()
	ssSeq = content[span[0]+5:span[1]-1]
	seqSet.append(ssSeq)
	
	return(seqSet)

def buildDataSet(aaSeq, ssSeq, cpHelix, structType):
	amino = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '-']
	kdValues = [1.8, -4.0, -3.5, -3.5, 2.5, -3.5,  -3.5, 0.4, -3.2, 4.5, 3.8, -3.9, 1.9, 2.8, -1.6, -0.8, -0.7, -0.9, -1.3, 4.2]
	
	tSample = []
	tLabel = []
	label = -1
	
	window = 7
	sideLen = int(window/2)
	left = 0
	right = sideLen + 1
	
	for i, aa in enumerate(aaSeq):

		set = [amino.index(aa)]
		#set.append(kdValues[amino.index(aa)])
		set.append(cpHelix[amino.index(aa)])
		
		fragment = ''
		if(i < sideLen):
			fragment = '-' * (sideLen - left) + aaSeq[(i-left):(i+right)]
			left += 1 
		elif(i > (len(aaSeq) - (sideLen + 1))):
			right -= 1
			fragment = aaSeq[(i-left):(i+right)] + '-' * ((sideLen + 1) - right)
		else:
			fragment = aaSeq[(i-left):(i+right)]

		#set += getNeighboringAminos(fragment, window, amino)
		
		set.append(getAvgNeighboringCPHelix(fragment, window, cpHelix, amino))
		
		tSample.append(set)
		
		if(structType == 'H'):
			label = isAlphaHelix(ssSeq[i])
		else:
			label = isBetaSheet(ssSeq[i])
		
		tLabel.append(label)
	
	return([tSample, tLabel])
	
def buildPredictedSequence(result, structType):
	predSeq = ''
	for r in result:
		if(r == 1):
			predSeq += structType
		else:
			predSeq += '-'
	return(predSeq)

##### end of functions #####

dir = 'rs126\\' # a directory where the training data files are stored.

trainingSet = buildTrainingSet(dir, 'H') # build the training set to predict if a given amino acid in a protein sequence will be part of an alpha helix or not.
trainingSamples = trainingSet[0] # first element of the list returned from 'buildTrainingSet' is an list of training samples.
trainingLabels = trainingSet[1] # second element is the list of labels for each training sample. '1' indicates that the amino acid was recorded as 'helix-forming', '0' indicates otherwise.
cpHelix = trainingSet[2]

#print(trainingSamples)
#print(trainingLabels)

##### training the machine learning model #####
clf = svm.SVC(gamma=0.5)
#clf = svm.SVC(C=1.5, gamma=1.0)
clf.fit(trainingSamples, trainingLabels) 
##### 'supported vector machine' model is used ##### 


##### predicting helix formation for a given set of test samples
#filenames = os.listdir('helix\\')
#for file in filenames:
	
##### preparing a test dataset #####
testContent = readFile('1lmb3.concise') # this file is used as a test data ('unknown' protein sequence to the machine).
testExtracted = extractSequenceSet(testContent)
testSeq = testExtracted[0]
testSecStructSeq = testExtracted[1]
testSet = buildDataSet(testSeq, testSecStructSeq, cpHelix, 'H')


result = clf.predict(testSet[0])

print(result)
predictedSecStructSeq = buildPredictedSequence(result, 'H')

print(testSeq)
print(testSecStructSeq)
print(predictedSecStructSeq)	
