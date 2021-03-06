'''
p3DataBuilder.py
The library contains functions relevant to building training and test datasets.
'''

import aminoHelper as aHelper
import featureHelper as fHelper

def buildTrainingDataset(dataset, features):
	samples = []
	labels = []
	for data in dataset:

		learningSet = buildDataSet(data[0], data[1], features)
		samples += learningSet[0]
		labels += learningSet[1]
	
	return([samples, labels])


def buildDataSet(aaSeq, ssSeq, features):
	amino = aHelper.aminoAcids
	amino.append('-')
	kdValues = aHelper.kdValues
	
	tSample = []
	tLabel = []
	label = -1
	
	window = 7
	sideLen = int(window/2)
	left = 0
	right = sideLen + 1
	
	for i, aa in enumerate(aaSeq):

		set = [amino.index(aa)]
		set.append(features[0][amino.index(aa)])
		set.append(features[1][amino.index(aa)])
		
		fragment = ''
		if(i < sideLen):
			fragment = '-' * (sideLen - left) + aaSeq[(i-left):(i+right)]
			left += 1 
		elif(i > (len(aaSeq) - (sideLen + 1))):
			right -= 1
			fragment = aaSeq[(i-left):(i+right)] + '-' * ((sideLen + 1) - right)
		else:
			fragment = aaSeq[(i-left):(i+right)]

		set.append(fHelper.getAvgNeighboringSecStructScore(fragment, window, features[0], amino))
		set.append(fHelper.getAvgNeighboringSecStructScore(fragment, window, features[1], amino))
		set.append(fHelper.getAvgNeighboringSecStructScore(fragment, window, kdValues, amino))
		tSample.append(set)
		
		# Build the labels to be used for training.
		label = fHelper.isStructureResidue(ssSeq[i])
		tLabel.append(label)
	
	return([tSample, tLabel])

def buildTestDataset(aaSeq, features):
	amino = aHelper.aminoAcids
	amino.append('-')
	kdValues = aHelper.kdValues
	
	testSample = []
	
	window = 7
	sideLen = int(window/2)
	left = 0
	right = sideLen + 1
	
	for i, aa in enumerate(aaSeq):

		set = [amino.index(aa)]
		set.append(features[0][amino.index(aa)])
		set.append(features[1][amino.index(aa)])
		
		fragment = ''
		if(i < sideLen):
			fragment = '-' * (sideLen - left) + aaSeq[(i-left):(i+right)]
			left += 1 
		elif(i > (len(aaSeq) - (sideLen + 1))):
			right -= 1
			fragment = aaSeq[(i-left):(i+right)] + '-' * ((sideLen + 1) - right)
		else:
			fragment = aaSeq[(i-left):(i+right)]

		set.append(fHelper.getAvgNeighboringSecStructScore(fragment, window, features[0], amino))
		set.append(fHelper.getAvgNeighboringSecStructScore(fragment, window, features[1], amino))
		set.append(fHelper.getAvgNeighboringSecStructScore(fragment, window, kdValues, amino))
		testSample.append(set)
	
	return(testSample)
	
def buildPredictedSequence(result):
	predSeq = ''
	for r in result:
		if(r == 1):
			predSeq += 'H'
		elif(r == 2):
			predSeq += 'E'
		else:
			predSeq += '-'
	return(predSeq)