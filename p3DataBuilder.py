import aminoHelper as aHelper
import featureHelper as fHelper

def buildTrainingDataset(dataset, features):
	samples = []
	labels = []
	for data in dataset:

		learningSet = buildDataSet(data[0], data[1], features, 'H')
		samples += learningSet[0]
		labels += learningSet[1]
	
	return([samples, labels])


def buildDataSet(aaSeq, ssSeq, features, structType):
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
		#set.append(kdValues[amino.index(aa)])
		set.append(features[0][amino.index(aa)])
		
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
		
		set.append(fHelper.getAvgNeighboringHelixScore(fragment, window, features[0], amino))
		
		tSample.append(set)
		
		if(structType == 'H'):
			label = fHelper.isAlphaHelix(ssSeq[i])
		else:
			label = fHelper.isBetaSheet(ssSeq[i])
		
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