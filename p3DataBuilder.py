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


def buildDataSet(aaSeq, ssSeq, features, training=1):
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

		#set += getNeighboringAminos(fragment, window, amino)
		set.append(fHelper.getAvgNeighboringSecStructScore(fragment, window, features[0], amino))
		set.append(fHelper.getAvgNeighboringSecStructScore(fragment, window, features[1], amino))
		tSample.append(set)
		
		if(training):
			label = fHelper.isStructureResidue(ssSeq[i])
			tLabel.append(label)
	
	return([tSample, tLabel])	

	
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