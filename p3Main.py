from sklearn import svm # scikit learn module for machine learning
import pandas
import p3DataReader as reader
import p3DataBuilder as builder
import aminoHelper as aHelper
import featureHelper as fHelper

features = []

rs126 = reader.getRS126Dataset('rs126\\')
aminoHelper = aHelper.AminoHelper()
helixScores = aminoHelper.calculateHelixScore(rs126)

features.append(helixScores)
print(aminoHelper.totalResidue)
print(aminoHelper.totalHelixResidue)
print(aHelper.aminoAcids)
print(aminoHelper.aminoCount)
print(aminoHelper.aminoHelixCount)
print(helixScores)

##### build a training dataset #####
trainingSet = builder.buildTrainingDataset(rs126, features) # build the training set to predict if a given amino acid in a protein sequence will be part of an alpha helix or not.
trainingSamples = trainingSet[0] # first element of the list returned from 'buildTrainingSet' is an list of training samples.
trainingLabels = trainingSet[1] # second element is the list of labels for each training sample. '1' indicates that the amino acid was recorded as 'helix-forming', '0' indicates otherwise.

##### training the machine learning model #####
clf = svm.SVC(gamma=0.5)
clf = svm.SVC(C=1.5, gamma=1.0)
clf.fit(trainingSamples, trainingLabels)  
	
##### preparing a test dataset #####
testContent = reader.readFile('1lmb3.concise') # this file is used as a test data ('unknown' protein sequence to the machine).
testExtracted = reader.extractSequenceSet(testContent)
testSeq = testExtracted[0]
testSecStructSeq = testExtracted[1]
testSet = builder.buildDataSet(testSeq, testSecStructSeq, features, 'H')

##### predict the secondary structure #####
result = clf.predict(testSet[0])
predictedSecStructSeq = builder.buildPredictedSequence(result, 'H')

##### print the results
print(testSeq)
print(testSecStructSeq)
print(predictedSecStructSeq)	

