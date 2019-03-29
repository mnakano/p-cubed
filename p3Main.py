from sklearn import svm # scikit learn module for machine learning
from joblib import dump, load # module used to save and load trained svm models
import p3DataReader as reader
import p3DataBuilder as builder
import aminoHelper as aHelper
import featureHelper as fHelper

features = []

#data = reader.getRS126Dataset('rs126\\')
data = reader.getTD9078Dataset('data.csv', 500)
aminoHelper = aHelper.AminoHelper()
scores = aminoHelper.calculateSecStructScore(data)
features += scores

##### preparing a test dataset #####
testContent = reader.readFile('1mcpl.concise') # this file is used as a test data ('unknown' protein sequence to the machine).
testExtracted = reader.extractSequenceSet(testContent)
testSeq = testExtracted[0]
testSecStructSeq = testExtracted[1]
testSet = builder.buildDataSet(testSeq, testSecStructSeq, features, 0)

##### predict the secondary structure #####
clf = load('clf-td9078-500.joblib')
#clf = load('clf-rs126.joblib')
result = clf.predict(testSet[0])
print(result)
predictedSecStructSeq = builder.buildPredictedSequence(result)

##### print the results #####
print("prediction complete")
print(testSeq)
print(testSecStructSeq)
print(predictedSecStructSeq)	

