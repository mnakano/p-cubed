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
print("Total residues: ", aminoHelper.totalResidue)
print("Total helix: ", aminoHelper.totalHelixResidue)
print(aHelper.aminoAcids)
print(aminoHelper.aminoCount)
print(aminoHelper.aminoHelixCount)
print("Helix: ", scores[0])
print("Sheet: ", scores[1])

##### build a training dataset #####
trainingSet = builder.buildTrainingDataset(data, features) # build the training set to predict if a given amino acid in a protein sequence will be part of an alpha helix or not.
trainingSamples = trainingSet[0] # first element of the list returned from 'buildTrainingSet' is an list of training samples.
trainingLabels = trainingSet[1] # second element is the list of labels for each training sample. '1' indicates that the amino acid was recorded as 'helix-forming', '0' indicates otherwise.

##### training the machine learning model #####
clf = svm.SVC(decision_function_shape='ovo', gamma='auto')
clf.fit(trainingSamples, trainingLabels)  

dump(clf, 'clf-td9078-500.joblib')

print('training complete')