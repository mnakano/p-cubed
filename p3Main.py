from sklearn import svm # scikit learn module for machine learning
'''
p3Main.py
This is the main progam which predicts secondary structure of unknonw amino acid sequences using a pre-trained model.
'''

from joblib import dump, load # module used to save and load trained svm models
import p3DataReader as reader
import p3DataBuilder as builder
import aminoHelper as aHelper
import featureHelper as fHelper

features = []

##### Calculate feature scores used for the test data.
#data = reader.getRS126Dataset('rs126\\')
data = reader.getTD9078Dataset('td9078.csv', 500)
aminoHelper = aHelper.AminoHelper()
scores = aminoHelper.calculateSecStructScore(data)
features += scores

##### preparing a test dataset #####
# Uncomment this block to perform prediction with a single test data.
'''
testContent = reader.readFile('1mcpl.concise') # this file is used as a test data ('unknown' protein sequence to the machine).
testExtracted = reader.extractSequenceSet(testContent)
testSeq = testExtracted[0]
testSecStructSeq = testExtracted[1]
testSet = builder.buildDataSet(testSeq, testSecStructSeq, features, 0)

##### predict the secondary structure #####
clf = load('clf-td9078-500.joblib')
#clf = load('models\\clf-rs126.joblib')
result = clf.predict(testSet[0])
print(result)
predictedSecStructSeq = builder.buildPredictedSequence(result)

##### print the results #####
print("prediction complete")
print(testSeq)
print(testSecStructSeq)
print(predictedSecStructSeq)
'''

# Uncomment this block to perform prediction with 10 selected data files.
clf = load('models\\clf-td9078-500.joblib')

totalResidue = 0
totalHelix = 0
totalSheet = 0
totalCorrect = 0
totalHelixCorrect = 0
totalSheetCorrect = 0

testDataset = reader.getRS126Dataset("test\\")
for data in testDataset:
	testSet = builder.buildDataSet(data[0], data[1], features, 0)
	result = clf.predict(testSet[0])
	predictedSecStructSeq = builder.buildPredictedSequence(result)
	print("Seq: ", data[0])
	print("Act: ", data[1])
	print("Pre: ", predictedSecStructSeq, "\n")
	
	##### Calcute the prediction stats
	totalResidue += len(data[0])
	for i in range(0, len(data[0])):
		if(data[1][i] == 'H'):
			totalHelix += 1
		elif(data[1][i] == 'E'):
			totalSheet += 1
			
		if(data[1][i] == predictedSecStructSeq[i]):
			totalCorrect += 1
			if(predictedSecStructSeq[i] == 'H'):
				totalHelixCorrect += 1
			elif(predictedSecStructSeq[i] == 'E'):
				totalSheetCorrect += 1

print("Total number of residues:", totalResidue)
print("Total number of Helix residues:", totalHelix)
print("Total number of Sheet residues:", totalSheet)
print("Total number of correct predcition:", totalCorrect)
print("Total number of correct helix prediction:", totalHelixCorrect)
print("Total number of correct sheet prediction:", totalSheetCorrect)
print("Pct total correct prediction:", float(totalCorrect/totalResidue) * 100)
print("Pct correct helix prediction:", float(totalHelixCorrect/totalHelix) * 100)
print("Pct correct sheet prediction:", float(totalSheetCorrect/totalSheet) * 100)
