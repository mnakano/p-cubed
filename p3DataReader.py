'''
p3DataReader.py
The library contains functions related to reading amino acid sequence data into the program.
'''

import os # Python's directory navigation module
import re # Python's regex module
import pandas

def readFile(filename):
	with open(filename) as f:
		lines = f.readlines()
	return(''.join(lines))

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

def getRS126Dataset(dir):
	filenames = os.listdir(dir)
	rs126Dataset = []
	for filename in filenames:
		content = readFile(dir + filename)
		rs126Dataset.append(extractSequenceSet(content))
	return(rs126Dataset)
	
def getTD9078Dataset(filename, numSamples):
	td9078Dataset = []
	df = pandas.read_csv(filename)
	df = df[df['has_nonstd_aa'] == False][['seq', 'sst3']].head(numSamples)
	for i in range(0, len(df.index)):
		list = []
		list.append(df.iloc[i,0])
		list.append(df.iloc[i,1])
		td9078Dataset.append(list)
	return(td9078Dataset)