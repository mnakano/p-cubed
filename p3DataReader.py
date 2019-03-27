import os # Python's directory navigation module
import re # Python's regex module

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
	return rs126Dataset
	
def getTD9078Dataset(filename):
	df = pandas.read_csv(filename)
	return(df[df['has_nonstd_aa'] == False][['seq', 'sst3']])