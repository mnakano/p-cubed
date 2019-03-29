aminoAcids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
kdValues = [1.8, -4.0, -3.5, -3.5, 2.5, -3.5,  -3.5, 0.4, -3.2, 4.5, 3.8, -3.9, 1.9, 2.8, -1.6, -0.8, -0.7, -0.9, -1.3, 4.2]

class AminoHelper:
	
	def __init__(self):
		self.aminoCount = [0] * 20 # Ni
		self.totalResidue = 0 # Nt
		
		# helix related variables
		self.aminoHelixCount = [0] * 20 # nij
		self.totalHelixResidue = 0 # Nj
		self.freqHelix = 0.0 # fj
		self.helixScoreList = []
		
		# sheet related variables
		self.aminoSheetCount = [0] * 20
		self.totalSheetResidue = 0
		self.freqSheet = 0.0
		self.sheetScoreList = []
	
	def calculateSecStructScore(self, seqDataset): # seqDataset contains a list of sequences and corresponding secondary structure sequence
		
		sequences = []
		ssSequences = []
		
		for data in seqDataset:
			seq = data[0]
			ssSeq = data[1]
			sequences.append(seq)
			ssSequences.append(ssSeq)
			residueList = list(seq)
			ssResidueList = list(ssSeq)
			
			# add the number of residues to the total residues count.
			self.totalResidue += len(residueList)
			# add the number of helix forming resiudes to total helix resiude count.
			self.totalHelixResidue += ssResidueList.count('H')
			# add the number of sheet forming resiudes to total sheet resiude count.
			self.totalSheetResidue += ssResidueList.count('E')
			
			# count the occurence of each amino acid in the evaluated aa sequence.
			for aa in aminoAcids:
				self.aminoCount[aminoAcids.index(aa)] += residueList.count(aa)
			
			# count the occurence of helix forming amino acid in the evaluated aa sequence.
			for i, r in enumerate(ssResidueList):
				if(r == 'H'):
					self.aminoHelixCount[aminoAcids.index(residueList[i])] += 1
				if(r == 'E'):
					self.aminoSheetCount[aminoAcids.index(residueList[i])] += 1
		
		self.freqHelix = self.totalHelixResidue / self.totalResidue # fj
		self.freqSheet = self.totalSheetResidue / self.totalResidue
		
		# calculate helix score for each amino acid and add it to helixScoreList
		for aa in aminoAcids:
			freqHelixAmino = self.aminoHelixCount[aminoAcids.index(aa)] / self.aminoCount[aminoAcids.index(aa)] # fij
			self.helixScoreList.append(freqHelixAmino / self.freqHelix) # Pij = fij / fj
			freqSheetAmino = self.aminoSheetCount[aminoAcids.index(aa)] / self.aminoCount[aminoAcids.index(aa)]
			self.sheetScoreList.append(freqSheetAmino / self.freqSheet)
		
		self.helixScoreList.append(0.0) # Pij for '-'
		self.sheetScoreList.append(0.0)
		
		return([self.helixScoreList, self.sheetScoreList])