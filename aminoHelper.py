aminoAcids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

class AminoHelper:
	
	def __init__(self):
		self.aminoCount = [0] * 20 # Ni
		self.totalResidue = 0 # Nt
		
		# helix related variables
		self.aminoHelixCount = [0] * 20 # nij
		self.totalHelixResidue = 0 # Nj
		self.freqHelix = 0.0 # fj
		self.helixScoreList = []
	
	def calculateHelixScore(self, seqDataset): # seqDataset contains a list of sequences and corresponding secondary structure sequence
		
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
			
			# count the occurence of each amino acid in the evaluated aa sequence.
			for aa in aminoAcids:
				self.aminoCount[aminoAcids.index(aa)] += residueList.count(aa)
			
			# count the occurence of helix forming amino acid in the evaluated aa sequence.
			for i, h in enumerate(ssResidueList):
				if(h == 'H'):
					self.aminoHelixCount[aminoAcids.index(residueList[i])] += 1
		
		self.freqHelix = self.totalHelixResidue / self.totalResidue # fj
		
		# calculate helix score for each amino acid and add it to helixScoreList
		for aa in aminoAcids:
			freqHelixAmino = self.aminoHelixCount[aminoAcids.index(aa)] / self.aminoCount[aminoAcids.index(aa)] # fij
			self.helixScoreList.append(freqHelixAmino / self.freqHelix) # Pij = fij / fj
		
		self.helixScoreList.append(0.0) # Pij for '-'
		
		return(self.helixScoreList)