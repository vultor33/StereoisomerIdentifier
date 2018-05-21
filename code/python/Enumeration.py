import itertools



class Enumeration:
	"""Class to enumerate stereoisomers"""

	def __init__(self):
		self.__allRot = []
		self.__substratalG = []
		self.__substratalR = []
		self.__substratalS = []
		self.__geoFileName = {}
		self._defineGeoFileName()
		

	def setRotations(self, allRotations):
		self.__allRot = allRotations
		
	def setSubstratal(self, substratal, letter):
		self._correctSubstratal(substratal,letter)

	def readRotations(self, geoCode):
		fileName = "rotations//rotations-" + str(geoCode) + ".txt"
		rotInput = open(fileName,"r")
		rotFileStream = rotInput.read().splitlines()
		for line in rotFileStream:
			rotI = []
			for iNum in line.split():
				rotI.append(int(float(iNum)))
			self.__allRot.append(rotI)
			
	def readScaffold(self, geoCode):
		fileName = self.__geoFileName[geoCode]
		scaInput = open(fileName,"r")
		scaFileStream = scaInput.read().splitlines()
		i = 2
		iLetter = 'G'
		substratalG = []
		substratalR = []
		substratalS = []
		while i < len(scaFileStream):
			if scaFileStream[i] == '':
				break
			if scaFileStream[i] == 'R':
				iLetter = 'R'
				i+=1
				continue
			if scaFileStream[i] == 'S':
				iLetter = 'S'
				i+=1
				continue
			if iLetter == 'G':
				substratalG.append(self._takePermutation(scaFileStream[i]))
			elif iLetter == 'R':
				substratalR.append(self._takePermutation(scaFileStream[i]))
			elif iLetter == 'S':
				substratalS.append(self._takePermutation(scaFileStream[i]))
			i+=1
		
		self._correctSubstratal(substratalG, 'G')
		self._correctSubstratal(substratalR, 'R')
		self._correctSubstratal(substratalS, 'S')
	

	def makeEnumeration(self, colorFirst, chelateFirst, rcw, stereoRemoved): # REFAZEEER
		i = 0
		while i < len(self.__substratalG) - 1:
			j = i + 1
			if i in stereoRemoved:
				i += 1
				continue
			while j < len(self.__substratalG):
				if j in stereoRemoved:
					j += 1
					continue
				if self._compareWithRotation(#inside compare:
				self._painting(colorFirst,self.__substratalG[i]),
				self._chelateRepositon(chelateFirst,self.__substratalG[i]),
				self._painting(colorFirst,self.__substratalG[j]),
				self._chelateRepositon(chelateFirst,self.__substratalG[j])):
					rcw.append(i)
					stereoRemoved.append(j)
				j+=1
			i+=1


		chirals = []
		achirals = []
		
		i = 0
		while i < len(self.__substratalR):
			if self._compareWithRotation(
			self._painting(colorFirst,self.__substratalR[i]),
			self._chelateRepositon(chelateFirst,self.__substratalR[i]),
			self._painting(colorFirst,self.__substratalS[i]),
			self._chelateRepositon(chelateFirst,self.__substratalS[i])):
				rcw.append(i)
				achirals.append(i)
				stereoRemoved.append(i)
			else:
				chirals.append(i)
			i+=1
		
		compAchirals = [achirals]
		compAchirals.append(achirals)
		compChirals = [chirals]
		compChirals.append(chirals)

		listCompareAchirals = list(itertools.product(*compAchirals))
		listCompareChirals = list(itertools.product(*compChirals))
		del compAchirals
		del compChirals
		
		achiralsRemoved = []
		for iCompare in listCompareAchirals:
			if iCompare[0] >= iCompare[1]:
				continue
			
			if iCompare[1] in achiralsRemoved:
				continue

			if self._compareWithRotation(
			self._painting(colorFirst,self.__substratalR[iCompare[0]]),
			self._chelateRepositon(chelateFirst,self.__substratalR[iCompare[0]]),
			self._painting(colorFirst,self.__substratalR[iCompare[1]]),
			self._chelateRepositon(chelateFirst,self.__substratalR[iCompare[1]])):
				rcw.append(iCompare[0])
				achiralsRemoved.append(iCompare[1])

		chiralsRemoved = []
		for iCompare in listCompareChirals:
			if iCompare[0] >= iCompare[1]:
				continue
			
			if iCompare[1] in chiralsRemoved:
				continue
			
			if self._compareWithRotation(
			self._painting(colorFirst,self.__substratalR[iCompare[0]]),
			self._chelateRepositon(chelateFirst,self.__substratalR[iCompare[0]]),
			self._painting(colorFirst,self.__substratalR[iCompare[1]]),
			self._chelateRepositon(chelateFirst,self.__substratalR[iCompare[1]])):
				rcw.append(iCompare[0])
				rcw.append(iCompare[0])
				chiralsRemoved.append(iCompare[1])
				chiralsRemoved.append(iCompare[1])

			if self._compareWithRotation(
			self._painting(colorFirst,self.__substratalR[iCompare[0]]),
			self._chelateRepositon(chelateFirst,self.__substratalR[iCompare[0]]),
			self._painting(colorFirst,self.__substratalS[iCompare[1]]),
			self._chelateRepositon(chelateFirst,self.__substratalS[iCompare[1]])):
				rcw.append(iCompare[0])
				rcw.append(iCompare[0])
				chiralsRemoved.append(iCompare[1])
				chiralsRemoved.append(iCompare[1])

		#no stereo estao os S que foram removidos
		#no achiral estao os R que foram removidos e que sao achirais
		#a lista dos chirais removem o R e o S
		print(stereoRemoved)
		print(achiralsRemoved)
		print(chiralsRemoved)

		stereoRemoved += achiralsRemoved + chiralsRemoved



		

	def printInfo(self):
		print('sca G: ',self.__substratalG)
		print('sca R: ',self.__substratalR)
		print('sca S: ',self.__substratalS)


	#Rotating permutations
	def _applyRotationPermut(self, lisPermut, rotation):
		rotatedPermut = list(lisPermut)
		i = 0
		while i < len(lisPermut):
			rotatedPermut[rotation[i]] = lisPermut[i]
			i+=1
		return rotatedPermut
	def _applyRotationChelate(self, chel, rotation):
		chelRotated = [];
		for k in chel:
			chelRotated.append(rotation[k])
		return chelRotated
		
	#Painting
	def _painting(self, colors, lisPermut):
		coloredPermut = []
		for iPer in lisPermut:
			coloredPermut.append(colors[iPer])
		return coloredPermut
	def _chelateRepositon(self, initialChel,lisPermut):
		lisChelate = []
		for chel in initialChel:
			auxChel = []
			for iChel in chel:
				auxChel.append(lisPermut.index(iChel))
			lisChelate.append(auxChel)
		return lisChelate

	#COMPARING
	def _compareWithRotation(self, lisPermut1, lisChelate1, lisPermut2, lisChelate2):
		for iChel in range(len(lisChelate1)):
			lisChelate1[iChel].sort()
			lisChelate2[iChel].sort()
		
		if lisPermut1 == lisPermut2:
			if lisChelate1 == lisChelate2:
				return True

		for rotation in self.__allRot:
			permutRotI = self._applyRotationPermut(lisPermut2,rotation)
			if permutRotI == lisPermut1:
				equal = True
				i = 0
				chelRotI = []
				for i in range(len(lisChelate2)):
					rotI = self._applyRotationChelate(lisChelate2[i],rotation)
					rotI.sort()
					chelRotI.append(rotI)
				for i in range(len(chelRotI)):
					if not chelRotI[i] in lisChelate1:
						equal = False
						break
				if equal:
					return equal
	
		return False

	
	def _takePermutation(self, line):
		i = 0
		firstBrac = False
		permut = []
		while i < len(line):
			if line[i] == '[' and firstBrac:
				i+=1
				while line[i] != ']':
					if line[i] != ' ':
						permut.append(int(float(line[i])))
					i+=1
			elif line[i] == '[':
				firstBrac = True
			i+=1
		return permut
			
		
	
	def _correctSubstratal(self, substratal, letter):
		for listPermut in substratal:
			newPermutList = []
			for iPer in listPermut:
				newPermutList.append(iPer - 1)
			if letter == 'G':
				self.__substratalG.append(newPermutList)
			if letter == 'R':
				self.__substratalR.append(newPermutList)
			if letter == 'S':
				self.__substratalS.append(newPermutList)
	
	def _defineGeoFileName(self):
	
		self.__geoFileName[40] = 'Stereoisomerlist//CN4//T-4//T-4-Mabcd.csv'	
		self.__geoFileName[41] = 'Stereoisomerlist//CN4//SP-4//SP-4-Mabcd.csv'	
		self.__geoFileName[42] = 'Stereoisomerlist//CN4//SS-4//SS-4-Mabcd.csv'	
		self.__geoFileName[43] = 'Stereoisomerlist//CN4//vTBPY-4//vTBPY-4-Mabcd.csv'	

		self.__geoFileName[50] = 'Stereoisomerlist//CN5//TBPY-5//TBPY-5-Mabcde.csv'	
		self.__geoFileName[51] = 'Stereoisomerlist//CN5//SPY-5//SPY-5-Mabcde.csv'	
		self.__geoFileName[53] = 'Stereoisomerlist//CN5//PP-5//PP-5-Mabcde.csv'	

		self.__geoFileName[60] = 'Stereoisomerlist//CN6//OC-6//OC-6-Mabcdef.csv'	
		self.__geoFileName[61] = 'Stereoisomerlist//CN6//TPR-6//TPR-6-Mabcdef.csv'	
		self.__geoFileName[62] = 'Stereoisomerlist//CN6//HP-6//HP-6-Mabcdef.csv'	
		self.__geoFileName[63] = 'Stereoisomerlist//CN6//PPY-6//PPY-6-Mabcdef.csv'	
		
		self.__geoFileName[70] = 'Stereoisomerlist//CN7//COC-7//COC-7-Mabcdefg.csv'	
		self.__geoFileName[71] = 'Stereoisomerlist//CN7//PBPY-7//PBPY-7-Mabcdefg.csv'	
		self.__geoFileName[72] = 'Stereoisomerlist//CN7//CTPR-7//CTPR-7-Mabcdefg.csv'	
		self.__geoFileName[73] = 'Stereoisomerlist//CN7//HPY-7//HPY-7-Mabcdefg.csv'	
		self.__geoFileName[74] = 'Stereoisomerlist//CN7//HP-7//HP-7-Mabcdefg.csv'	

		self.__geoFileName[80] = 'Stereoisomerlist//CN8//SAPR-8//SAPR-8-Mabcdefgh.csv'	
		self.__geoFileName[81] = 'Stereoisomerlist//CN8//TDD-8//TDD-8-Mabcdefgh.csv'	
		self.__geoFileName[82] = 'Stereoisomerlist//CN8//BTPR-8//BTPR-8-Mabcdefgh.csv'	
		self.__geoFileName[83] = 'Stereoisomerlist//CN8//HBPY-8//HBPY-8-Mabcdefgh.csv'	
		self.__geoFileName[84] = 'Stereoisomerlist//CN8//CU-8//CU-8-Mabcdefgh.csv'	
		self.__geoFileName[85] = 'Stereoisomerlist//CN8//ETBPY-8//ETBPY-8-Mabcdefgh.csv'	
		self.__geoFileName[86] = 'Stereoisomerlist//CN8//HPY-8//HPY-8-Mabcdefgh.csv'	
		self.__geoFileName[87] = 'Stereoisomerlist//CN8//OP-8//OP-8-Mabcdefgh.csv'	
	
	


if __name__ ==  "__main__":
	print("Module to enumerate stereoisomers")
