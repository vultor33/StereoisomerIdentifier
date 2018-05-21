
class Enumeration:
	"""Class to enumerate stereoisomers"""

	def __init__(self):
		self.__allRot = []
		self.__substratalG = []
		self.__substratalR = []
		self.__substratalS = []

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
		fileName = 'Stereoisomerlist//CN4//SP-4//SP-4-Mabcd.csv'  # LEEEEEEEEEEEEEEEEEEEEERRRR A PARTIR DO GEO CODE
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
				#print(i, ' - ', j)
				if self._compareWithRotation(#inside compare:
				self._painting(colorFirst,self.__substratalG[i]),
				self._chelateRepositon(chelateFirst,self.__substratalG[i]),
				self._painting(colorFirst,self.__substratalG[j]),
				self._chelateRepositon(chelateFirst,self.__substratalG[j])):
					#print('removed: ', j)
					rcw.append(i)
					stereoRemoved.append(j)
				j+=1
			i+=1
		

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
	

if __name__ ==  "__main__":
	print("Module to enumerate stereoisomers")
