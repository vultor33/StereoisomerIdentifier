
class Enumeration:
	"""Class to enumerate stereoisomers"""

	def __init__(self, allRotations, substratal):
		self.__allRot = allRotations
		self.__substratal = []
		self._correctSubstratal(substratal)

	def makeEnumeration(self, colorFirst, chelateFirst, rcw, stereoRemoved):
		i = 0
		while i < len(self.__substratal) - 1:
			j = i + 1
			if i in stereoRemoved:
				i += 1
				continue
			while j < len(self.__substratal):
				if j in stereoRemoved:
					j += 1
					continue
				if self._compareWithRotation(
					self._painting(colorFirst,self.__substratal[i]),
					self._chelateRepositon(chelateFirst,self.__substratal[i]),
					self._painting(colorFirst,self.__substratal[j]),
					self._chelateRepositon(chelateFirst,self.__substratal[j])):
					rcw.append(i)
					stereoRemoved.append(j)
					print(i, ' equal ', j)
				j+=1
			i+=1

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
		for chel in lisChelate1:
			chel.sort()
		for chel in lisChelate2:
			chel.sort()
		
		if lisPermut1 == lisPermut2:
			if lisChelate1 == lisChelate2:
				return True

		for rotation in self.__allRot:
			permutRotI = self._applyRotationPermut(lisPermut2,rotation)
			if permutRotI == lisPermut1:
				equal = True
				i = 0
				while i < len(lisChelate2):
					chelRotI = self._applyRotationChelate(lisChelate2[i],rotation)
					chelRotI.sort()
					if lisChelate1[i] != chelRotI: # WARNING - nao sei se a ordem importa - talve substituir por 'chelRotI in lisChelate'
						equal = False
						break
					i+=1
				if equal:
					return equal
	
		return False

	
	def _correctSubstratal(self, substratal):
		for listPermut in substratal:
			newPermutList = []
			for iPer in listPermut:
				newPermutList.append(iPer - 1)
			self.__substratal.append(newPermutList)
	

if __name__ ==  "__main__":
	print("Module to enumerate stereoisomers")
