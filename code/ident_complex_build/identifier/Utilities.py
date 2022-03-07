import os
from identifier.ErrorMessages import ErrorMessages

class Utilities:
	def __init__(self):
		self.__ErrorMessages_ = ErrorMessages()

		#parameters
		self.__LINE_NOT_FOUND = 'This line wasnt found on file'
		self.__allMetals = ["Sc","Ti","Cr","Mn","Fe","Co","Ni","Cu","Zn",
		"Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
		"La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
		"Hf","Ta","Re","Os","Ir","Pt","Au","Hg",
		"Ac","Th","Pa","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr",
		"Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Uub",
		"Al","Ga","Ge","In","Sn","Sb","Tl","Pb","Bi","Po",
		"V","Y","W","U"]

	def compareTwoSetOfLinesInOrder(self,lines1, lines2):
		for i in range(len(lines1)):
			if lines1[i] != lines2[i]:
				print('line: ', str(i), ' is different')
				print('1: ',lines1[i])
				print('2: ',lines2[i])
				return False

		return True		

	def compareTwoSetOfLinesOutOfOrder(self,lines1, lines2):
		for line1I in lines1:
			if not line1I in lines2:
					print(self.__LINE_NOT_FOUND + ' 2: ',line1I)
					return False
		for line2I in lines2:
			if not line2I in lines1:
					print(self.__LINE_NOT_FOUND + ' 1: ',line2I)
					return False

		return True		

	def isOrderedFilesEqual(self,fileName1, fileName2):
		if not os.path.isfile(fileName1):
			print(self.__ErrorMessages_.getIOError() + fileName1)
			return False
		if not os.path.isfile(fileName2):
			print(self.__ErrorMessages_.getIOError() + fileName2)
			return False
		
		file1 = open(fileName1, "r")
		file2 = open(fileName2, "r")
		fileStream1 = file1.read().splitlines()
		fileStream2 = file2.read().splitlines()

		areEqual = self.compareTwoSetOfLinesInOrder(fileStream1,fileStream2)
		return areEqual


	def isOutOfOrderFilesEqual(self,fileName1, fileName2):
		if not os.path.isfile(fileName1):
			print(self.__ErrorMessages_.getIOError() + fileName1)
			return False
		if not os.path.isfile(fileName2):
			print(self.__ErrorMessages_.getIOError() + fileName2)
			return False

		file1 = open(fileName1, "r")
		file2 = open(fileName2, "r")
		fileStream1 = file1.read().splitlines()
		fileStream2 = file2.read().splitlines()

		areEqual = self.compareTwoSetOfLinesOutOfOrder(fileStream1,fileStream2)
		return areEqual

	def dictFindKeyByValue(self, inputDictionary, search_value):
		for key, value in inputDictionary.items():
			if value == search_value:
				return key

	def justLetters(self, enterString):
		return ''.join([i for i in enterString if i.isalpha()])

	def untilPoint(self, enterString):
		return enterString.partition(".")[0]
	
	def isMetal(self, metalString):
		return self.untilPoint(metalString) in self.__allMetals
		
		
		
		
		
