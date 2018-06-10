import csv
from DataAllFormulas import allFormList

class StatisticalAnalysis:
	"""Class to analyze identification results"""


	def __init__(self, fileName):
		self._fileName = fileName
		self.fileCodes = open('CSD-codes-with-errors.txt','w')
		self.__metalStatistics = {}
		


	def analyze(self):
		k = 0
		l = 0
		erdk = 0
		enl = 0
		egrap = 0
		efor = 0
		e0met = 0
		enol = 0
		epol = 0
		nmetalsmax = 0
		metMaxCode = ''
		mixl = 0
		ntotalMix = 0


		kTl = 0	
		with open(self._fileName, 'r') as csvfile:
			spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')
			for auxRow in spamreader:
				row = auxRow[0:auxRow.index('')]
				if row[0] == 'CSD':
					continue

				cont = False
				for elem in row:
					if "E." in elem:
						cont = True
						if "E.rdkit" in elem:
							erdk += 1
						elif "E.NL>8" in elem:
							enl += 1
						elif "E.Graph" in elem:
							egrap += 1
						elif "E.Formula" in elem:
							efor += 1
						elif "E.0Metals" in elem:
							e0met += 1
						elif "E.NoLigands" in elem:
							enol += 1
						elif "E.Polyedron" in elem:
							epol += 1
						break
				if cont:
					continue
					
				self._analyzeAllMetalsRow(row)
			
				if row[1] == 'Mix':
					mixl += 1
					nmet = (len(row)-2)/4
					ntotalMix += nmet
					if nmet > nmetalsmax:
						nmetalsmax = nmet
						metMaxCode = row[0]
				
			
				if len(row) > 3:
					if row[1] == 'Dimetal':
						stereoID1 = row[4].split('-')
						stereoID2 = row[8].split('-')
						if len(stereoID1) > 4 and len(stereoID2) > 4:
							k +=1
							if stereoID1[0] == stereoID2[0]:
								l+=1
								#print('1:  ',stereoID1[4], '  2:  ',stereoID2[4])
								#print('1:  ',stereoID1[3], '  2:  ',stereoID2[3])
							else:
								continue
								#print(row)
		
	
	
		print('mix:  ',mixl)
		print('numero de metais avulsos:  ',ntotalMix)
		print("DIMETALS")
		print('total:  ',k)
		print('same code: ',l)
		
		csvMetalStatistics = open('csvMetalStatistics.csv','w')		
		for key in self.__metalStatistics:
			csvMetalStatistics.write(key + ';')
			for value in self.__metalStatistics[key]:
				csvMetalStatistics.write(str(value) + ';')
			csvMetalStatistics.write('\n')
		csvMetalStatistics.close()

		#print(nmetalsmax)
		#print(metMaxCode)
		#ERRORS
		print("rdk:  ",erdk)
		print("enl:  ",enl)
		print("grap:  ",egrap)
		print("efor:  ",efor)
		print("emet:  ",e0met)
		print("enol:  ",enol)
		print("epol:   ",epol)
		print("soma:  ",erdk + enl + egrap + efor + e0met + enol + epol)



	def _analyzeAllMetalsRow(self,row):
		#dar sort no metal depois
		nmet = (len(row)-2)/4
		i = 0

		while i < nmet:
			metal = self._findMetal(row[4 * i + 2])
			chelN = self._chelateNumber(row[4 * i + 3])
			stereoID = row[4 * i + 4].split('-')
			letter = ''
			nCoord = ''
			geo = ''
			if stereoID[1] == 'TP' or stereoID[2] == '2' or stereoID[2] == '1':
				letter = 'G'
				nCoord = stereoID[2]
			else:
				nCoord = stereoID[2]
				letter = stereoID[3]
			
			#montar uma estatistica de cada metal.
			# 
		
			self._generateExcelStatisticsForMetals(
				metal,
				nCoord,
				letter,self.__allGeo.index(stereoID[1] + '-' + stereoID[2]),
				self.__allPossibleChel.index(chelN))
		
			i+=1
	

	def _generateExcelStatisticsForMetals(self,metal,nCoord,letterI,geoI,chelI):
		stat = []
		if metal in self.__metalStatistics:
			stat = self.__metalStatistics[metal]
		else:
			stat = [0] * 63
		
		stat[int(float(nCoord)) - 1] += 1
		if letterI == 'G':
			stat[8] +=1
		elif letterI == 'R':
			stat[9] += 1
		elif letterI == 'S':
			stat[10] += 1
		stat[11+geoI] += 1
		stat[41+chelI] += 1
		
		self.__metalStatistics[metal] = stat
		
	


	def _findMetal(self,metal):
		for m in self.__allMetals:
			if m in metal:
				return m

	def _chelateNumber(self,formula):
		for lineFor in allFormList:
			for iFormula in lineFor:
				if iFormula == 0:
					continue
				if iFormula[1] == formula:
					chel = iFormula[4]
					nChelates = []
					for iChel in chel:
						nChelates.append(len(iChel))
					nChelates.sort()
					return nChelates

	def _generateAllPossibleChel(self):
		allChel = []
		for lineFor in allFormList:
			for iFormula in lineFor:
				if iFormula == 0:
					continue
				chel = iFormula[4]
				nChelates = []
				for iChel in chel:
					nChelates.append(len(iChel))
				nChelates.sort()
				if not nChelates in allChel:
					allChel.append(nChelates)
		print(allChel)
	
	def _rowToCsvString(self,row):
		strOut = ''
		for elem in row:
			strOut += elem + ';'
		return strOut
	
	
	

	#one letter metals must come last
	__allMetals = ["Sc","Ti","Cr","Mn","Fe","Co","Ni","Cu","Zn",
		"Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
		"La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
		"Hf","Ta","Re","Os","Ir","Pt","Au","Hg",
		"Ac","Th","Pa","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr",
		"Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Uub",
		"Al","Ga","Ge","In","Sn","Sb","Tl","Pb","Bi","Po",
		"V","Y","W","U"]
	__allPossibleChel = [[], [2], [3], [4], [2, 2], [5], [2, 3], [6], [3, 3], [2, 2, 2], [2, 4], [7], [2, 5], [3, 4], [2, 2, 3], [8], [4, 4], [2, 2, 2, 2], [2, 6], [2, 3, 3], [3, 5], [2, 2, 4]]
	__allGeo = ["L-1",
		"A-2","L-2",
		"TP-3","TPY-3","TS-3",
		"SP-4","SS-4","T-4","vTBPY-4",
		"PP-5","SPY-5","TBPY-5",
		"HP-6","OC-6","PPY-6","TPR-6",
		"COC-7","CTPR-7","HP-7","HPY-7","PBPY-7",
		"BTPR-8","CU-8","ETBPY-8","HBPY-8","HPY-8","OP-8","SAPR-8","TDD-8"]
	__allLetters = ["G","R","S"]

if __name__ ==  "__main__":
	print("Module to analize stereoisomer identification")
