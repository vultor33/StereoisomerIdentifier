import os.path
import itertools
from collections import Counter

class Enumeration:
    """Class to enumerate stereoisomers"""

    def __init__(self, geoCode):
        self.__allRot = []
        self.__substratalG = []
        self.__substratalR = []
        self.__substratalS = []
        self.__geoFileName = {}
        self.__geoCodeString = {}
        self._defineGeoFileName()
        self.__geoCode = geoCode
        self._readRotations()
        self._readScaffold()

    def setRotations(self, allRotations):
        self.__allRot = allRotations
        
    def setSubstratal(self, substratal, letter):
        self._correctSubstratal(substratal,letter)

    def getEnumerationFileName(self, formula):
        return self.__geoCodeString[self.__geoCode] + '-' + formula + ".csv"

    def _readRotations(self):
        fileName = os.path.join("rotations", "rotations-" + str(self.__geoCode) + ".txt")
        rotInput = open(fileName,"r")
        rotFileStream = rotInput.read().splitlines()
        for line in rotFileStream:
            if line == '':
                break
            rotI = []
            for iNum in line.split():
                rotI.append(int(float(iNum)))
            self.__allRot.append(rotI)
            
    def _readScaffold(self):
        fileName = self.__geoFileName[self.__geoCode]
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
    

    def makeEnumeration(self, formula, colorFirst, chelateFirst, counting):
        enumerationFileName = self.getEnumerationFileName(formula)
        enumOut = open(enumerationFileName, "w")
        referenceLine = self.__geoCodeString[self.__geoCode] + '  '
        referenceLine += formula + '  NC: ' + str(len(colorFirst)) + '   types: '
        for color in colorFirst:
            referenceLine += str(color) + ' '
        referenceLine += '  chelates: ' + str(len(chelateFirst)) + ' '
        for chel in chelateFirst:
            referenceLine += 'cI-length: ' + str(len(chel)) + ' cI: '
            for chelI in chel:
                referenceLine += str(chelI) + ' '

        enumOut.write("{}\n".format(referenceLine))

        # direct counting
        if len(self.__substratalG) != 0:
            i = 0
            stereoRemoved = []
            rcwStereo = []
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
                        rcwStereo.append(i)
                        stereoRemoved.append(j)
                    j+=1
                i+=1

            stereos = list(range(len(self.__substratalG)))
            for remStereo in stereoRemoved:
                stereos.remove(remStereo)
            del stereoRemoved
            rcwStereoCounter = Counter(rcwStereo)
            del rcwStereo
            
            #PRINTING
            enumOut.write("G  {}\n".format(len(stereos)))
            i = 0
            while i < len(self.__substratalG):
                if i in stereos:
                    for isoI in self.__substratalG[i]:
                        enumOut.write("{} ".format(str(isoI)))
                    if i in rcwStereoCounter:
                        enumOut.write(";{}\n".format(
                        str(rcwStereoCounter[i] + 1)))
                    else:
                        enumOut.write(";1\n")
                i+=1
            enumOut.write("R  0\nS  0\n")            
            enumOut.close()
            counting.append(len(stereos))
            return
        

        # R AND S ISOMERS
        chirals = []
        achirals = []
        i = 0
        while i < len(self.__substratalR):
            if self._compareWithRotation(
            self._painting(colorFirst,self.__substratalR[i]),
            self._chelateRepositon(chelateFirst,self.__substratalR[i]),
            self._painting(colorFirst,self.__substratalS[i]),
            self._chelateRepositon(chelateFirst,self.__substratalS[i])):
                achirals.append(i)
            else:
                chirals.append(i)
            i+=1
        
        compAchirals = [achirals]
        compAchirals.append(achirals)
        listCompareAchirals = list(itertools.product(*compAchirals))
        del compAchirals
        achiralsRemoved = []
        rcwAchiral = []
        for iCompare in listCompareAchirals:
            if iCompare[0] >= iCompare[1]:
                continue
            
            if iCompare[0] in achiralsRemoved:
                continue
            if iCompare[1] in achiralsRemoved:
                continue


            if self._compareWithRotation(
            self._painting(colorFirst,self.__substratalR[iCompare[0]]),
            self._chelateRepositon(chelateFirst,self.__substratalR[iCompare[0]]),
            self._painting(colorFirst,self.__substratalR[iCompare[1]]),
            self._chelateRepositon(chelateFirst,self.__substratalR[iCompare[1]])):
                rcwAchiral.append(iCompare[0])
                achiralsRemoved.append(iCompare[1])

        del listCompareAchirals
        for remAchiral in achiralsRemoved:
            achirals.remove(remAchiral)
        del achiralsRemoved
        rcwAchiralCounter = Counter(rcwAchiral)
        del rcwAchiral

        #PRINTING
        enumOut.write("G  {}\n".format(len(achirals)))
        i = 0
        while i < len(self.__substratalR):
            if i in achirals:
                for isoI in self.__substratalR[i]:
                    enumOut.write("{} ".format(str(isoI)))
                if i in rcwAchiralCounter:
                    enumOut.write(";{}\n".format(
                    str(2*rcwAchiralCounter[i] + 2)))
                else:
                    enumOut.write(";2\n")
            i+=1
        counting.append(len(achirals))
        del achirals
        del rcwAchiralCounter

        compChirals = [chirals]
        compChirals.append(chirals)
        listCompareChirals = list(itertools.product(*compChirals))
        del compChirals
        chiralsRemoved = []
        rcwChiral = []
        for iCompare in listCompareChirals:
            if iCompare[0] >= iCompare[1]:
                continue
            
            if iCompare[0] in chiralsRemoved:
                continue
            if iCompare[1] in chiralsRemoved:
                continue
            
            if self._compareWithRotation(
            self._painting(colorFirst,self.__substratalR[iCompare[0]]),
            self._chelateRepositon(chelateFirst,self.__substratalR[iCompare[0]]),
            self._painting(colorFirst,self.__substratalR[iCompare[1]]),
            self._chelateRepositon(chelateFirst,self.__substratalR[iCompare[1]])):
                rcwChiral.append(iCompare[0])
                chiralsRemoved.append(iCompare[1])

            if self._compareWithRotation(
            self._painting(colorFirst,self.__substratalR[iCompare[0]]),
            self._chelateRepositon(chelateFirst,self.__substratalR[iCompare[0]]),
            self._painting(colorFirst,self.__substratalS[iCompare[1]]),
            self._chelateRepositon(chelateFirst,self.__substratalS[iCompare[1]])):
                rcwChiral.append(iCompare[0])
                chiralsRemoved.append(iCompare[1])

        del listCompareChirals
        for remChiral in chiralsRemoved:
            chirals.remove(remChiral)
        del chiralsRemoved
        rcwChiralCounter = Counter(rcwChiral)
        del rcwChiral

        #PRINTING
        enumOut.write("R  {}\n".format(len(chirals)))
        i = 0
        while i < len(self.__substratalR):
            if i in chirals:
                for isoI in self.__substratalR[i]:
                    enumOut.write("{} ".format(str(isoI)))
                if i in rcwChiralCounter:
                    enumOut.write(";{}\n".format(
                    str(rcwChiralCounter[i] + 1)))
                else:
                    enumOut.write(";1\n")
            i+=1
        enumOut.write("S  {}\n".format(len(chirals)))
        i = 0
        while i < len(self.__substratalS):
            if i in chirals:
                for isoI in self.__substratalS[i]:
                    enumOut.write("{} ".format(str(isoI)))
                if i in rcwChiralCounter:
                    enumOut.write(";{}\n".format(
                    str(rcwChiralCounter[i] + 1)))
                else:
                    enumOut.write(";1\n")
            i+=1

        enumOut.close()
        counting.append(len(chirals))


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
        
        lisChelate1.sort()
        lisChelate2.sort()
        
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
        folderScaffoldName = "SubstratalScaffolds"
        self.__geoFileName[40] = os.path.join(folderScaffoldName,'T-4-Mabcd.csv')
        self.__geoFileName[41] = os.path.join(folderScaffoldName,'SP-4-Mabcd.csv')
        self.__geoFileName[42] = os.path.join(folderScaffoldName,'SS-4-Mabcd.csv')
        self.__geoFileName[43] = os.path.join(folderScaffoldName,'vTBPY-4-Mabcd.csv')    
        self.__geoFileName[50] = os.path.join(folderScaffoldName,'TBPY-5-Mabcde.csv')    
        self.__geoFileName[51] = os.path.join(folderScaffoldName,'SPY-5-Mabcde.csv')    
        self.__geoFileName[53] = os.path.join(folderScaffoldName,'PP-5-Mabcde.csv')    
        self.__geoFileName[60] = os.path.join(folderScaffoldName,'OC-6-Mabcdef.csv')    
        self.__geoFileName[61] = os.path.join(folderScaffoldName,'TPR-6-Mabcdef.csv')    
        self.__geoFileName[62] = os.path.join(folderScaffoldName,'HP-6-Mabcdef.csv')    
        self.__geoFileName[63] = os.path.join(folderScaffoldName,'PPY-6-Mabcdef.csv')    
        self.__geoFileName[70] = os.path.join(folderScaffoldName,'COC-7-Mabcdefg.csv')    
        self.__geoFileName[71] = os.path.join(folderScaffoldName,'PBPY-7-Mabcdefg.csv')    
        self.__geoFileName[72] = os.path.join(folderScaffoldName,'CTPR-7-Mabcdefg.csv')    
        self.__geoFileName[73] = os.path.join(folderScaffoldName,'HPY-7-Mabcdefg.csv')    
        self.__geoFileName[74] = os.path.join(folderScaffoldName,'HP-7-Mabcdefg.csv')    
        self.__geoFileName[80] = os.path.join(folderScaffoldName,'SAPR-8-Mabcdefgh.csv')    
        self.__geoFileName[81] = os.path.join(folderScaffoldName,'TDD-8-Mabcdefgh.csv')    
        self.__geoFileName[82] = os.path.join(folderScaffoldName,'BTPR-8-Mabcdefgh.csv')    
        self.__geoFileName[83] = os.path.join(folderScaffoldName,'HBPY-8-Mabcdefgh.csv')    
        self.__geoFileName[84] = os.path.join(folderScaffoldName,'CU-8-Mabcdefgh.csv')    
        self.__geoFileName[85] = os.path.join(folderScaffoldName,'ETBPY-8-Mabcdefgh.csv')    
        self.__geoFileName[86] = os.path.join(folderScaffoldName,'HPY-8-Mabcdefgh.csv')    
        self.__geoFileName[87] = os.path.join(folderScaffoldName,'OP-8-Mabcdefgh.csv')
        
        self.__geoCodeString[40] = 'T-4'    
        self.__geoCodeString[41] = 'SP-4'    
        self.__geoCodeString[42] = 'SS-4'    
        self.__geoCodeString[43] = 'vTBPY-4'    

        self.__geoCodeString[50] = 'TBPY-5'    
        self.__geoCodeString[51] = 'SPY-5'    
        self.__geoCodeString[53] = 'PP-5'    

        self.__geoCodeString[60] = 'OC-6'    
        self.__geoCodeString[61] = 'TPR-6'    
        self.__geoCodeString[62] = 'HP-6'    
        self.__geoCodeString[63] = 'PPY-6'    
        
        self.__geoCodeString[70] = 'COC-7'    
        self.__geoCodeString[71] = 'PBPY-7'    
        self.__geoCodeString[72] = 'CTPR-7'    
        self.__geoCodeString[73] = 'HPY-7'    
        self.__geoCodeString[74] = 'HP-7'    

        self.__geoCodeString[80] = 'SAPR-8'    
        self.__geoCodeString[81] = 'TDD-8'    
        self.__geoCodeString[82] = 'BTPR-8'    
        self.__geoCodeString[83] = 'HBPY-8'    
        self.__geoCodeString[84] = 'CU-8'    
        self.__geoCodeString[85] = 'ETBPY-8'    
        self.__geoCodeString[86] = 'HPY-8'    
        self.__geoCodeString[87] = 'OP-8'
        
    