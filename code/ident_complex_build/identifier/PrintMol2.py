REESCALE = 3.0
ATOMLABELS = ['Au','C','Se','Na','O','Ca','I','P','N','He','B']  
FLAGS = [
    '#   generated by: ComplexBuild\n' +
    '\n\n' +
    '@<TRIPOS>MOLECULE\n' +
    'name\n',
    'SMALL\n' +
    'NO_CHARGES\n' +
    '\n\n' +
    '@<TRIPOS>ATOM\n',
    '@<TRIPOS>BOND\n',
    'M  END\n']

def printMol2(fileName, colors, chelates, refGeometry):
    printMol = defineMOL2(colors, chelates, refGeometry)
    file = open(fileName, 'w')
    file.write(printMol)
    file.close()

def defineMOL2(colors, chelates, refGeometry):
    natoms, printAtoms = definePrintAtomsMOL2(colors, refGeometry)
    nbonds, printBonds = definePrintAllBondsMOL2(chelates, refGeometry)
    printMol2 = ''
    printMol2 += FLAGS[0]
    printMol2 += str(natoms) + ' ' + str(nbonds) + '\n'
    printMol2 += FLAGS[1]
    printMol2 += printAtoms
    printMol2 += FLAGS[2]
    printMol2 += printBonds
    printMol2 += FLAGS[3]
    return printMol2

#####################################################################
#####################################################################

def definePrintAtomsMOL2(atomColors, refGeometry):
    printAllAtoms = ''
    atomXyz = [0,0,0]
    atomColor = -1
    i = 1
    printAllAtoms += atomlineMOL2(i, atomColor, atomXyz)
    for i in range(2, len(refGeometry) + 2):
        printAllAtoms += atomlineMOL2(i, atomColors[i-2], refGeometry[i-2])
    natoms = i
    return natoms, printAllAtoms
    
def atomlineMOL2(i, atomColor, atomXyz):
    istr = str(i)
    aflag = 'A' + istr
    printatom = ' ' + istr + ' ' + aflag + '   '
    printatom += str(atomXyz[0] * REESCALE) + '   '
    printatom += str(atomXyz[1] * REESCALE) + '   '
    printatom += str(atomXyz[2] * REESCALE) + '  '
    printatom += ATOMLABELS[atomColor + 1] + '\n'
    return printatom


#####################################################################
#####################################################################

def definePrintAllBondsMOL2(chelatePermutation, refGeometry):
    printAllBonds = ''
    for i in range(1, len(refGeometry) + 1):
        printAllBonds += metalbondsMOL2(i)
    
    i2 = len(refGeometry) + 1
    for chelate in chelatePermutation:
        i2, printline = isomerChelateMOL2(i2, chelate)
        printAllBonds += printline

    nbonds = i2 - 1
    return nbonds, printAllBonds

def metalbondsMOL2(i):
    printline = ' ' + str(i) + ' 1  ' + str(i+1) + ' 1\n'
    return printline

def isomerChelateMOL2(i, chelate):
    i2 = i
    printline = ''
    for iso1 in range(len(chelate) - 1):
        for iso2 in range(iso1 + 1,len(chelate)):
            chelI = chelate[iso1]
            chelJ = chelate[iso2]
            printline += isomerBidantateMOL2(i2, chelate)
            i2 += 1
    return i2, printline

def isomerBidantateMOL2(i, chelate):
    printline = ' ' + str(i) + ' '
    printline += str(chelate[0] + 2) + '  ' + str(chelate[1] + 2)
    printline += ' 1\n'
    return printline