from Utilities import Utilities
from AllMolecularFormulasGenerator import AllMolecularFormulasGenerator
import os.path

def test_AllMolecularFormulasGenerator():
    useAllMolecularFormulasGenerator(6)
    util_ = Utilities()
    allFormulas_ = AllMolecularFormulasGenerator()
    outputFileName = allFormulas_.getOutputFileName()
    testFilePath = os.path.join("TestFiles", outputFileName + '-reference')
    assert(util_.isOutOfOrderFilesEqual(outputFileName,testFilePath))


def useAllMolecularFormulasGenerator(nCoordinationMax):
    allFormulas_ = AllMolecularFormulasGenerator()
    fileName = allFormulas_.getOutputFileName()
    allFormFile = open(fileName, "w")
    allFormFile.write("allFormList = [[0]]\n")
    for i in range(1,nCoordinationMax + 1):
        objG = AllMolecularFormulasGenerator()
        objG.generateAllFormulas(i)
        objG.printAllFormulasToStream(allFormFile)
    allFormFile.close()
