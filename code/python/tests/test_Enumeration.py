import os
import os.path
from Utilities import Utilities
from Enumeration import Enumeration
from DataAllFormulas import allFormList

# TODO -- alem do identificador, um arquivo counting e gerado. O counting tambem precisa ser testado aqui.

def test_Enumeration():
    geoCode = 60
    ncoord = 6
    i = 6
    formula = allFormList[ncoord][i][2]
    colorFirst = allFormList[ncoord][i][3]
    chelation = allFormList[ncoord][i][4]
    enumeration = Enumeration(geoCode)
    fileName = enumeration.getEnumerationFileName(formula)
    if os.path.isfile(fileName):
        os.remove(fileName)

    counting = []
    enumeration.makeEnumeration(formula, colorFirst, chelation, counting)  
    util_ = Utilities()
    referenceFile = os.path.join("tests", "test_files")
    referenceFile = os.path.join(referenceFile, fileName)
    assert(util_.isOutOfOrderFilesEqual(fileName,referenceFile))
    os.remove(fileName)
    
