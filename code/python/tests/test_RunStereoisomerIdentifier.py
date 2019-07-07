from RunStereoisomerIdentifier import RunStereoisomerIdentifier
from Utilities import Utilities
import os.path

def test_RunStereoisomerIdentifier():
    testFilePath = os.path.join("tests", "test_files")
    calculatingFile = os.path.join(testFilePath, "calculating.csv")
    if os.path.isfile(calculatingFile):
        os.remove(calculatingFile)
    pathInput = testFilePath + '\\'
    pathOutput = testFilePath + '\\'

    extension = ".search1.mol2"
    calcFilesTemp = [
            "AVIYAH",
            "TALWOS",
            "BOGNIV01",
            "FAHWAM",
            "ITADAK",
            "LILHIZ",
            "NAGWIB",
            "OTETOY",
            "PIDYOR",
            "QEZCUU",
            "VACGUC01"]
    run_ = RunStereoisomerIdentifier(pathInput,pathOutput, extension)
    run_.runFromList(calcFilesTemp)
    util_ = Utilities()
    assert(util_.isOrderedFilesEqual(run_.getOutputFileName(),run_.getOutputFileName() + '-reference'))
    os.remove(calculatingFile)


