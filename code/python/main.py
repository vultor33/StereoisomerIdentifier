import sys
from RunStereoisomerIdentifier import RunStereoisomerIdentifier

typedCommands = sys.argv

if len(typedCommands) != 2:
    print('Input format:')
    print('python main.py [fileName]')
else:
    pathInput = ""
    pathOutput = ""
    extension = ""
    calcFilesTemp = [str(typedCommands[1])]
    run_ = RunStereoisomerIdentifier(pathInput,pathOutput, extension)
    run_.activateKeepIdentifierFiles()
    run_.runFromList(calcFilesTemp)






