#from RunStereoisomerIdentifier import testRunStereoisomerIdentifier
#testRunStereoisomerIdentifier()

#from AllMolecularFormulasGenerator import testAllMolecularFormulasGenerator
#testAllMolecularFormulasGenerator()

#from RunStereoisomerIdentifier import RunStereoisomerIdentifier
#pathInput = "G:\\!CSD-database\\"
#pathOutput = "TestFiles\\"
#extension = ".search1.mol2"
#run_ = RunStereoisomerIdentifier(pathInput,pathOutput, extension)
#run_.runAllFilesFromInputDirectoryWithLimits(306,307)


from RunStereoisomerIdentifier import RunStereoisomerIdentifier
pathInput = "TestFiles\\"
pathOutput = "TestFiles\\"
extension = ".search1.mol2"
calcFilesTemp = [
			"ATIWIK"]
run_ = RunStereoisomerIdentifier(pathInput,pathOutput, extension)
run_.activateKeepIdentifierFiles()
run_.runFromList(calcFilesTemp)
