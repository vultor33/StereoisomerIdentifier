import glob


def splitIsomerFile(fileName): 
	if fileName.find('counting') > -1:
		return
	isomerFile = open(fileName,"r")
	fileStream = isomerFile.read().splitlines()
	isomerFile.close()
	listR = []
	listS = []
	listG = []
	i = 1
	while i < len(fileStream):
		if fileStream == '':
			continue
		if i < len(fileStream)-1:
			if fileStream[i+1] == '':
				listG.append(fileStream[i])
				i+=2
			else:
				listR.append(fileStream[i])
				listS.append(fileStream[i+1])
				i+=3
		else:
			listG.append(fileStream[i])
			i+=1

	isomerFile = open(fileName, "w")
	isomerFile.write("{}\n".format(fileStream[0]))
	isomerFile.write("G\n")
	for isomer in listG:
		isomerFile.write("{}\n".format(isomer))
	isomerFile.write("R\n")
	for isomer in listR:
		isomerFile.write("{}\n".format(isomer))
	isomerFile.write("S\n")
	for isomer in listS:
		isomerFile.write("{}\n".format(isomer))
	isomerFile.close()
	return


stereoList = glob.glob('StereoisomerList//*')
for cn in stereoList:
	allGeo = glob.glob(cn + '//*')
	for geo in allGeo:
		allIsoFiles = glob.glob(geo + '//*')
		for isoFile in allIsoFiles:
			splitIsomerFile(isoFile)
			


if __name__ ==  "__main__":
	print("This module transform old format of stereoisomer into RSG format")