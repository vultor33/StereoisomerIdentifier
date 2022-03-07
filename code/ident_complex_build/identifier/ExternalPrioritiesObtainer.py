from rdkit import Chem
from identifier.ErrorMessages import ErrorMessages

class ExternalPrioritiesObtainer:
	"""Class to calculate CIP rules from a mol file

    Parameters
    -------
		File name in mol format

    Returns
    -------
		CIP priority rules applied to each of its atoms\n
		Raises ErrorMessages().getCipApplicationError if it fails
	"""

	def __init__(self, fileName, molstream = None):
		self.__errorMessages_ = ErrorMessages()
		self.__mol = None
		self.__priorities = []

		self._loadMolFormat(fileName, molstream)
		if self._isPrioritiesSucessfulGenerated():
			self._applyCipRules()

	def getPriorities(self):
		if not self._isPrioritiesSucessfulGenerated():
			raise Exception(self.__errorMessages_.getCipApplicationError())
		return self.__priorities

	def _isPrioritiesSucessfulGenerated(self):
		return not self.__mol is None

	def _loadMolFormat(self,fileName, molstream):
		if molstream != None:
			self.__mol = Chem.MolFromMolBlock(molstream, removeHs = False)	
		else:
			self.__mol = Chem.MolFromMolFile(fileName + '.mol', removeHs = False)

	def _applyCipRules(self):
		Chem.AssignStereochemistry(self.__mol, flagPossibleStereoCenters=True)
		self.__priorities = [int(a.GetProp('_CIPRank')) for a in self.__mol.GetAtoms()]
		self._reverseNumbers()
		
	def _reverseNumbers(self):
		maxPriority = max(self.__priorities)
		for i in range(len(self.__priorities)):
			self.__priorities[i] = maxPriority - self.__priorities[i] 


