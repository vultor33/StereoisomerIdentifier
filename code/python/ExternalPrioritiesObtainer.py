from rdkit import Chem
from ErrorMessages import ErrorMessages

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

	def __init__(self, fileName):
		self.__errorMessages_ = ErrorMessages()
		self.__mol = None
		self.__priorities = []

		self._loadMolFormat(fileName)
		if self._isPrioritiesSucessfulGenerated():
			self._applyCipRules()

	def getPriorities(self):
		if not self._isPrioritiesSucessfulGenerated():
			raise Exception(self.__errorMessages_.getCipApplicationError())
		return self.__priorities

	def _isPrioritiesSucessfulGenerated(self):
		return not self.__mol is None

	def _loadMolFormat(self,fileName):
		self.__mol = Chem.MolFromMolFile(fileName + '.mol', removeHs = False)
		
	def _loadMol2Format(self,fileName):		
		self.__mol = Chem.MolFromMol2File(fileName + '.mol2', removeHs = False)

	def _applyCipRules(self):
		Chem.AssignStereochemistry(self.__mol, flagPossibleStereoCenters=True)
		self.__priorities = [int(a.GetProp('_CIPRank')) for a in self.__mol.GetAtoms()]
		
	def _applySchneiderRules(self):
		self.__priorities = list(Chem.CanonicalRankAtoms(self.__mol, breakTies=False))
		


