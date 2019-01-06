#include "StereoisomerIdentifier.h"

#include <sstream>
#include <iostream>
#include <algorithm>

#include "Geometries.h"
#include "MarquesEnantiomers.h"
#include "ReadWriteFormats.h"

using namespace std;

StereoisomerIdentifier::StereoisomerIdentifier(){}

StereoisomerIdentifier::~StereoisomerIdentifier(){}

void StereoisomerIdentifier::identify(const string &fileName)
{
	DELTA_CUT = 0.15e0;

	cppOut_.open((fileName + ".log").c_str());
	
	vector<CoordXYZ> coordMol = readInput(fileName);

	std::vector<CoordXYZ> idealGeo = findShape(coordMol);

	bool isStereoisomerAlreadyIdentified = calculateEndConditionsWriteIfIsFinished(
		fileName,
		coordMol);
	if (isStereoisomerAlreadyIdentified)
		return;

	string isomerLine = findStereoisomer(
		idealGeo, //add types and chelates
		coordMol//remove metal, add types and chelates
	);

	writeOutput(fileName, isomerLine);

	cppOut_.close();
}


void StereoisomerIdentifier::writeOutput(const std::string &fileName, std::string &isomerLine)
{
	if (isomerLine == "")
	{
		cppOut_ << fileName << endl
			<< "rmsdfailed" << endl
			<< -1 << endl;
	}
	else
	{
		cppOut_ << fileName << endl
			<< geo_.sizeToGeometryCode(geoCode) << "-" << stereoLetter << "-" << steroisomerIndex << endl
			<< rmsdShape << endl;
	}
}

bool StereoisomerIdentifier::calculateEndConditionsWriteIfIsFinished(
	const std::string &fileName, 
	const std::vector<CoordXYZ> &coordMol)
{
	if (geoCode == 20)
	{
		double angleCoord = 180.0e0 - auxMath_.angleFrom3Points(
			coordMol[1].x, coordMol[1].y, coordMol[1].z,
			coordMol[0].x, coordMol[0].y, coordMol[0].z,
			coordMol[2].x, coordMol[2].y, coordMol[2].z) * (180 / auxMath_._pi);

		if (angleCoord > 5.0e0)
		{
			cppOut_ << fileName << endl
				<< "A-2" << endl
				<< rmsdShape << endl;
		}
		else
		{
			cppOut_ << fileName << endl
				<< "L-2" << endl
				<< rmsdShape << endl;
		}
		cppOut_.close();
		return true;
	}
	else if (rmsdShape > DELTA_CUT)
	{
		cppOut_ << fileName << endl
			<< "failed" << endl
			<< rmsdShape << endl;
		cppOut_.close();
		return true;
	}
	else if (geoCode == 32)
	{
		cppOut_ << fileName << endl
			<< geo_.sizeToGeometryCode(geoCode) << endl
			<< rmsdShape << endl;
		cppOut_.close();
		return true;
	}
	return false;
}


void StereoisomerIdentifier::generateAllMol(string &fileName, int geoCode)
{
	ReadWriteFormats rwf_;
	ifstream file_(fileName.c_str());
	vector<int> atomTypes;
	vector< vector<int> > chel;
	rwf_.readAtomTypesAndChelates(file_, atomTypes, chel);

	Geometries geo_;
	vector<CoordXYZ> idealGeo;
	double dummy;
	vector<int> reflecDummy;
	geo_.selectGeometry(
		geoCode,
		idealGeo,
		dummy,
		reflecDummy);

	string line;
	string letterType;
	int stereoIndexLine = 0;
	while (getline(file_, line))
	{
		if ((line[0] == 'G') || (line[0] == 'R') || (line[0] == 'S'))
		{
			letterType = line[0];
			stereoIndexLine = 0;
			continue;
		}
		stereoIndexLine++;

		stringstream convertLine;
		convertLine << line;
		vector<int> permutation(atomTypes.size());
		for (size_t i = 0; i < atomTypes.size(); i++)
			convertLine >> permutation[i];

		for (size_t k = 0; k < idealGeo.size(); k++)
			idealGeo[k].atomlabel = setLabel(atomTypes[permutation[k]]);
		vector< vector<int> > chelPermuted = chel;
		for (size_t i = 0; i < chelPermuted.size(); i++)
		{
			for (size_t j = 0; j < chelPermuted[i].size(); j++)
				chelPermuted[i][j] = findIndexByValue(permutation, chelPermuted[i][j]);
		}

		printMol(
			idealGeo,
			chelPermuted,
			letterType,
			stereoIndexLine);
	}
}

string StereoisomerIdentifier::setLabel(int i)
{
	switch (i)
	{
	case 0:
		return "O";
	case 1:
		return "P";
	case 2:
		return "N";
	case 3:
		return "I";
	case 4:
		return "Na";
	case 5:
		return "Ca";
	case 6:
		return "Se";
	case 7:
		return "C";
	case 8:
		return "Ti";
	case 9:
		return "He";
	case 10:
		return "B";
	case 11:
		return "Sc";
	default:
		cout << "label not found - exiting" << endl;
		exit(1);
	}
}


void StereoisomerIdentifier::reescaleMetalLigandDistancesToOne(vector<CoordXYZ> &coord)
{
	// put metal at 0
	for (size_t i = 1; i < coord.size(); i++)
	{
		coord[i].x -= coord[0].x;
		coord[i].y -= coord[0].y;
		coord[i].z -= coord[0].z;
		double rInv = 1.0e0 / auxMath_.norm(coord[i].x, coord[i].y, coord[i].z);
		coord[i].x *= rInv;
		coord[i].y *= rInv;
		coord[i].z *= rInv;
	}
	coord[0].x = 0.0e0;
	coord[0].y = 0.0e0;
	coord[0].z = 0.0e0;
}


std::vector<CoordXYZ> StereoisomerIdentifier::readInput(const std::string &fileName)
{
	ifstream input_(fileName.c_str());
	string line;
	getline(input_, line);
	stringstream convert1;
	convert1 << line;
	convert1 >> molecularFormula;
	getline(input_, line);
	if (line != "")
	{
		stringstream convertChelates;
		convertChelates << line;
		ReadWriteFormats rwf_;
		chelates = rwf_.readChelates(convertChelates);
	}

	string dum1;
	vector<double> xInp, yInp, zInp;
	while (getline(input_, line))
	{
		if (line == "end")
			break;

		double auxX, auxY, auxZ;
		int auxRank;
		stringstream convert;
		convert << line;
		convert >> dum1 >> auxX >> auxY >> auxZ >> auxRank;
		xInp.push_back(auxX);
		yInp.push_back(auxY);
		zInp.push_back(auxZ);
		atomTypesCahnIngoldPrelog.push_back(auxRank);
	}

	size_t nCoord = xInp.size();
	vector<CoordXYZ> coordMol;
	for (size_t i = 0; i < atomTypesCahnIngoldPrelog.size(); i++)
	{
		CoordXYZ auxMol2;
		auxMol2.atomlabel = "H";
		auxMol2.x = xInp[i];
		auxMol2.y = yInp[i];
		auxMol2.z = zInp[i];
		coordMol.push_back(auxMol2);
	}
	reescaleMetalLigandDistancesToOne(coordMol);
	return coordMol;
}

std::vector<CoordXYZ> StereoisomerIdentifier::findShape(const std::vector<CoordXYZ> &coordMol)
{
	Geometries geo_;
	vector<int> avaibleGeometries = geo_.avaibleGeometries(coordMol.size() - 1);
	CoordXYZ metal;
	metal.atomlabel = "H";
	metal.x = 0.0e0;
	metal.y = 0.0e0;
	metal.z = 0.0e0;

	int iMin = -1;
	double rmsdMin = 1.0e99;
	MarquesEnantiomers mrq_;
	for (size_t i = 0; i < avaibleGeometries.size(); i++)
	{
		vector<CoordXYZ> coord;
		vector<CoordXYZ> coord2 = coordMol;
		double dummy;
		vector<int> reflecDummy;
		geo_.selectGeometry(
			avaibleGeometries[i],
			coord,
			dummy,
			reflecDummy);
		coord.insert(coord.begin(), metal);

		//printStereoMol(coord2, "EXPERIMENTAL", 1);//activate for debug
		//printStereoMol(coord, "X", avaibleGeometries[i]); //activate for debug

		for (size_t i = 0; i < coord.size(); i++)
			coord[i].atomlabel = "H";
		double rmsd = mrq_.marquesRmsdEqualMass(coord, coord2); //rmsd is different from previous value because the metal was removed
		
		if (rmsd < rmsdMin)
		{
			iMin = i;
			rmsdMin = rmsd;
		}
	}
	geoCode = avaibleGeometries[iMin];
	rmsdShape = rmsdMin;

	vector<CoordXYZ> coord;
	double dummy;
	vector<int> reflecDummy;
	geo_.selectGeometry(
		avaibleGeometries[iMin],
		coord,
		dummy,
		reflecDummy);
	return coord;

}


std::string StereoisomerIdentifier::findStereoisomer(
	//input
	std::vector<CoordXYZ> &idealGeo,//add types and chelates
	//CSD
	std::vector<CoordXYZ> &coordMol) //remove metal, add types and chelates
{
	Geometries geo_;
	ReadWriteFormats rwf_;
	MarquesEnantiomers mqRmsd_;

	coordMol.erase(coordMol.begin());
	atomTypesCahnIngoldPrelog.erase(atomTypesCahnIngoldPrelog.begin());//metal removed
	for (size_t i = 0; i < coordMol.size(); i++)
		coordMol[i].atomlabel = setLabel(atomTypesCahnIngoldPrelog[i]);

	string fileIsomerPath = filePath(coordMol.size(), geo_.sizeToGeometryCode(geoCode));
	fileIsomerPath += geo_.sizeToGeometryCode(geoCode) + "-" + molecularFormula + ".csv";
	ifstream fileIsomers_(fileIsomerPath.c_str());
	rwf_.readAtomTypesAndChelates(
		fileIsomers_,
		idealTypes,
		idealChelates);

	addChelate(coordMol, chelates); //add dummy atoms to coordMol

	int stereoIndexLine = 0;
	double minimumRmsd = 1.0e99;
	double auxRmsd;
	string line;
	string letterType = "";
	string minimumLine;
	while (getline(fileIsomers_, line))
	{
		if ((line[0] == 'G') || (line[0] == 'R') || (line[0] == 'S'))
		{
			letterType = line[0];
			stereoIndexLine = 0;
			continue;
		}
		stereoIndexLine++;

		vector<CoordXYZ> outGeometry = coordMol;
		vector<CoordXYZ> molI = idealGeo;

		stringstream convertLine;
		convertLine << line;
		vector<int> permutation(idealTypes.size());
		for (size_t i = 0; i < idealTypes.size(); i++)
			convertLine >> permutation[i];

		for (size_t k = 0; k < molI.size(); k++)
			molI[k].atomlabel = setLabel(idealTypes[permutation[k]]);
		vector< vector<int> > chelPermuted = idealChelates;
		for (size_t i = 0; i < chelPermuted.size(); i++)
		{
			for (size_t j = 0; j < chelPermuted[i].size(); j++)
				chelPermuted[i][j] = findIndexByValue(permutation, chelPermuted[i][j]);
		}
		addChelate(molI, chelPermuted);

		//printStereoMol(molI, letterType, stereoIndexLine); activate for debug

		auxRmsd = mqRmsd_.marquesRmsdEqualMass(outGeometry, molI);
		if (auxRmsd == -1.0e0)
		{
			continue;
			//minimumRmsd = -1.0e0;
			//break;
		}

		if (auxRmsd < minimumRmsd)
		{
			minimumRmsd = auxRmsd;
			minimumLine = line;
			steroisomerIndex = stereoIndexLine;
			stereoLetter = letterType;
		}
	}
	fileIsomers_.close();


	/* USE FOR TESTS 
	vector<int> permutationI(idealGeo.size());
	stringstream convertLine;
	convertLine << minimumLine;
	for (size_t i = 0; i < idealTypes.size(); i++)
		convertLine >> permutationI[i];
	for (size_t k = 0; k < idealGeo.size(); k++)
		idealGeo[k].atomlabel = setLabel(idealTypes[permutationI[k]]);
	for (size_t i = 0; i < idealChelates.size(); i++)
	{
		for (size_t j = 0; j < idealChelates[i].size(); j++)
			idealChelates[i][j] = findIndexByValue(permutationI, idealChelates[i][j]);
	}
	addChelate(idealGeo, idealChelates);
	string idealFile = fileName + "-ideal.xyz";
	string csdFile = fileName + "-csd.xyz";
	printMol(
		idealGeo,
		idealChelates,
		stereoLetter,
		idealFile);
	printMol(
		coordMol,
		chelates,
		stereoLetter,
		csdFile);
	*/


	if (minimumRmsd == -1.0e0)
		return "";
	else
		return minimumLine;


}


void StereoisomerIdentifier::printMol(const std::vector<CoordXYZ> &mol, const std::string &fileName)
{
	ofstream out_(fileName.c_str());
	out_ << mol.size() << endl << endl;
	for (size_t i = 0; i < mol.size(); i++)
	{
		out_ << mol[i].atomlabel << "  "
			<< mol[i].x * 3.0e0 << "  "
			<< mol[i].y * 3.0e0 << "  "
			<< mol[i].z * 3.0e0 << endl;
	}
	out_.close();

}

void StereoisomerIdentifier::printMol(const std::vector<CoordXYZ> &mol)
{
	ofstream out_("test-mol.xyz");
	out_ << mol.size() << endl << endl;
	for (size_t i = 0; i < mol.size(); i++)
	{
		out_ << mol[i].atomlabel << "  "
			<< mol[i].x << "  "
			<< mol[i].y << "  "
			<< mol[i].z << endl;
	}
	out_.close();

}

void StereoisomerIdentifier::printMol(
	const std::vector<CoordXYZ> &mol,
	const std::vector< std::vector<int> > &chelates,
	std::string letter,
	const int index)
{
	string fileName = generateStereoisomerFileName(letter, index);

	//WATCH ATOM TYPES OF THE CHELATES POINTS
	ofstream out_;
	//out_.open(fileName.c_str(), std::ofstream::out | std::ofstream::app);
	out_.open(fileName.c_str());
	double reescale = 2.0e0;
	int countChel = 0;
	for (size_t i = 0; i < chelates.size(); i++)
		for (size_t j = 0; j < chelates[i].size() - 1; j++)
			for (size_t k = j + 1; k < chelates[i].size(); k++)
				countChel += 2;


	out_ << mol.size() + countChel << endl << letter << endl;
	for (size_t i = 0; i < mol.size(); i++)
	{
		out_ << mol[i].atomlabel << "  "
			<< mol[i].x * reescale << "  "
			<< mol[i].y * reescale << "  "
			<< mol[i].z * reescale << endl;
	}
	for (size_t i = 0; i < chelates.size(); i++)
	{
		for (size_t j = 0; j < chelates[i].size() - 1; j++)
		{
			for (size_t k = j + 1; k < chelates[i].size(); k++)
			{
				size_t chel1 = chelates[i][j];
				size_t chel2 = chelates[i][k];
				CoordXYZ meanPoint;
				meanPoint.x = reescale * (mol[chel2].x - mol[chel1].x);
				meanPoint.y = reescale * (mol[chel2].y - mol[chel1].y);
				meanPoint.z = reescale * (mol[chel2].z - mol[chel1].z);
				out_ << "H" << "  "
					<< reescale * mol[chel1].x + meanPoint.x * 0.25e0 << "  "
					<< reescale * mol[chel1].y + meanPoint.y * 0.25e0 << "  "
					<< reescale * mol[chel1].z + meanPoint.z * 0.25e0 << endl;
				out_ << "H" << "  "
					<< reescale * mol[chel1].x + meanPoint.x * 0.75e0 << "  "
					<< reescale * mol[chel1].y + meanPoint.y * 0.75e0 << "  "
					<< reescale * mol[chel1].z + meanPoint.z * 0.75e0 << endl;
			}
		}
	}

	out_.close();
}

string StereoisomerIdentifier::generateStereoisomerFileName(std::string letter, int index)
{
	stringstream convert;
	convert << index;
	return letter + "-" + convert.str() + ".xyz";
}

void StereoisomerIdentifier::printStereoMol(const std::vector<CoordXYZ> &mol, const std::string &letter, const int index)
{
	string fileName = generateStereoisomerFileName(letter, index);
	printMol(mol, fileName);
}

void StereoisomerIdentifier::addChelate(std::vector<CoordXYZ> &mol, std::vector< std::vector<int> > &chelates)
{
	for (size_t i = 0; i < chelates.size(); i++)
	{
		for (size_t j = 0; j < chelates[i].size() - 1; j++)
		{
			for (size_t k = j + 1; k < chelates[i].size(); k++)
			{
				size_t chel1 = chelates[i][j];
				size_t chel2 = chelates[i][k];
				CoordXYZ auxPoint1;
				auxPoint1.atomlabel = "H";
				auxPoint1.x = (mol[chel2].x - mol[chel1].x);
				auxPoint1.y = (mol[chel2].y - mol[chel1].y);
				auxPoint1.z = (mol[chel2].z - mol[chel1].z);
				CoordXYZ auxPoint2;
				auxPoint2.atomlabel = "H";
				auxPoint2.x = mol[chel1].x + auxPoint1.x * 0.25e0;
				auxPoint2.y = mol[chel1].y + auxPoint1.y * 0.25e0;
				auxPoint2.z = mol[chel1].z + auxPoint1.z * 0.25e0;
				auxPoint1.x = mol[chel1].x + auxPoint1.x * 0.75e0;
				auxPoint1.y = mol[chel1].y + auxPoint1.y * 0.75e0;
				auxPoint1.z = mol[chel1].z + auxPoint1.z * 0.75e0;
				mol.push_back(auxPoint1);
				mol.push_back(auxPoint2);
			}
		}
	}
}

std::vector< std::string > StereoisomerIdentifier::readAllPermutations(
	std::string fileName,
	std::string fileFolder,
	int coordination,
	std::vector<int> &atomTypes,
	std::vector<int> &bidentateChosen)
{
	// read input
	ReadWriteFormats rwf_;
	ifstream fileIsomers_((fileFolder + fileName).c_str());
	int nBidentates = 0;
	rwf_.readAtomTypesAndBidentateChosenFileWithLabels(
		fileIsomers_,
		atomTypes,
		bidentateChosen,
		coordination,
		nBidentates);
	string line;
	vector<string> allPermut;
	while (!fileIsomers_.eof())
	{
		vector<int> permutation = readCauchyNotationsEnantiomers(fileIsomers_, coordination);
		if (permutation.size() == 0)
			continue;
		string permtString = permutationToString0Correction(permutation);
		allPermut.push_back(permtString);
	}
	fileIsomers_.close();
	return allPermut;
}

vector<int> StereoisomerIdentifier::readCauchyNotationsEnantiomers(ifstream & openendFile_, int size)
{
	vector<int> notation;
	if (openendFile_.eof())
		return notation;

	string auxline;
	getline(openendFile_, auxline);
	if (auxline == "")
		return notation;
	notation.resize(size);

	size_t brack1Temp = auxline.find("]");
	size_t brack1 = auxline.find("[", brack1Temp + 1, 1);
	size_t brack2 = auxline.find("]", brack1Temp + 1, 1);
	string permString = auxline.substr(brack1 + 1, brack2 - brack1 - 1);
	stringstream line;
	line << permString;
	for (int i = 0; i < size; i++)
	{
		line >> notation[i];
	}
	return notation;
}

vector<int> StereoisomerIdentifier::readCauchyNotationsEnantiomers(string & line, int size)
{
	vector<int> notation;

	if (line == "")
		return notation;
	notation.resize(size);

	size_t brack1Temp = line.find("]");
	size_t brack1 = line.find("[", brack1Temp + 1, 1);
	size_t brack2 = line.find("]", brack1Temp + 1, 1);
	string permString = line.substr(brack1 + 1, brack2 - brack1 - 1);
	stringstream convert;
	convert << permString;
	for (int i = 0; i < size; i++)
	{
		convert >> notation[i];
	}
	return notation;
}

string StereoisomerIdentifier::permutationToString0Correction(vector<int> &permutation)
{
	stringstream permt;
	for (size_t i = 0; i < permutation.size(); i++)
		permt << (permutation[i] - 1) << " ";
	return permt.str();
}

vector<int> StereoisomerIdentifier::stringToPermutation(string entryString, size_t size)
{
	stringstream convert;
	convert << entryString;
	vector<int> permutation(size);
	for (size_t i = 0; i < size; i++)
	{
		convert >> permutation[i];
	}
	return permutation;
}

std::string StereoisomerIdentifier::filePath(int coordination, const std::string &shape)
{
	string file = "Stereoisomerlist/";
	stringstream convert;
	string coordString;
	convert << coordination;
	convert >> coordString;
	file += "CN" + coordString + "/";
	file += shape + "/";
	return file;
}

std::vector<std::string> StereoisomerIdentifier::findCountingLine(int coordination,
	const std::string &shape,
	const std::string &formula)
{
	string file = filePath(coordination, shape);
	file += shape + "-counting.csv";
	ifstream count_(file.c_str());
	string line;
	vector<string> numbersLine;
	while (getline(count_, line))
	{
		string auxFormula = line.substr(0,line.find(";"));
		if (auxFormula == formula)
		{
			numbersLine.push_back(line);
			string line2;
			while (getline(count_, line2))
			{
				if (line2[0] != ';')
					break;
				numbersLine.push_back(line2);
			}
			break;
		}
	}
	return numbersLine;
}

int StereoisomerIdentifier::findIndexByValue(std::vector<int> & vec, int value)
{
	vector<int>::iterator it = find(vec.begin(), vec.end(), value);

	if (it != vec.end())
		return (int)distance(vec.begin(), it);
	else
	{
		cout << "Error on:  StereoisomerIdentifier::findIndexByValue -- exiting" << endl;
		exit(1);
	}
}
