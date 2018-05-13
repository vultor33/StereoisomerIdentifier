#include "StereoisomerIdentifier.h"

#include <sstream>
#include <iostream>

#include "Geometries.h"
#include "MarquesEnantiomers.h"
#include "ReadWriteFormats.h"

StereoisomerIdentifier::StereoisomerIdentifier(){}

StereoisomerIdentifier::~StereoisomerIdentifier(){}

using namespace std;

void StereoisomerIdentifier::identify(const string &fileName)
{
	vector<CoordXYZ> coordMol = readInput(fileName);

	int geoCode;
	double rmsd;
	std::vector<CoordXYZ> idealGeo = findShape(coordMol,geoCode,rmsd);

	vector<CoordXYZ> inputGeometryLigands = coordMol;
	inputGeometryLigands.erase(inputGeometryLigands.begin());
	for (size_t i = 0; i < inputGeometryLigands.size(); i++)
		inputGeometryLigands[i].atomlabel = setLabel(i);

	Geometries geo_;
	cout << "SHAPE:  " << geo_.sizeToGeometryCode(geoCode) << "  -  " << rmsd;

	std::vector<int> atomTypes;
	std::vector<int> bidentateChosen;
	vector<string> allPerm = readAllPermutations(
		geo_.sizeToGeometryCode(geoCode) + ".csv",
		"",
		inputGeometryLigands.size(),
		atomTypes,
		bidentateChosen);


	double minimumDist = 1.0e99;
	int minimumPermut = -1;
	double auxDist;
	MarquesEnantiomers mqRmsd_;
	for (size_t i = 0; i < allPerm.size(); i++)
	{
		vector<CoordXYZ> outGeometry = inputGeometryLigands;
		vector<int> bidentatePermutationRotated = bidentateChosen;
		vector<int> permutationI = stringToPermutation(allPerm[i], atomTypes.size());
		vector<CoordXYZ> molI = idealGeo;
		for (size_t i = 0; i < molI.size(); i++)
			molI[i].atomlabel = setLabel(atomTypes[permutationI[i]]);

		auxDist = mqRmsd_.marquesRmsd(outGeometry, molI);

		cout << "rmsd " << auxDist << endl;

		if (auxDist < minimumDist)
		{
			minimumDist = auxDist;
			minimumPermut = i;
		}
	}

	cout << "rmsd isomero:  " << minimumDist << endl
		<< "isomero:  " << minimumPermut << endl;




}

string StereoisomerIdentifier::setLabel(int i)
{
	switch (i)
	{
	case 0:
		return "H";
	case 1:
		return "He";
	case 2:
		return "Li";
	case 3:
		return "Be";
	case 4:
		return "B";
	case 5:
		return "C";
	case 6:
		return "N";
	case 7:
		return "O";
	case 8:
		return "F";
	case 9:
		return "Na";
	case 10:
		return "Mg";
	case 11:
		return "Al";
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
	vector<double> xInp, yInp, zInp;
	vector<int> rank;
	while (getline(input_, line))
	{
		if (line == "end")
			break;

		double auxX, auxY, auxZ;
		int auxRank;
		stringstream convert;
		convert << line;
		convert >> auxX >> auxY >> auxZ >> auxRank;
		xInp.push_back(auxX);
		yInp.push_back(auxY);
		zInp.push_back(auxZ);
		rank.push_back(auxRank);
	}
	vector<int> instructions = auxMath_.vector_ordering(rank);
	auxMath_.vector_ordering_with_instructions(xInp, instructions);
	auxMath_.vector_ordering_with_instructions(yInp, instructions);
	auxMath_.vector_ordering_with_instructions(zInp, instructions);
	size_t nCoord = xInp.size();
	vector<CoordXYZ> coordMol;
	CoordXYZ auxMol;
	auxMol.x = xInp[nCoord - 1];
	auxMol.y = yInp[nCoord - 1];
	auxMol.z = zInp[nCoord - 1];
	coordMol.push_back(auxMol);
	for (size_t i = 0; i < rank.size() - 1; i++)
	{
		CoordXYZ auxMol2;
		auxMol2.x = xInp[i];
		auxMol2.y = yInp[i];
		auxMol2.z = zInp[i];
		coordMol.push_back(auxMol2);
	}
	for (size_t i = 0; i < coordMol.size(); i++)
		coordMol[i].atomlabel = "H";
	reescaleMetalLigandDistancesToOne(coordMol);
	return coordMol;
}

std::vector<CoordXYZ> StereoisomerIdentifier::findShape(
	const std::vector<CoordXYZ> &coordMol,
	int &geoCode,
	double &rmsdShape)
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
		for (size_t i = 0; i < coord.size(); i++)
			coord[i].atomlabel = "H";
		double rmsd = mrq_.marquesRmsd(coord, coord2);
		if (rmsd < rmsdMin)
		{
			iMin = i;
			rmsdMin = rmsd;
		}

		ofstream out_;
		out_.open(".teste.xyz", std::ofstream::out | std::ofstream::app);
		out_ << coordMol.size() << endl << endl;
		for (size_t i = 0; i < coordMol.size(); i++)
		{
			out_ << coord[i].atomlabel << "  "
				<< coord[i].x << "  "
				<< coord[i].y << "  "
				<< coord[i].z << endl;
		}
		out_.close();
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