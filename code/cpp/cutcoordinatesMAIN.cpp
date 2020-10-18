#include <iostream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <fstream>

#include "ReadInput.h"
#include "StereoisomerIdentifier.h"
#include "Geometries.h"
#include "Coordstructs.h"

using namespace std;

int main(int argc, char *argv[])
{
        int geoCode;
        stringstream convert0;
        convert0 << argv[1];
        convert0 >> geoCode;

	Geometries geom_;
	vector<double> rotations;
	vector<CoordXYZ> mol;
	vector<int> reflec;
	double cutangle;
	vector< vector<int> > allreflections;
	ofstream cppOut_;
	cppOut_.open("geometries.json", ofstream::app);

	cppOut_ << "\"" << geom_.sizeToGeometryCode(geoCode) << "\":{" << endl;
	cppOut_ << "  \"coordinates\":[" << endl;
	rotations = geom_.selectGeometry(geoCode,mol,cutangle,reflec);
	for(size_t i = 0; i < mol.size() -1; i++){
		cppOut_ << "    [" <<  mol[i].x << "," << mol[i].y << ", " << mol[i].z << "]," << endl;
	}
	int i = mol.size() - 1;
	cppOut_ << "    [" <<  mol[i].x << "," << mol[i].y << ", " << mol[i].z << "]]," << endl;

	geom_.selectGeometrySymmetries(geoCode, allreflections);
	cppOut_ << "  \"reflections\":{" << endl;
	for(size_t i = 0; i < allreflections.size(); i++){		
		string symmflag = geom_.selectGeometrySymmetriesFlag(geoCode, i, 1);
		cppOut_ << "    \"" << symmflag << "\":[";
		for(size_t j = 0; j < allreflections[i].size() - 1; j++){
			cppOut_ << allreflections[i][j] << ",";
		}
		int j = allreflections[i].size() - 1;
		cppOut_ << allreflections[i][j] << "]";
		if( i == allreflections.size() - 1)
			cppOut_ << "},";
		else
			cppOut_ << ",";
		cppOut_ << endl;
	}

	cppOut_ << "  \"rotations\":{" << endl;
	for(size_t i = 0; i < allreflections.size() - 1; i++){		
		string symmflag = geom_.selectGeometrySymmetriesFlag(geoCode, i, 0);
		cppOut_ << "    \"" << symmflag << "\":[]";
		if( i == allreflections.size() - 2)
			cppOut_ << "}}," << endl;
		else
			cppOut_ << ",";
		cppOut_ << endl;
	}









	return 0;

	ReadInput readInp_ = ReadInput(argc, argv);
	StereoisomerIdentifier stereo;
	string fileName = readInp_.getFileName();
	if (readInp_.getOption() == "Identify")
	{
		stereo.identify(fileName);
	}
	else if (readInp_.getOption() == "generateAllMol")
	{
		stereo.generateAllMol(fileName, readInp_.getGeoCode());
	}
	else
	{
		cout << "Run options not found, check input command-line for errors";
	}
	return 0;
}
