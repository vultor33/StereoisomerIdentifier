#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>

#include "Coordstructs.h"
#include "Geometries.h"
#include "MarquesEnantiomers.h"
#include "AuxMath.h"
#include "StereoisomerIdentifier.h"
#include "ReadWriteFormats.h"

using namespace std;

int main(int argc, char *argv[])
{
	StereoisomerIdentifier stereo;
	stereo.generateAllMol();
	return 0;




	string fileName;
	if (argc > 1)
	{
		stringstream convert0;
		convert0 << argv[1];
		convert0 >> fileName;
	}
	else
	{
		//fileName = "ACPNEU-cpp.inp";
		cout << "Input file not found - exiting";
		return 1;
	}
	stereo.identify(fileName);
	return 0;
}
