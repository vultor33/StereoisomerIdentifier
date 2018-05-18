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

using namespace std;

int main(int argc, char *argv[])
{
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
	StereoisomerIdentifier stereo;
	stereo.identify(fileName);
	return 0;
}
