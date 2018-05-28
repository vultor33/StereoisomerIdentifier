#include <iostream>
#include <string>
#include <sstream>
#include <stdlib.h>

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
		fileName = "ABELUO.search1-cpp.inp";
		//cout << "Input file not found - exiting";
		//return 1;
	}
	StereoisomerIdentifier stereo;
	stereo.identify(fileName);
	return 0;
}
