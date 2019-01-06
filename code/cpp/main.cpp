#include <iostream>
#include <string>
#include <sstream>
#include <stdlib.h>

#include "ReadInput.h"
#include "StereoisomerIdentifier.h"

using namespace std;

int main(int argc, char *argv[])
{
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
