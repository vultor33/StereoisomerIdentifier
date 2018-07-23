#include "ReadInput.h"

#include <sstream>

using namespace std;

ReadInput::ReadInput(int argc, char *argv[])
{
	if (argc == 3)
	{
		option = "generateAllMol";
		int geoCode;
		stringstream convert0;
		convert0 << argv[1] << "  " << argv[2];
		convert0 >> fileName >> geoCode;
	}
	else if (argc == 2)
	{
		stringstream convert0;
		convert0 << argv[1];
		convert0 >> fileName;
		option = "Identify";
	}
	else
	{
		//fileName = "OC-6-a4bc.csv";
		//option = "generateAllMol";
		//geoCode = 60;
		fileName = "V1-LEMCAH.search1-cpp.inp";
		option = "Identify";
	}




}


