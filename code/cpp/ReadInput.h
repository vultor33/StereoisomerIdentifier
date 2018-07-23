#ifndef READINPUT_H
#define READINPUT_H

#include <string>

class ReadInput
{
public:
	ReadInput(int argc, char *argv[]);

	std::string getOption() { return option; }

	std::string getFileName() { return fileName; }

	int getGeoCode() { return geoCode; }


private:
	std::string option;
	std::string fileName;
	int geoCode;
};


#endif
