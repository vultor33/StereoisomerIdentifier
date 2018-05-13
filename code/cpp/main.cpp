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

int main()
{
	// ATENCAO  --- preciso centralizar o muffin

	StereoisomerIdentifier stereo;
	stereo.identify("ACPNEU-cpp.inp");
	return 0;
}
