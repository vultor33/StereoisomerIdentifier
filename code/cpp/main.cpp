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
	/* ATENCAO  --- preciso centralizar o muffin
	Sobre o "Get Your Atoms in Order An Open-Source Implementation of a Novel and Robust Molecular Canonicalization Algorithm"
	a mudanca no indice muda o reordamento canonico, se o atomo nao sair do lugar, a mudanca do reordamento muda o isomero.

	*/
	string fileName;
	if (argc > 1)
	{
		stringstream convert0;
		convert0 << argv[1];
		convert0 >> fileName;
	}
	else
	{
		fileName = "ACPNEU-cpp.inp";
		//cout << "Input file not found - exiting";
		//return 1;
	}
	StereoisomerIdentifier stereo;
	stereo.identify(fileName);
	return 0;
}
