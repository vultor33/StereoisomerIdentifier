#pragma once
#ifndef STEREOISOMERIDENTIFIER_H
#define STEREOISOMERIDENTIFIER_H

#include <vector>
#include <string>
#include <fstream>
#include <sstream>

#include "AuxMath.h"
#include "Coordstructs.h"
#include "Geometries.h"
#include "MarquesEnantiomers.h"
#include "ReadWriteFormats.h"



class StereoisomerIdentifier
{
public:
	StereoisomerIdentifier();
	~StereoisomerIdentifier();

	void identify(const std::string &fileName_in);

	void generateAllMol(std::string &fileName, int geoCode);

private:
	AuxMath auxMath_;
	double DELTA_CUT;

	// Read input
	Geometries geo_;
	std::ofstream cppOut_;
	std::string molecularFormula;
	std::vector<int> atomTypesCahnIngoldPrelog;
	std::vector< std::vector<int> > chelates;

	//Find shape
	int geoCode;
	double rmsdShape;

	//Find stereoisomer
	int steroisomerIndex;
	std::string stereoLetter;
	std::vector<int> idealTypes;
	std::vector< std::vector<int> > idealChelates;


	std::vector<CoordXYZ> readInput(const std::string &fileName);

	std::vector<CoordXYZ> findShape(const std::vector<CoordXYZ> &coordMol);

	bool calculateEndConditionsWriteIfIsFinished(
		const std::string &fileName,
		const std::vector<CoordXYZ> &coordMol);

	std::string findStereoisomer(
		//input
		std::vector<CoordXYZ> &idealGeo,//add types and chelates
		//CSD
		std::vector<CoordXYZ> &coordMol);//remove metal, add types and chelates
		//output

	void writeOutput(const std::string &fileName, std::string &isomerLine);

	std::string setLabel(int i);

	void reescaleMetalLigandDistancesToOne(std::vector<CoordXYZ> &coord);

	std::vector< std::string > readAllPermutations(
		std::string fileName,
		std::string fileFolder,
		int coordination,
		std::vector<int> &atomTypes,
		std::vector<int> &bidentateChosen);

	std::vector<int> readCauchyNotationsEnantiomers(
		std::ifstream & openendFile_,
		int size);

	std::vector<int> readCauchyNotationsEnantiomers(
		std::string & openendFile_,
		int size);

	void addChelate(std::vector<CoordXYZ> &mol, std::vector< std::vector<int> > &chelates);

	std::vector<std::string> findCountingLine(int coordination, 
		const std::string &shape, 
		const std::string &formula);

	std::string permutationToString0Correction(std::vector<int> &permutation);

	std::vector<int> stringToPermutation(std::string entryString, size_t size);

	std::string filePath(int coordination, const std::string &shape);

	void printMol(const std::vector<CoordXYZ> &mol);

	void printMol(const std::vector<CoordXYZ> &mol, const std::string &fileName);

	void printMol(
		const std::vector<CoordXYZ> &mol, 
		const std::vector< std::vector<int> > &chelates,
		std::string letter,
		const int index);

	void printStereoMol(const std::vector<CoordXYZ> &mol, const std::string &letter, const int index);

	int findIndexByValue(std::vector<int> & vec, int value);
	
	std::string generateStereoisomerFileName(std::string letter, int index);

};

#endif