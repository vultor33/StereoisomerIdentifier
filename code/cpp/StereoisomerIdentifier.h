#pragma once
#ifndef STEREOISOMERIDENTIFIER_H
#define STEREOISOMERIDENTIFIER_H

#include <vector>
#include <string>
#include <fstream>
#include <sstream>

#include "AuxMath.h"
#include "Coordstructs.h"



class StereoisomerIdentifier
{
public:
	StereoisomerIdentifier();
	~StereoisomerIdentifier();

	void identify(const std::string &fileName_in);

	void identifyMonoVersion(const std::string &fileName);

	void generateAllMol(std::string &fileName, int geoCode);

private:
	std::string fileName;

	std::string setLabel(int i);

	void reescaleMetalLigandDistancesToOne(std::vector<CoordXYZ> &coord);

	AuxMath auxMath_;

	std::vector<CoordXYZ> readInput(
		const std::string &fileName, 
		std::string &formula,
		std::vector<int> &atomTypesCahnIngoldPrelog,
		std::vector< std::vector<int> > & chelates);

	std::vector<CoordXYZ> readInputMonoVersion(
		const std::string &fileName,
		std::string &formula,
		std::vector<int> &atomTypesCahnIngoldPrelog);


	std::vector<CoordXYZ> findShape(
		const std::vector<CoordXYZ> &coordMol,
		int &geoCode,
		double &rmsdShape);

	std::string findStereoisomer(
		//input
		std::vector<CoordXYZ> &idealGeo,//add types and chelates
		const std::string &molecularFormula,
		const int geoCode,
		//CSD
		std::vector<CoordXYZ> &coordMol,//remove metal, add types and chelates
		std::vector<int> &atomTypesCahnIngoldPrelog,
		std::vector< std::vector<int> > &chelates,
		//output
		std::vector<int> &idealTypes,
		std::vector< std::vector<int> > &idealChelates,
		int &steroisomerIndex,
		std::string &stereoLetter);

	std::string findStereoisomerMonoVersion(
		const std::string &molecularFormula,
		const int geoCode,
		int &minimumLine,
		std::string &stereoLetter,
		const std::vector<CoordXYZ> &idealGeo,
		std::vector<CoordXYZ> &coordMol,
		std::vector<int> &atomTypes,
		std::vector<int> &atomTypesCahnIngoldPrelog);


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
		std::string title,
		const std::string &fileName);

	int findIndexByValue(std::vector<int> & vec, int value);
};

#endif