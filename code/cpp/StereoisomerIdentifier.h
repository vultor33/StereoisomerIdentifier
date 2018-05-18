#pragma once
#ifndef STEREOISOMERIDENTIFIER_H
#define STEREOISOMERIDENTIFIER_H

#include <vector>
#include <string>
#include <fstream>

#include "AuxMath.h"
#include "Coordstructs.h"



class StereoisomerIdentifier
{
public:
	StereoisomerIdentifier();
	~StereoisomerIdentifier();

	void identify(const std::string &fileName);

private:
	std::string setLabel(int i);

	void reescaleMetalLigandDistancesToOne(std::vector<CoordXYZ> &coord);

	AuxMath auxMath_;

	std::vector<CoordXYZ> readInput(
		const std::string &fileName, 
		std::string &formula,
		std::vector<int> &atomTypesCahnIngoldPrelog);

	std::vector<CoordXYZ> findShape(
		const std::vector<CoordXYZ> &coordMol,
		int &geoCode,
		double &rmsdShape);

	std::string findStereoisomer(
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

	std::vector<std::string> findCountingLine(int coordination, 
		const std::string &shape, 
		const std::string &formula);

	std::string permutationToString0Correction(std::vector<int> &permutation);

	std::vector<int> stringToPermutation(std::string entryString, size_t size);

	std::string filePath(int coordination, const std::string &shape);

	void printMol(const std::vector<CoordXYZ> &mol);

	void printMol(const std::vector<CoordXYZ> &mol, const std::string &fileName);


};

#endif