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

	std::vector<CoordXYZ> readInput(const std::string &fileName);

	std::vector<CoordXYZ> findShape(
		const std::vector<CoordXYZ> &coordMol,
		int &geoCode,
		double &rmsdShape);

	void printMol(const std::vector<CoordXYZ> &mol);

	std::vector< std::string > readAllPermutations(
		std::string fileName,
		std::string fileFolder,
		int coordination,
		std::vector<int> &atomTypes,
		std::vector<int> &bidentateChosen);

	std::vector<int> readCauchyNotationsEnantiomers(
		std::ifstream & openendFile_,
		int size);

	std::string permutationToString0Correction(std::vector<int> &permutation);

	std::vector<int> stringToPermutation(std::string entryString, size_t size);

};

#endif