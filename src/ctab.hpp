/*
 * ctab.h
 *
 *  Created on: Nov 12, 2013
 *      Author: dmitry
 */

#pragma once

#include <string>
#include <vector>
#include <map>
#include <istream>
#include <ostream>
#include "periodic.hpp"

//one atom entry in a MOL file
struct AtomEntry{
	double x,y,z;
	Code code;
	AtomEntry(){}
	AtomEntry(double x_, double y_, double z_, int h, Code code_):
		x(x_), y(y_), z(z_), code(code_){}
};

//one bound entry in a MOL file
struct BoundEntry{
	int a1, a2; // atom indices
	int type;   // bound type
	BoundEntry(){}
	BoundEntry(int a, int b, int type_):
		a1(a), a2(b), type(type_){}
};

struct CTab{
// header block
	std::string name;
	std::string descr;
	std::string comment;
//
	int atomLists, chiral;
// atom block
	std::vector<AtomEntry> atoms;
// bounds block
	std::vector<BoundEntry> bounds;

	CTab() :
		name(), descr(), comment(), atomLists(0), chiral(0), atoms(), bounds(){}
};

//a single entry of SDF database file
struct SDF{
	CTab mol;
	std::map<std::string, std::vector<std::string>> props;
};

CTab readMol(std::istream& inp);
void writeMol(CTab& tab, std::ostream& out);
std::vector<SDF> readSdf(std::istream& inp);
