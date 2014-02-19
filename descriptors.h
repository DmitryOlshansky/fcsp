#pragma once

#include <vector>
#include <istream>
#include <unordered_map>
#include "periodic.h"

struct LevelOne{
	Code center;
	int valence;
	int dc;
	LevelOne(Code center_, int valence_, int dc_):
		center(center_), valence(valence_), dc(dc_){}
	bool operator<(const LevelOne& rhs) const
	{
		return center < rhs.center || (center == rhs.center && valence < rhs.valence);
	}
};

struct Linked{
	Code atom;
	int bondType;
	Linked(){}
	Linked(Code atom_, int bondType_):
		atom(atom_), bondType(bondType_){}
};

struct LevelTwo{
	Code center;
	int valence;
	std::vector<Linked> bonds;
	int dc;
	LevelTwo(Code center_, int valence_, std::vector<Linked> linked, int dc_) :
		center(center_), valence(valence_), bonds(linked), dc(dc_){}
	bool operator<(const LevelOne& rhs) const
	{
		return center < rhs.center || (center == rhs.center && valence < rhs.valence);
	}
};

//read
void read1stOrder(std::istream& inp, std::vector<LevelOne> &dest);
void read2ndOrder(std::istream& inp, std::vector<LevelTwo> &dest);
