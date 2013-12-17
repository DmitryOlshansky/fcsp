#pragma once

#include <vector>
#include <istream>
#include <unordered_map>
#include "periodic.h"

struct LevelOne{
	Code center;
	std::vector<int> valence;
	int index;
	LevelOne(Code center_, std::vector<int> valence_, int index_):
		center(center_), valence(valence_), index(index_){}
};

enum {
	Z = -1,
	R = -2
};

struct Linked{
	int atom;
	int bondType;
};

struct LevelTwo{
	Code center;
	std::vector<Linked> bonds;
	int valence;
	int index;
};

std::unordered_map<int, LevelOne> read1stOrder(std::istream& inp);
std::unordered_map<int, LevelTwo> read2ndOrder(std::istream& inp);
