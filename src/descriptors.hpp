#pragma once

#include <vector>
#include <istream>
#include <unordered_map>
#include "periodic.hpp"
#include "chemgraph.hpp"

struct LevelOne{
    Code center;
    int valence;
    int dc;
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
    int start;
    int valence;
    bool monolith;
    bool replOnly; // if only allowed in replacements
    std::vector<Linked> bonds;
    int dc;
    bool operator<(const LevelTwo& rhs) const
    {
        return center < rhs.center || (center == rhs.center && valence < rhs.valence);
    }
};


struct Replacement{
    ChemGraph piece;
    int a1, a2; //vertices of replacements
    int dc, coupling;

    Replacement(ChemGraph g, int dc_, int coupling_);
};


//read
std::vector<LevelOne> read1stOrder(std::istream& inp);
std::vector<LevelTwo> read2ndOrder(std::istream& inp);
std::vector<Replacement> readReplacements(std::istream& inp);
