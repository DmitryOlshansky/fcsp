/*
 * fcsp.hpp
 *
 *  Created on: Dec 16, 2013
 *      Author: dmitry
 */
#pragma once

#include <vector>
#include <iostream>
#include <memory>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/graphviz.hpp>
#include "periodic.h"
#include "ctab.h"
#include "descriptors.h"

struct Atom{
	double x, y, z;
	Code code;
	boost::default_color_type color; //used by DFS algorithm
	int path;
	int valence;
	int piE;
	bool inAromaCycle; //part of aromatic cycle
	Atom(){}
	Atom(double x_, double y_, double z_, Code code_) :
		x(x_), y(y_), z(z_), code(code_), path(0), inAromaCycle(false){}
};

struct Bound{
	int type; //
	Bound(){}
	Bound(int type_) :type(type_){}
};

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, Atom, Bound> ChemGraph;

struct Replacement{
	ChemGraph piece;
	int a1, a2; //vertices of replacements
	int dc, coupling;

	Replacement(ChemGraph g, int dc_, int coupling_) :
		piece(std::move(g)), dc(dc_), coupling(coupling_)
	{
		auto asym = Code("A1");
		auto bsym = Code("A2");
		auto r = vertices(piece);
		for (auto i = r.first; i != r.second; i++)
		{
			if (piece[*i].code == asym)
			{
				a1 = *i;
				piece[*i].code = Code("R");
			}
			if (piece[*i].code == bsym)
			{
				a2 = *i;
				piece[*i].code = Code("R");
			}
		}
	}
};

struct FCSP {
	FCSP(std::vector<LevelOne> first, std::vector<LevelTwo> second, std::vector<Replacement> replacements);
	void load(std::istream& inp);
	void dumpGraph(std::ostream& dot);
	void process(std::ostream& out);
	~FCSP();
private:
	struct Impl;
	std::unique_ptr<Impl> pimpl;
};

ChemGraph toGraph(CTab& tab);
void dumpGraph(ChemGraph& graph, std::ostream& out);