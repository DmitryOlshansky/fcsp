/*
 * fcsp.hpp
 *
 *  Created on: Dec 16, 2013
 *      Author: dmitry
 */
#pragma once

#include <vector>
#include <string>
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

struct AtomVertex{
	Code code;
// temporaries for DFS that is used to find linear descriptors
	boost::default_color_type color;
	int path;
// deduced as part of FCSP algorithm
	int valence; // effective valence
	int piE; // number of PI-electrons
	bool inAromaCycle; // is part of aromatic cycle?
	AtomVertex(){}
	AtomVertex(Code code_):
		code(code_), path(0), valence(0), piE(0), inAromaCycle(false){}
};

struct Bound{
	int type; //
	Bound(){}
	Bound(int type_) :type(type_){}
};

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, AtomVertex, Bound> ChemGraph;

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

enum FCSPFMT {
	JSON, // array of JSON arrays with pairs : (code,bindings)
	CSV, // CSV - 2 columns: file-name,codes
	TXT // TXT - line per file, whitespace separated codes
};

struct FCSPOptions{
	std::vector<LevelOne> first;
	std::vector<LevelTwo> second;
	std::vector<Replacement> replacements;
	bool long41;
	FCSPFMT format;
};

struct FCSP {
	FCSP(FCSPOptions opts);
	void load(std::istream& inp);
	void dumpGraph(std::ostream& dot);
	void process(std::ostream& out, std::string filename="");
	~FCSP();
private:
	struct Impl;
	std::unique_ptr<Impl> pimpl;
};

ChemGraph toGraph(CTab& tab);
void dumpGraph(ChemGraph& graph, std::ostream& out);