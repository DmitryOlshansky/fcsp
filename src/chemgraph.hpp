#pragma once

#include <ostream>
#include <boost/graph/adjacency_list.hpp>
#include "periodic.hpp"
#include "ctab.hpp" // MOL file format (aka CTable)

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

using ChemGraph = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, AtomVertex, Bound>;

ChemGraph toGraph(CTab& tab);
void dumpGraph(ChemGraph& graph, std::ostream& out);
