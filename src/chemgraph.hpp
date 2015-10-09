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

enum {
	SINGLE = 1,
	DOUBLE,
	TRIPPLE,
	AROMATIC
};

struct Bound{
	int type; //
	Bound(){}
	Bound(int type_) :type(type_){}
};

using ChemGraph = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, AtomVertex, Bound>;
using vd = ChemGraph::vertex_descriptor;
using ed = ChemGraph::edge_descriptor;

ChemGraph toGraph(CTab& tab);
ChemGraph& addHydrogen(ChemGraph& graph);
int getValence(ChemGraph& graph, ChemGraph::vertex_descriptor vertex);
void dumpGraph(ChemGraph& graph, std::ostream& out);

template<class T, class EdgeMap>
std::vector<vd> cycleToChain(std::vector<T>& ic, EdgeMap&& mapper)
{
	typename std::vector<vd> vc;
	auto seed = mapper(ic.front());
	vc.push_back(seed.first);
	vc.push_back(seed.second);
	while (vc.front() != vc.back())
	{
		bool found = false;
		for (auto e : ic)
		{
			auto p = mapper(e);
			//has common vertex with back of chain and not == second one
			if (p.first == vc.back() && p.second != vc[vc.size() - 2])
				vc.push_back(p.second);
			else if (p.second == vc.back() && p.first != vc[vc.size() - 2])
				vc.push_back(p.first);
			//has common vertex with front of chain and not == second one
			else if (p.first == vc[0] && p.second != vc[1])
				vc.insert(vc.begin(), p.second);
			else if (p.second == vc[0] && p.first != vc[1])
				vc.insert(vc.begin(), p.first);
			else
				continue;
			found = true;
			break;
		}
		assert(found);
	}
	vc.pop_back();
	return vc;
}

// transitional helper for Cycle
struct edge_less
{
	bool operator()(std::pair<vd, vd> lhs, std::pair<vd, vd> rhs)
	{
		return lhs.first < rhs.first ||
			(lhs.first == rhs.first && lhs.second < rhs.second);
	}
};


struct Cycle{
public:
	std::vector<vd> chain;
	std::vector<std::pair<vd, vd>> edges; // ordered vertices (first<second)
	bool aromatic_;
public:
	bool aromatic()const{ return aromatic_; }
	// chemical notion of size - number of edges
	size_t size()const{ return edges.size(); }
	Cycle(std::vector<std::pair<vd, vd>> edges_);
	friend std::ostream& operator<<(std::ostream& stream, const Cycle& cycle);
	// sets aromatic flags on atoms and cycle itself iff aromatic
	Cycle& markAromatic(ChemGraph& graph);
};

std::ostream& operator<<(std::ostream& stream, const Cycle& cycle);

// Obtain minimal cycle basis
std::vector<Cycle> minimalCycleBasis(ChemGraph& graph);