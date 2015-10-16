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

template<typename V, typename E>
using Graph = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, V, E>;

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


template<class T>
inline bool operator<(const std::pair<T,T>& a, std::pair<T,T>& b){
    return a.first < b.first || (a.first == b.first && a.second < b.second);
}

// transitional helper for Cycle
struct edge_less
{
	bool operator()(std::pair<vd, vd> lhs, std::pair<vd, vd> rhs)
	{
		return lhs < rhs;
	}
};


struct Cycle{
public:
	std::vector<vd> chain;
	std::vector<std::pair<vd, vd>> edges; // ordered vertices (first<second)
	bool aromatic_;
public:
	// From set of edges
	Cycle(std::vector<std::pair<vd, vd>> edges_);
	// From chain of vertices
	Cycle(std::vector<vd> chain);
	// True if there is intersection of this and that cycle
	bool intersects(const Cycle& that)const;
	// Get set of edges of the intersection of this and that
	std::vector<std::pair<vd,vd>> intersection(const Cycle& c)const;
	// True if this cylce is aromatic
	bool aromatic()const{ return aromatic_; }
	// Chemical notion of size - number of edges
	size_t size()const{ return edges.size(); }
	// sets aromatic flags on atoms and cycle itself iff aromatic
	Cycle& markAromatic(ChemGraph& graph);
	// output as chain
	friend std::ostream& operator<<(std::ostream& stream, const Cycle& cycle);
};

std::ostream& operator<<(std::ostream& stream, const Cycle& cycle);

// Obtain minimal cycle basis
std::vector<Cycle> minimalCycleBasis(ChemGraph& graph);

// Impl class
template<class Vertex, class Edge>
struct ConnectedComponents{
    using G = Graph<Vertex,Edge>;
    ConnectedComponents(G& graph):
    components(), g(graph),labels(boost::num_vertices(graph)), label(0){
        findComponents();
    }
    std::vector<G> components;
private:
    void findComponents(){
    		using namespace boost;
    		auto const V = num_vertices(g);
        for(size_t v=0; v<V; v++){
            if(!labels[v]){
                ++label; //start new component
                dfs(v);
            }
        }
        if(label == 1){
            components.push_back(g);
            return;
        }
        components.resize(label);
        // for each component map global --> component
        std::vector<std::vector<size_t>> g2c(label);
        for(auto & vec : g2c)
            vec.resize(V);

        // map vertices to components and record positions
        for(size_t v=0; v<V; v++){
            int lbl = labels[v] - 1;
            g2c[lbl][v] = add_vertex(g[v], components[lbl]);
        }
        // use map to copy over edges
        for(size_t v=0; v<V; v++){
        		auto adj = adjacent_vertices(v, g);
            for(auto p=adj.first; p!=adj.second; p++){
            		auto w = *p;
                if (w < v){ //deduplicate
                    continue;
                }
                int lbl = labels[v] - 1;
                auto e = edge(v, w, g).first;
                add_edge(g2c[lbl][v], g2c[lbl][w], g[e], components[lbl]);
            }
        }
    }

    void dfs(size_t v){
        labels[v] = label;
        auto adj = adjacent_vertices(v, g);
        for(auto p = adj.first; p != adj.second; p++){
        		auto w = *p;
            if(!labels[w]){
                dfs(w);
            }
        }
    }
    G& g;
    std::vector<int> labels;
    int label;
};


// Will copy g if it's the only component
template<class V, class E>
std::vector<Graph<V,E>> connectedComponents(Graph<V,E> &g)
{
    return ConnectedComponents<V,E>(g).components;
}

