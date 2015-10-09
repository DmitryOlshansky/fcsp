#include <boost/graph/graphviz.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/visitors.hpp>
#include "chemgraph.hpp"
#include "ctab.hpp"
#include "log.hpp"

using namespace std;
using namespace boost;

ChemGraph toGraph(CTab& tab)
{
	ChemGraph graph;
	//add_vertex
	for_each(tab.atoms.begin(), tab.atoms.end(), [&graph](const AtomEntry& a)
	{
		add_vertex(AtomVertex(a.code), graph);
	});
	auto vrange = vertices(graph);
	//add edges
	for_each(tab.bounds.begin(), tab.bounds.end(), [&graph, &vrange](const BoundEntry& b)
	{
		add_edge(vrange.first[b.a1] - 1, vrange.first[b.a2] - 1, Bound(b.type), graph);
	});
	return graph;
}


class CodeWriter {
public:
	CodeWriter(ChemGraph& graph):g(graph){}
	template <class VertexOrEdge>
	void operator()(std::ostream& out, const VertexOrEdge& v) const {
	  out << "[label=\"" << g[v].code.symbol() << " [" << v << "]" << "\"]";
	}
private:
	ChemGraph& g;
};


void dumpGraph(ChemGraph& graph, ostream& out)
{
	write_graphviz(out, graph, CodeWriter(graph));
}

ostream& operator<<(ostream& stream, const Cycle& cycle)
{
	for (auto & e : cycle.edges)
	{
		stream << e.first << "--" << e.second << '\n';
	}
	return stream;
}


void logCycle(vector<pair<vd,vd>>& c)
{
	for (auto & e : c)
	{
		LOG(DEBUG) << e.first << "--" << e.second << '\n';
	}
	LOG(DEBUG) << endline;
}


struct markLoops : public dfs_visitor<>{
	vector<pair<vd,vd>> path;
	vector<vector<pair<vd, vd>>>& cycles;
	markLoops(vector<vector<pair<vd, vd>>> &cycles_) :
		cycles(cycles_){}

	template<class Edge, class Graph>
	void tree_edge(Edge e, Graph& g)
	{
		auto a = source(e, g);
		auto b = target(e, g);
		//back-track to the start of this edge
		while (!path.empty() && path.back().second != a){
			path.pop_back();
		}
		path.emplace_back(a, b);
		LOG(TRACE) << a << "-->" << b << endline;
	}

	template<class Edge, class Graph>
	void back_edge(Edge e, Graph& g)
	{
		auto a = source(e, g);
		auto b = target(e, g);
		auto c = make_pair(b, a);
		LOG(TRACE) << a << "++>" << b << endline;
		if (find_if(path.begin(), path.end(), [c](const pair<int, int> &p){
			return p.first == c.first && p.second == c.second;
		}) == path.end())
		{
			auto end = find_if(path.rbegin(), path.rend(), [a](const pair<int, int> &p){
				return p.second == a;
			});
			auto len = path.rend() - end;
			auto start = find_if(path.begin(), path.begin() + len, [b](const pair<int, int> &p){
				return p.first == b;
			});
			vector<pair<vd, vd>> v{ start, path.begin()+len };
			v.emplace_back(a, b);
			LOG(DEBUG) << "~~~~~" << endline;
			logCycle(v);
			cycles.push_back(std::move(v));
			LOG(DEBUG) << "*****" << endline;
		}
	}
};



// just array of sorted pairs
vector<vector<pair<vd,vd>>> cycleBasis(ChemGraph& graph)
{
	vector<vector<pair<vd,vd>>> cycles;
	depth_first_search(graph, markLoops(cycles), get(&AtomVertex::color, graph));
	for (auto & c : cycles)
	{
		transform(c.begin(), c.end(), c.begin(), [](const pair<vd, vd> &p){ 
			return p.first < p.second ? p : make_pair(p.second, p.first);
		});
		sort(c.begin(), c.end(), edge_less());
	}
	return cycles;
}

Cycle::Cycle(vector<pair<vd, vd>> edges_):edges(edges_)
{
	chain = cycleToChain(edges, [](const pair<vd,vd>&p){ return p; });
}

Cycle& Cycle::markAromatic(ChemGraph& g)
{
	LOG(TRACE) << "CHAIN: ";
	for (int k : chain)
		LOG(TRACE) << k << " - ";
	LOG(TRACE) << endline;
	int piEl = 0;
	for_each(chain.begin(), chain.end(), [&g, &piEl](int n){
		LOG(TRACE) << "Atom # " << n << " " << g[n].code.symbol() << " pi E = " << g[n].piE << endline;
		if (g[n].piE > 0)
			piEl += g[n].piE;
		//TODO: add debug trace for < 0
	});
	aromatic_ = ((piEl - 2) % 4 == 0); // Hukkel rule 4n + 2
	LOG(TRACE) << "PI E " << piEl << " aromatic? : " << aromatic_ << endline;
	//assign cyclic-only DC	
	if (aromatic_)
	{
		for (auto n : chain)
		{
			g[n].inAromaCycle = true;
		}
	}
	return *this;
}

vector<Cycle> minimalCycleBasis(ChemGraph& graph)
{
	auto cycles = cycleBasis(graph);
	LOG(DEBUG) << "BEFORE SPLIT" << endline;
	for_each(cycles.begin(), cycles.end(), &logCycle);
	vector<pair<vd, vd>> t, t2, uc, xc;
	//TODO: need some solid proof and potentially incorrect in complex cases:
	// what happens after 2 cycles are replaced with XOR or U of original pair?
	for (int q = 0; q < 10; q++)//hack
	for (size_t i = 0; i < cycles.size(); i++)
	for (size_t j = i + 1; j < cycles.size(); j++)
	{
		auto& c1 = cycles[i];
		auto& c2 = cycles[j];
		t.resize(min(c1.size(), c2.size()));
		t2.resize(max(c1.size(), c2.size()));
		auto tend = set_intersection(c1.begin(), c1.end(), c2.begin(), c2.end(), t.begin(), edge_less());
		if (t.begin() == tend) // empty intersection
			continue;
		uc.resize(c1.size() + c2.size());
		xc.resize(c1.size() + c2.size());
		auto xcend = set_symmetric_difference(c1.begin(), c1.end(), c2.begin(), c2.end(), xc.begin(), edge_less());
		//uc.resize(ucend - uc.begin());
		xc.resize(xcend - xc.begin());
		size_t len1 = c1.size(), len2 = c2.size();
		size_t /*lenU = uc.size(),*/ lenX = xc.size();
		pair<size_t, vector<pair<vd, vd>>*> result[4];
		result[0] = make_pair(len1, &c1);
		result[1] = make_pair(len2, &c2);			
		int m = 2;
		/*if (uc != c1)
			result[m++] = make_pair(lenU, &uc);
		else*/ if (xc != c2)
			result[m++] = make_pair(lenX, &xc);
		for (int k = 0; k < m; k++)
		{
			LOG(TRACE) << result[k].first << " ";
		}
		LOG(TRACE) << endline;
		sort(result, result + m, [](pair<size_t, vector<pair<vd, vd>>*> a, pair<size_t, vector<pair<vd, vd>>*> b){
			return a.first < b.first;
		});

		if ((&c1 == result[0].second && &c2 == result[1].second) ||
			(&c1 == result[1].second && &c2 == result[0].second))
			continue; //same cycles have won
		auto n1 = *result[1].second;
		auto n2 = *result[0].second;
		LOG(TRACE) << n1.size() << " " << n2.size() << endline;
		c1 = n1;
		c2 = n2;
	}
	LOG(DEBUG) << "AFTER SPLIT" << endline;
	for_each(cycles.begin(), cycles.end(), &logCycle);
	vector<Cycle> minCycles;
	transform(cycles.begin(), cycles.end(), back_inserter(minCycles), [](vector<pair<vd,vd>>& cyc){
		return Cycle(cyc);
	});
	return minCycles;
}