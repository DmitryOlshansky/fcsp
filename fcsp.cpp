/*
 * fcsp.cpp
 *
 *  Created on: Dec 16, 2013
 *      Author: dmitry
 */
#include <algorithm>
#include <iomanip>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/graphviz.hpp>

#include "fcsp.hpp"
#include "ctab.h"
#include "descriptors.h"

using namespace std;
using namespace boost;

void printCycle(vector<pair<int, int>>& v)
{
	for (auto & e : v)
	{
		cout << e.first << "--" << e.second << endl;
	}
	cout << endl;
}

inline bool couplingLink(int type)
{
	return type != 1; // not a single link - aromatic or double, triple
}

struct Atom{
	double x,y,z;
	Code code;
	default_color_type color; //used by DFS algorithm
	int path;
	bool inCycle; //part of cycle
	Atom(){}
	Atom(double x_, double y_, double z_, Code code_):
		x(x_), y(y_), z(z_), code(code_), path(0){}
};

struct Bound{
	int type; //
	Bound(){}
	Bound(int type_):type(type_){}
};

typedef adjacency_list<vecS, vecS, undirectedS, Atom, Bound> ChemGraph;

struct markLoops : public dfs_visitor<>{
	vector<pair<int,int>> path;
	vector<vector<pair<int, int>>>& cycles;
	markLoops(vector<vector<pair<int, int>>> &cycles_) :
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
		cout << a << "-->" << b << endl;
	}

	template<class Edge, class Graph>
	void back_edge(Edge e, Graph& g)
	{
		auto a = source(e, g);
		auto b = target(e, g);
		auto c = make_pair((int)b, (int)a);
		cout << a << "++>" << b << endl;
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
			vector<pair<int, int>> v{ start, path.begin()+len };
			v.emplace_back(a, b);
			cout << "~~~~~" << endl;
			printCycle(v);
			cycles.push_back(std::move(v));
			cout << "*****" << endl;
		}
	}
};

bool edge_less(pair<int, int> lhs, pair<int, int> rhs)
{
	return lhs.first == rhs.first ?
		lhs.second < rhs.second : lhs.first < rhs.first;
}


ChemGraph toGraph(CTab& tab)
{
	ChemGraph graph;
	//add_vertex
	for_each(tab.atoms.begin(), tab.atoms.end(), [&graph](const AtomEntry& a)
	{
		add_vertex(Atom(a.x, a.y, a.z, a.code), graph);
	});
	auto vrange = vertices(graph);
	//add edges
	for_each(tab.bounds.begin(), tab.bounds.end(), [&graph, &vrange](const BoundEntry& b)
	{
		add_edge(vrange.first[b.a1] - 1, vrange.first[b.a2] - 1, Bound(b.type), graph);
	});
	return graph;
}

void processCycles(ChemGraph& graph)
{
	vector<vector<pair<int, int>>> cycles;
	depth_first_search(graph, markLoops(cycles), get(&Atom::color, graph));
	for (size_t i = 0; i < cycles.size(); i++)
	for (size_t j = i + 1; j < cycles.size(); j++)
	{
		//split cycles (XOR)
		//pick pick the shortest two of (c1, c2, c1 XOR c2)
		auto& c1 = cycles[i];
		auto& c2 = cycles[j];
		for (auto p = c1.begin(); p != c1.end(); p++)
		{
			auto q = find(c2.begin(), c2.end(), *p);
			if (q != c2.end()) // c1[p] is found in c2 at it
			{
				//compute maximum length of overlap
				auto len = min(c2.end() - q, c1.end() - p);
				//get iterators up to the first non-matching edge
				auto m = mismatch(p, p + len, q);
				//c1[p..m.second] and c2[it..m.first]
				assert(m.first - p == m.second - q);
				len = m.first - p;
				auto lc1 = c1.size() - len;
				auto lc2 = c2.size() - len;
				cout << "C1: " << c1.size() << ", C2: " << c2.size() << ", C1 ^ C2: " << lc1 + lc2 << endl;
				if (c1.size() > c2.size())
				{
					if (c1.size() > lc1 + lc2)
					{
						//kill intersection
						auto pos = p - c1.begin();
						c1.erase(p, p + len);
						//reserve space for lc2
						c1.insert(c1.begin() + pos, lc2, pair<int, int>(-1, -1));
						p = c1.begin() + pos;
						p = reverse_copy(c2.begin(), q, p);
						reverse_copy(q + len, c2.end(), p);
					}
					//no split XOR cycle is the smalest
				}
				else
				{
					//c2 >= c1
					if (c2.size() > lc1 + lc2)
					{
						//split c2
						//kill intersection
						auto pos = q - c2.begin();
						c2.erase(q, q + len);
						//reserve space for lc2
						c2.insert(c2.begin() + pos, lc1, pair<int, int>(-1, -1));
						q = c2.begin() + pos;
						q = reverse_copy(c1.begin(), p, q);
						reverse_copy(p + len, c1.end(), q);

					}
					//no split XOR cycle is the smalest
				}
				break;
			}
		}
	}
	cout << "AFTER SPLIT" << endl;
	for_each(cycles.begin(), cycles.end(), &printCycle);
}

struct TrackPath: public default_bfs_visitor {
	ChemGraph& g;
	bool& coupled;
	ChemGraph::vertex_descriptor tgt;
	TrackPath(ChemGraph& graph, ChemGraph::vertex_descriptor t, bool& c):
		g(graph), coupled(c), tgt(t){ coupled = false; }

	template<typename Vertex>
	void initialize_vertex(Vertex v, const ChemGraph&) const
	{
		//auto a = target(e, g);
		if (g[v].code != C)
			g[v].color = default_color_type::black_color;
		g[v].path = 0;
	}

	template<typename Edge>
	void tree_edge(Edge e, const ChemGraph&) const
	{
		auto d = target(e, g);
		auto s = source(e, g);
		if(couplingLink(g[e].type)){
			coupled = true;
		}
		//cout << s << "-->" << d << endl;
		g[d].path = g[s].path + 1;
	}
};

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

struct FCSP::Impl{
	Impl(std::vector<LevelOne> f, std::vector<LevelTwo> s) :
		order1(std::move(f)), order2(std::move(s)){}

	void load(istream& inp)
	{
		tab = readMol(inp);
		graph = toGraph(tab);
	}

	void dumpGraph(ostream& out)
	{
		write_graphviz(out, graph, CodeWriter(graph));
	}

	void linear(ostream& out)
	{
		using vd = ChemGraph::vertex_descriptor;
		vector<pair<vd, int>> dcs;
		auto vrtx = vertices(graph);
		for (auto i = vrtx.first; i != vrtx.second; i++)
		{
			//TODO: pattern match 2nd order
			auto edges = out_edges(*i, graph);
			int valency = 0;
			for (auto p = edges.first; p != edges.second; p++)
			{
				valency += graph[*p].type;
			}
			cout << "VALENCY " << valency << " for " << graph[*i].code.symbol() << endl;
			LevelOne t(graph[*i].code, valency, 0);
			auto range = equal_range(order1.begin(), order1.end(), t);
			auto j = range.first;
			for (; j != range.second; ++j)
			{
				//out << *i << ": " << atomSymbol(j->second.center)
				//		<< " DC:" << j->second.index << endl;
				dcs.emplace_back(*i, j->dc);
			}
		}
		out << endl;
		sort(dcs.begin(), dcs.end(), [](const pair<vd, int>& a, const pair<vd, int>& b){
			return a.second < b.second;
		});

		queue<ChemGraph::vertex_descriptor> buf;
		for (size_t i = 0; i < dcs.size(); i++)
		for (size_t j = i  + 1; j < dcs.size(); j++)
		{
			bool coupled = false;
			breadth_first_search(graph, dcs[j].first, buf,
					TrackPath(graph, dcs[i].first, coupled), get(&Atom::color, graph));
			auto vtx = vertices(graph);
			if(graph[dcs[i].first].path){ // not zero - connected
				/*{
					for (auto i = vtx.first; i != vtx.second; i++)
						cout << *i << ": " << graph[*i].path << endl;
				}*/
				out << setfill('0') << setw(2) << dcs[i].second
					<< setfill('0') << setw(2) << graph[dcs[i].first].path - 1
					<< setfill('0') << setw(2) << dcs[j].second
					<< (coupled ? 1 : 0) << endl;
			}
		}
	}

	void cyclic(ostream& out)
	{

	}
private:
	std::vector<LevelOne> order1;
	std::vector<LevelTwo> order2;
	ChemGraph graph;
	CTab tab;
};

FCSP::FCSP(std::vector<LevelOne> first, std::vector<LevelTwo> second) :
	pimpl(new FCSP::Impl(std::move(first), std::move(second))){}

void FCSP::load(std::istream& inp)
{
	pimpl->load(inp);
}

void FCSP::dumpGraph(std::ostream& dot)
{
	pimpl->dumpGraph(dot);
}

void FCSP::linear(std::ostream& out)
{
	pimpl->linear(out);
}

FCSP::~FCSP(){}
