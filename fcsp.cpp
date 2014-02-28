/*
 * fcsp.cpp
 *
 *  Created on: Dec 16, 2013
 *      Author: dmitry
 */
#include <algorithm>
#include <iomanip>
#include "ullmann.h"
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

void addHydrogen(CTab& tab, ChemGraph& graph)
{	
	/*auto c = Code("C");
	auto h = Code("H");
	auto& atoms = tab.atoms;
	for (size_t i = 0; i < atoms.size(); i++)
	{
		
			//for (int k = 0; k < atoms[i].implicitH; k++)
			{
				auto j = add_vertex(Atom(0.0, 0.0, 0.0, h), graph);
				add_edge(j, i, Bound(1), graph);
			}
		}
	}*/
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

void dumpGraph(ChemGraph& graph, ostream& out)
{
	write_graphviz(out, graph, CodeWriter(graph));
}

int getValence(ChemGraph& graph, ChemGraph::vertex_descriptor vertex)
{
	auto edges = out_edges(vertex, graph);
	int valence = 0;
	for (auto p = edges.first; p != edges.second; p++)
	{
		valence += graph[*p].type;
	}
	return valence;
}

struct FCSP::Impl{
	Impl(std::vector<LevelOne> f, std::vector<LevelTwo> s, std::vector<Replacement> r) :
		order1(std::move(f)), order2(std::move(s)), repls(std::move(r)){}

	void load(istream& inp)
	{
		CTab tab = readMol(inp);
		graph = toGraph(tab);
		addHydrogen(tab, graph);
	}

	void dumpGraph(ostream& out)
	{
		::dumpGraph(graph, out);
	}

	void locateDCs()
	{
		//TODO: add implicit Hydrogen!!! (but only for FCSP molecule, not in toGraph)
		dcs.clear();
		using vd = ChemGraph::vertex_descriptor;
		auto vrtx = vertices(graph);
		for (auto i = vrtx.first; i != vrtx.second; i++)
		{
			auto edges = out_edges(*i, graph);
			int valency = 0;
			for (auto p = edges.first; p != edges.second; p++)
			{
				valency += graph[*p].type;
			}
			cout << "VALENCY " << valency << " for " << graph[*i].code.symbol() << endl;
			LevelOne t(graph[*i].code, valency, 0);
			auto range = equal_range(order1.begin(), order1.end(), t);
			for (auto j = range.first; j != range.second; ++j)
			{
				//out << *i << ": " << atomSymbol(j->second.center)
				//		<< " DC:" << j->second.index << endl;
				dcs.emplace_back(*i, j->dc);
			}
			auto range2 = make_pair(order2.begin(), order2.end());

			for (auto j = range2.first; j != range2.second; j++)
			{
				if (edges.second - edges.first != j->bonds.size())
					continue;
				if (!j->center.matches(graph[*i].code.code()))
					continue;
				cout << "Matched CENTER " << j->center.symbol() << endl;
				//if (j->valence != valency)
				//	continue;
				vector<bool> matched(j->bonds.size());
				auto p = edges.first;
				for (; p != edges.second; p++)
				{
					size_t k;
					for (k = 0; k < j->bonds.size(); k++)
					{
						if (matched[k])
							continue;
						if (graph[*p].type != j->bonds[k].bondType)
							continue;
						auto v = target(*p, graph);
						if (!j->bonds[k].atom.matches(graph[v].code.code()))
							continue;
						cout << "   matched " << j->bonds[k].atom.symbol()
							<< " vs " << graph[v].code.symbol() << endl;
						matched[k] = true;
						break;
					}
					if (k == j->bonds.size())
						break;
				}
				if (p == edges.second)
				{
					cout << "DONE." << endl;
					dcs.emplace_back(*i, j->dc);
				}
			}
		}
		cout << endl;
		sort(dcs.begin(), dcs.end(), [](const pair<vd, int>& a, const pair<vd, int>& b){
			return a.second < b.second;
		});
		for (auto& dc : dcs)
		{
			cout << dc.second << endl;
		}
	}

	void linear(ostream& out)
	{
		using vd = ChemGraph::vertex_descriptor;
		queue<vd> buf;
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

	//find where repPos is mapped in mapping m
	int mapVertex(size_t repPos, const vector<pair<size_t, size_t>>& m)
	{
		//locate mapping of 1st DC
		auto dcVIt = find_if(m.begin(), m.end(), [repPos](const pair<size_t, size_t>& p){
			return p.first == repPos;
		});
		assert(dcVIt != m.end()); // mapping must include it
		return dcVIt->second;
	}

	// -1 if no descriptor present (just some non-H atom)
	int mapDC(size_t dcV)
	{
		auto dcIt = find_if(dcs.begin(), dcs.end(), [dcV](const pair<size_t, int> & p){
			return p.first == dcV;
		});
		if (dcIt == dcs.end())
		{
			//debug warning
			cerr << "Mapped replacement but first DC has no match" << endl;
			return -1;
		}
		return dcIt->second;
	}

	void replacement(ostream& out)
	{
		//cout << "REPLACEMENTS!" << endl;
		for (auto& r : repls)
		{
			vector<vector<pair<size_t, size_t>>> mappings;
			auto& g = graph;
			ullmann_all(r.piece, graph, [&r, &g](ChemGraph::vertex_descriptor a, ChemGraph::vertex_descriptor b){
				return r.piece[a].code.matches(g[b].code.code());
			}, [&r, &g](ChemGraph::edge_descriptor a, ChemGraph::edge_descriptor b){
				return r.piece[a].type == g[b].type;
			}, mappings);
			vector<pair<size_t, size_t>> used_pairs;
			for (auto& m : mappings)
			{
				int fV = mapVertex(r.a1, m);
				int sV = mapVertex(r.a2, m);
				// check if this pair of vertices was used before
				if (used_pairs.end() != find_if(used_pairs.begin(), used_pairs.end(), [fV, sV](const pair<size_t, size_t>& p){
					return (p.first == fV && p.second == sV) || (p.first == sV && p.second == fV);
				}))
					continue;
				int fdc = mapDC(fV);
				if (fdc < 0)
					continue;
				int sdc = mapDC(sV);
				if (sdc < 0)
					continue;
				used_pairs.emplace_back(fV, sV);
				// the usual rule of smaller DC first
				if (fdc > sdc)
					swap(fdc, sdc);
				cout << setfill('0') << setw(2) << fdc
					<< setfill('0') << setw(2) << r.dc
					<< setfill('0') << setw(2) << sdc
					<< r.coupling << endl;
			}
		}
	}
private:
	std::vector<LevelOne> order1;
	std::vector<LevelTwo> order2;
	std::vector<Replacement> repls;
	ChemGraph graph;
	vector<pair<ChemGraph::vertex_descriptor, int>> dcs;
};

FCSP::FCSP(std::vector<LevelOne> first, std::vector<LevelTwo> second, std::vector<Replacement> repls) :
	pimpl(new FCSP::Impl(std::move(first), std::move(second), std::move(repls))){}

void FCSP::load(std::istream& inp)
{
	pimpl->load(inp);
}


void FCSP::dumpGraph(std::ostream& dot)
{
	pimpl->dumpGraph(dot);
}

void FCSP::locateDCs()
{
	pimpl->locateDCs();
}

void FCSP::cyclic(std::ostream& out)
{
	pimpl->cyclic(out);
}

void FCSP::replacement(ostream& out)
{
	pimpl->replacement(out);
}

void FCSP::linear(std::ostream& out)
{
	pimpl->linear(out);
}

FCSP::~FCSP(){}
