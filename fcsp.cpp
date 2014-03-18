/*
 * fcsp.cpp
 *
 *  Created on: Dec 16, 2013
 *      Author: dmitry
 */
#include <algorithm>
#include <iomanip>
#include <set>
#include "ullmann.h"
#include "fcsp.hpp"
#include "ctab.h"
#include "descriptors.h"

enum { NON_PASSABLE = 10000 };
using namespace std;
using namespace boost;

using vd = ChemGraph::vertex_descriptor;

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

struct edge_less
{
	bool operator()(pair<int, int> lhs, pair<int, int> rhs)
	{
		return lhs.first < rhs.first ||
			(lhs.first == rhs.first && lhs.second < rhs.second);
	}
};

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

struct TrackPath: public default_bfs_visitor {
	ChemGraph& g;
	ChemGraph::vertex_descriptor start, tgt;
	TrackPath(ChemGraph& graph, ChemGraph::vertex_descriptor s, ChemGraph::vertex_descriptor t) :
		g(graph), start(s), tgt(t){ }

	template<typename Vertex>
	void initialize_vertex(Vertex v, const ChemGraph&) const
	{
		//auto a = target(e, g);
		g[v].path = 0;
	}

	template<typename Edge>
	void tree_edge(Edge e, const ChemGraph&) const
	{
		auto d = target(e, g);
		auto s = source(e, g);
		//cout << s << "-->" << d << endl;
		if (g[d].inAromaCycle && d != tgt) //all paths  of aromatic cycle  are inpassable
			g[d].path = NON_PASSABLE;
		else if (g[d].inAromaCycle && g[s].inAromaCycle)
			g[d].path = NON_PASSABLE;
		else if (g[d].code != C && d != tgt)
			g[d].path = NON_PASSABLE;
		else
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

pair<int,int> multiCount(ChemGraph& graph, ChemGraph::vertex_descriptor vertex)
{
	auto edges = out_edges(vertex, graph);
	int dual = 0, tripple = 0;
	for (auto p = edges.first; p != edges.second; p++)
	{
		if (graph[*p].type == 2)
			dual++;
		else if (graph[*p].type == 3)
			tripple++;
	}
	return make_pair(dual, tripple);
}

struct FCSP::Impl{
	Impl(std::vector<LevelOne> f, std::vector<LevelTwo> s, std::vector<Replacement> r) :
		order1(std::move(f)), order2(std::move(s)), repls(std::move(r)){}

	void load(istream& inp)
	{
		CTab tab = readMol(inp);
		graph = toGraph(tab);
		//TODO: add implicit Hydrogen!!! (but only for FCSP molecule, not in toGraph)
		addHydrogen(tab, graph);
	}

	void process(ostream& out)
	{
		locateDCs();
		locateCycles();
		cyclic(out);
		linear(out);
		replacement(out);
	}

	void dumpGraph(ostream& out)
	{
		::dumpGraph(graph, out);
	}

	void locateDCs()
	{
		dcs.clear();
		using vd = ChemGraph::vertex_descriptor;
		auto vrtx = vertices(graph);
		for (auto i = vrtx.first; i != vrtx.second; i++)
		{
			auto edges = out_edges(*i, graph);
			//pre-calculate per-atom properties
			int valency = getValence(graph, *i);
			graph[*i].valence = valency;
			auto dual_tripple = multiCount(graph, *i);
			int piE = countPiElectrons(graph[*i].code, valency, dual_tripple.first, dual_tripple.second);
			graph[*i].piE = piE > 0 ? piE : 0;
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
				if (j->valence != valency)
					continue;
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
	}

	static set<int> cycleToVertices(const vector<pair<int, int>>& cycle)
	{
		set<int> vertices;
		for (auto e : cycle)
		{
			vertices.insert(e.first);
			vertices.insert(e.second);
		}
		return vertices;
	}

	static int chainAt(const vector<int>& chain, int idx)
	{
		if (idx < 0)
			idx += (int)chain.size();
		if (idx >(int)chain.size())
			idx -= (int)chain.size();
		return idx;
	}

	static bool heteroatom(const Atom& atom)
	{
		return !atom.code.matches(C) && !atom.code.matches(H);
	}

	void locateCycles()
	{
		cycles.clear();
		depth_first_search(graph, markLoops(cycles), get(&Atom::color, graph));
		for (auto & c : cycles)
		{
			sort(c.begin(), c.end(), edge_less());
		}
		cout << "BEFORE SPLIT" << endl;
		for_each(cycles.begin(), cycles.end(), &printCycle);
		vector<pair<int, int>> t, uc, xc;
		//TODO: need some solid proof and potentially incorrect in complex cases:
		// what happens after 2 cycles are replaced with XOR or U of original pair?
		for (size_t i = 0; i < cycles.size(); i++)
		for (size_t j = i + 1; j < cycles.size(); j++)
		{
			auto& c1 = cycles[i];
			auto& c2 = cycles[j];
			t.resize(min(c1.size(), c2.size()));
			auto tend = set_intersection(c1.begin(), c1.end(), c2.begin(), c2.end(), t.begin(), edge_less());
			if (t.begin() == tend) // empty intersection
				continue;
			uc.resize(c1.size() + c2.size());
			xc.resize(c1.size() + c2.size());
			auto ucend = set_union(c1.begin(), c1.end(), c2.begin(), c2.end(), uc.begin(), edge_less());
			auto xcend = set_symmetric_difference(c1.begin(), c1.end(), c2.begin(), c2.end(), xc.begin(), edge_less());
			uc.resize(ucend - uc.begin());
			xc.resize(xcend - xc.begin());
			size_t len1 = c1.size(), len2 = c2.size();
			size_t lenU = uc.size(), lenX = xc.size();
			pair<size_t, vector<pair<int, int>>*> result[4];
			result[0] = make_pair(len1, &c1);
			result[1] = make_pair(len2, &c2);
			result[2] = make_pair(lenU, &uc);
			result[3] = make_pair(lenX, &xc);
			for (int k = 0; k < 4; k++)
			{
				cout << result[k].first << " ";
			}
			cout << endl;
			sort(result, result + 4, [](pair<size_t, vector<pair<int, int>>*> a, pair<size_t, vector<pair<int, int>>*> b){
				return a.first < b.first;
			});

			if ((&c1 == result[0].second && &c2 == result[1].second) ||
				(&c1 == result[1].second && &c2 == result[0].second))
				continue; //same cycles have won
			auto n1 = *result[1].second;
			auto n2 = *result[0].second;
			cout << n1.size() << " " << n2.size() << endl;
			c1 = n1;
			c2 = n2;
		}
		cout << "AFTER SPLIT" << endl;
		for (auto ic : cycles)
		{
			assert(ic.size() > 2);
			vector<int> vc;
			vc.push_back(ic.front().first);
			vc.push_back(ic.front().second);
			while (vc.front() != vc.back())
			{
				for (auto p : ic)
				{
					//has common vertex with back of chain and not == second one
					if (p.first == vc.back() && p.second != vc[vc.size()-2])
						vc.push_back(p.second);
					else if (p.second == vc.back() && p.first != vc[vc.size()-2])
						vc.push_back(p.first);
					//has common vertex with front of chain and not == second one
					else if(p.first == vc[0] && p.second != vc[1])
						vc.insert(vc.begin(), p.second);
					else if (p.second == vc[0] && p.first != vc[1])
						vc.insert(vc.begin(), p.first);
					else
						continue;
					break;
				}
			}
			vc.pop_back();
			cout << "CHAIN: ";
			for (int k : vc)
				cout << k << " - ";
			cout << endl;
			auto& g = graph;
			int piEl = 0;
			for_each(vc.begin(), vc.end(), [&g, &piEl](int n){
				cout << "Atom # " << n << " " << g[n].code.symbol() << " pi E = " << g[n].piE << endl;
				if (g[n].piE > 0)
					piEl += g[n].piE;
				//TODO: add debug trace for < 0
			});
			bool aromatic = ((piEl - 2) % 4 == 0);
			cout << "PI E " << piEl << " aromatic? : " << aromatic << endl;
			//assign cyclic-only DC
			if (aromatic)
			{
				for (auto n : vc)
				{
					g[n].inAromaCycle = true;
					//any atom in aromatic cycle
					dcs.emplace_back(n, 33);
					//hetero-atom in aromatic cycle - DC 34
					if (!g[n].code.matches(C))
						dcs.emplace_back(n, 34);
					auto idx = find(vc.begin(), vc.end(), n) - vc.begin();
					assert(idx < vc.size());
					if (heteroatom(graph[chainAt(vc, idx - 1)])
						|| heteroatom(graph[chainAt(vc, idx + 1)]))
					{
						dcs.emplace_back(n, 35);
					}
					else if (heteroatom(graph[chainAt(vc, idx - 2)])
						|| heteroatom(graph[chainAt(vc, idx + 2)]))
					{
						dcs.emplace_back(n, 36);
					}
					else if (heteroatom(graph[chainAt(vc, idx - 3)])
						|| heteroatom(graph[chainAt(vc, idx + 3)]))
					{
						dcs.emplace_back(n, 37);
					}
					//TODO: DC 40 - need reference table of valence, to test charged atoms
				}
			}
			using std::move;
			chains.push_back(move(vc));
		}
		typedef ChemGraph::vertex_descriptor vd;
		sort(dcs.begin(), dcs.end(), [](const pair<vd, int>& a, const pair<vd, int>& b){
			return a.second < b.second;
		});
		cout << "DCs:" << endl;
		for (auto& dc : dcs)
		{
			cout << dc.first << " -DC-> " << dc.second << endl;
		}
		cout << endl;
	}

	template<class Fn>
	void applyPath(vd start, vd end, Fn&& fn)
	{
		vd current = end;
		fn(current);
		while (current != start)
		{
			//cout << graph[current].path << " -->" << endl;
			auto nearby = adjacent_vertices(current, graph);
			auto i = nearby.first;
			for (; i != nearby.second; i++)
			{
				//cout << graph[*i].path << ".." ;
				if (graph[*i].path == graph[current].path - 1)
				{
					current = *i;
					fn(current);
					break;
				}
			}
			if (i == nearby.second)
				break;
		}
	}

	void linear(ostream& out)
	{
		
		queue<vd> buf;
		for (size_t i = 0; i < dcs.size(); i++)
		for (size_t j = i  + 1; j < dcs.size(); j++)
		{
			vd start = dcs[j].first;
			vd end = dcs[i].first;
			breadth_first_search(graph, start, buf,
				TrackPath(graph, start, end), get(&Atom::color, graph));
			
			if (graph[dcs[i].first].path && graph[dcs[i].first].path < NON_PASSABLE)
			{
				bool coupled = false;
				auto &g = graph;
				
				applyPath(start, end, [g, &coupled](vd v){
					//cout << "Apply: " << g[v].code.symbol() << endl;
					if (g[v].code.matches(C) && g[v].piE > 0)
						coupled = true;

				});
				out << setfill('0') << setw(2) << dcs[i].second
					<< setfill('0') << setw(2) << graph[dcs[i].first].path - 1
					<< setfill('0') << setw(2) << dcs[j].second
					<< (coupled ? 1 : 0) << endl;
			}
		}
	}

	void cyclic(ostream& out)
	{
		//write out simple cycles
		for (auto c : chains)
		{
			int piE = 0;
			if (graph[c.front()].inAromaCycle)
			{
				for (int v : c)
					piE += graph[v].piE;
			}
			cout << setfill('0') << setw(1) << c.size()
				<< "," 
				<< setfill('0') << setw(2) << piE
				<< endl;
		}
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
	//location of DCs in 'graph' and their numeric value
	vector<pair<ChemGraph::vertex_descriptor, int>> dcs;
	//sorted arrays of edges - basic cycles
	vector<vector<pair<int, int>>> cycles;
	//same basic cycles represented as chains of vertices
	vector<vector<int>> chains;
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

void FCSP::process(std::ostream& out)
{
	pimpl->process(out);
}

FCSP::~FCSP(){}
