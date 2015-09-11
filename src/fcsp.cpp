/*
 * fcsp.cpp
 *
 *  Created on: Dec 16, 2013
 *      Author: dmitry
 */
#include <algorithm>
#include <numeric>
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
using ed = ChemGraph::edge_descriptor;

void printCycle(vector<pair<int, int>>& v)
{
	for (auto & e : v)
	{
		cerr << e.first << "--" << e.second << endl;
	}
	cerr << endl;
}

void printJsArray(vector<int>& v, ostream& out)
{
	out << "[";
	sort(v.begin(), v.end());
	v.erase(unique(v.begin(), v.end()), v.end());
	for(size_t i=0;i<v.size();i++)
	{
		if(i != 0) out << ", ";
		out << v[i];
	}
	out << "]";
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
		cerr << a << "-->" << b << endl;
	}

	template<class Edge, class Graph>
	void back_edge(Edge e, Graph& g)
	{
		auto a = source(e, g);
		auto b = target(e, g);
		auto c = make_pair((int)b, (int)a);
		cerr << a << "++>" << b << endl;
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
			cerr << "~~~~~" << endl;
			printCycle(v);
			cycles.push_back(std::move(v));
			cerr << "*****" << endl;
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

struct TrackPath: public default_bfs_visitor {
	ChemGraph& g;
	vector<pair<vd, int>>& dcs;
	ChemGraph::vertex_descriptor start, tgt;
	bool pass_4546;
	TrackPath(ChemGraph& graph, ChemGraph::vertex_descriptor s, ChemGraph::vertex_descriptor t, 
		vector<pair<vd, int>>& dcsArr) :
		g(graph), start(s), tgt(t), dcs(dcsArr)
	{
		pass_4546 = find_if(dcs.begin(), dcs.end(), [&](pair<vd,int> p){
			return (p.first == start || p.first == tgt) && (p.second == 45 || p.second == 46);
		}) == dcs.end();
	}

	template<typename Vertex>
	void initialize_vertex(Vertex v, const ChemGraph&) const
	{
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
		else if(!pass_4546 && g[d].code == C)
		{
			bool hit4546 = find_if(dcs.begin(), dcs.end(), [d](pair<vd,int> p){
				return p.first == d && (p.second == 45 || p.second == 46);
			}) != dcs.end();
			if(hit4546)
				g[d].path = NON_PASSABLE;
			else
				g[d].path = g[s].path + 1;
		}
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

int singleCount(ChemGraph& graph, ChemGraph::vertex_descriptor vertex)
{
	auto edges = out_edges(vertex, graph);
	int cnt = 0;
	for (auto p = edges.first; p != edges.second; p++)
	{
		if (graph[*p].type == 1)
			cnt++;		
	}
	return cnt;
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
	}
	// Очистить все переменные состояния кодировщика
	void clear()
	{
		jsons.clear();
		dcs.clear();
		dcsAtoms.clear();
		cycles.clear();
		chains.clear();
	}

	void process(ostream& out)
	{
		clear();
		locatePiElectrons();
		locateDCs();
		locateCycles();
		cyclic(out);
		linear(out);
		replacement(out);
		out << "[";
		for(size_t i=0; i<jsons.size();i++){
			if(i != 0)
				out << ", ";
			out << jsons[i];
		}
		out << "]";
		out << endl;
	}

	void dumpGraph(ostream& out)
	{
		::dumpGraph(graph, out);
	}

	void outputJs(string code, vector<int>& atoms)
	{
		stringstream s;
		s << "{" << "\"code\" : \""<<code<<"\",";
		s << "\"place\": ";
		printJsArray(atoms, s); 
		s << "}";
		jsons.push_back(s.str());
	}

	void locatePiElectrons()
	{
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
		}
	}

	void locateDCs()
	{
		auto vrtx = vertices(graph);
		for (auto i = vrtx.first; i != vrtx.second; i++)
		{
			int valency = graph[*i].valence;
			auto edges = out_edges(*i, graph);
			LevelOne t(graph[*i].code, valency, 0);
			auto range = equal_range(order1.begin(), order1.end(), t);
			for (auto j = range.first; j != range.second; ++j)
			{
				//out << *i << ": " << atomSymbol(j->second.center)
				//		<< " DC:" << j->second.index << endl;
				dcs.emplace_back(*i, j->dc);
			}
			if(graph[*i].code == C && !graph[*i].inAromaCycle)
			{
				// check for 45 & 46
				auto edges = out_edges(*i, graph);
				for (auto p = edges.first; p != edges.second; p++)
				{
					auto tgt = target(*p, graph);
					if (graph[*p].type == 2 && graph[tgt].code == C && !graph[tgt].inAromaCycle)
					{
						auto tgt_edges = out_edges(tgt, graph);
						auto cnt = count_if(tgt_edges.first, tgt_edges.second, [&](ed edge){
							if(graph[edge].type == 2){
								auto t2 = target(edge, graph);
								if(graph[t2].code == C && !graph[t2].inAromaCycle)
									return true;
							}
							return false;
						});
						if(cnt == 2) // 2 double links both with C - 45 DC
						{ 
							dcs.emplace_back(*i, 45);
						}
						else // only one double link
						{
							dcs.emplace_back(*i, 46);
						}
					}
					else if (graph[*p].type == 3 && graph[tgt].code == C && !graph[tgt].inAromaCycle)
					{
						dcs.emplace_back(*i, 45);
					}
				}
			}
			auto range2 = make_pair(order2.begin(), order2.end());

			for (auto j = range2.first; j != range2.second; j++)
			{
				if (edges.second - edges.first != j->bonds.size())
					continue;
				if (!j->center.matches(graph[*i].code.code()))
					continue;
				//cout << "Matched CENTER " << j->center.symbol() << endl;
				if (j->valence != valency)
					continue;
				vector<bool> matched(j->bonds.size());
				vector<int> atoms; // atoms in this center
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
						//cout << "   matched " << j->bonds[k].atom.symbol()
						//	<< " vs " << graph[v].code.symbol() << endl;
						atoms.push_back(v);
						matched[k] = true;
						break;
					}
					if (k == j->bonds.size())
						break;
				}
				if (p == edges.second)
				{
					//cout << "DONE." << endl;
					dcs.emplace_back(*i, j->dc);
					dcsAtoms.insert(make_pair(*i, atoms));
				}
			}
		}
		sort(dcs.begin(), dcs.end(), [](const pair<vd, int>& a, const pair<vd, int>& b){
			return a.second < b.second;
		});
		cerr << "DCs:" << endl;
		for (auto& dc : dcs)
		{
			cerr << dc.first << " -DC-> " << dc.second << endl;
		}
		cerr << endl;
		//cout << endl;
	}

	template<class T>
	static const T& chainAt(const vector<T>& chain, int idx)
	{
		while (idx < 0)
			idx += (int)chain.size();
		while (idx >= (int)chain.size())
			idx -= (int)chain.size();
		return chain[idx];
	}

	static string keyatom(ChemGraph& g, vd v)
	{
		//numbers mean at least x links of given type
		//val == 0 - do not care
		struct Entry {
			char symbol;
			Code code;
			int valence;
			int singleLinks;
			int dualLinks;
			int trippleLinks;
		};
		Entry table[] = {
			{ 'M', Code(N), 0, 1, 1, 0 },
			{ 'N', Code(N), 0, 2, 0, 0 },
			{ 'Q', Code(O), 2, 0, 0, 0 },
			{ 'R', Code(O), 3, 0, 0, 0 },
			{ 'T', Code(S), 0, 1, 1, 0 },
			{ 'S', Code(S), 2, 2, 0, 0 }
		};
		char c = 0;
		for (auto& e : table)
		{
			if (g[v].code == e.code)
			{
				if (e.valence && g[v].valence != e.valence)
					continue;
				auto dt = multiCount(g, v);
				auto s = singleCount(g, v);
				if (e.singleLinks && e.singleLinks > s)
					continue; //not enough links
				if (e.dualLinks && e.dualLinks > dt.first)
					continue;
				if (e.trippleLinks && e.trippleLinks > dt.second)
					continue;
				return string(1, e.symbol);
			}
		}
		if (heteroatom(g[v]))
			return g[v].code.symbol();
		else
			return string();
	}

	static bool heteroatom(const Atom& atom)
	{
		return !atom.code.matches(C) && !atom.code.matches(H);
	}

	template<class T, class EdgeMap>
	static vector<int> cycleToChain(vector<T>& ic, EdgeMap&& mapper)
	{
		vector<int> vc;
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

	void locateCycles()
	{
		depth_first_search(graph, markLoops(cycles), get(&Atom::color, graph));
		for (auto & c : cycles)
		{
			transform(c.begin(), c.end(), c.begin(), [](const pair<int, int> &p){ 
				return p.first < p.second ? p : make_pair(p.second, p.first);
			});
			sort(c.begin(), c.end(), edge_less());
		}
		cerr << "BEFORE SPLIT" << endl;
		for_each(cycles.begin(), cycles.end(), &printCycle);
		vector<pair<int, int>> t, t2, uc, xc;
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
			pair<size_t, vector<pair<int, int>>*> result[4];
			result[0] = make_pair(len1, &c1);
			result[1] = make_pair(len2, &c2);			
			int m = 2;
			/*if (uc != c1)
				result[m++] = make_pair(lenU, &uc);
			else*/ if (xc != c2)
				result[m++] = make_pair(lenX, &xc);
			for (int k = 0; k < m; k++)
			{
				cerr << result[k].first << " ";
			}
			cerr << endl;
			sort(result, result + m, [](pair<size_t, vector<pair<int, int>>*> a, pair<size_t, vector<pair<int, int>>*> b){
				return a.first < b.first;
			});

			if ((&c1 == result[0].second && &c2 == result[1].second) ||
				(&c1 == result[1].second && &c2 == result[0].second))
				continue; //same cycles have won
			auto n1 = *result[1].second;
			auto n2 = *result[0].second;
			cerr << n1.size() << " " << n2.size() << endl;
			c1 = n1;
			c2 = n2;
		}
		cerr << "AFTER SPLIT" << endl;
		for_each(cycles.begin(), cycles.end(), &printCycle);
		for (auto ic : cycles)
		{
			assert(ic.size() > 2);
			auto vc = cycleToChain(ic, [](const pair<int, int>& p){ return p; });
			cerr << "CHAIN: ";
			for (int k : vc)
				cerr << k << " - ";
			cerr << endl;
			auto& g = graph;
			int piEl = 0;
			for_each(vc.begin(), vc.end(), [&g, &piEl](int n){
				cerr << "Atom # " << n << " " << g[n].code.symbol() << " pi E = " << g[n].piE << endl;
				if (g[n].piE > 0)
					piEl += g[n].piE;
				//TODO: add debug trace for < 0
			});
			bool aromatic = ((piEl - 2) % 4 == 0);
			cerr << "PI E " << piEl << " aromatic? : " << aromatic << endl;
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
	}

	void locateIrregular()
	{
		/*for (auto& r : irregularTempl)
		{
			vector<vector<pair<size_t, size_t>>> mappings;
			auto& g = graph;
			ullmann_all(r.piece, graph, [&r, &g](ChemGraph::vertex_descriptor a, ChemGraph::vertex_descriptor b){
				return r.piece[a].code == g[b].code;
			}, [&r, &g](ChemGraph::edge_descriptor a, ChemGraph::edge_descriptor b){
				return r.piece[a].type == g[b].type;
			}, mappings);
		}*/
	}

	void addDescriptorAtoms(vector<int>& fragment, vd dc)
	{
		if(dcsAtoms.find(dc) != dcsAtoms.end()){
			fragment.insert(fragment.end(), dcsAtoms[dc].begin(), dcsAtoms[dc].end());
		}
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
				TrackPath(graph, start, end, dcs), get(&Atom::color, graph));
			
			if (graph[dcs[i].first].path && graph[dcs[i].first].path < NON_PASSABLE)
			{
				bool coupled = false;
				auto &g = graph;
				vector<int> fragment;
				applyPath(start, end, [g, &coupled,&fragment](vd v){
					//cout << "Apply: " << g[v].code.symbol() << endl;
					if (g[v].code.matches(C) && g[v].piE > 0)
						coupled = true;
					fragment.push_back((int)v);
				});
				int dc1 = dcs[i].second; // end
				int dc2 = dcs[j].second; // start(!)
				int len = graph[dcs[i].first].path - 1;
				// SPECIAL case for 00 segment
				if(len == 0)
					coupled = true;
				// SPECIAL RULE for DC = 41 
				// then add 1 to length
				// Cannot generalize for all C not in aroma cycle 
				if (graph[start].code.matches(C) && dc2 == 41 && !graph[start].inAromaCycle)
					len += 1;
				if (graph[end].code.matches(C) && dc1 == 41 && !graph[end].inAromaCycle)
					len += 1;

				stringstream buffer;
				buffer << setfill('0') << setw(2) << dc1
					<< setfill('0') << setw(2) << len
					<< setfill('0') << setw(2) << dc2
					<< (coupled ? 1 : 0);
				addDescriptorAtoms(fragment, dcs[i].first);
				addDescriptorAtoms(fragment, dcs[j].first);
				outputJs(buffer.str(), fragment);
			}
		}
	}

	struct Edge{
		pair<int, int> e;
		int cycNum; //e принадлежит циклу cycles[cycNum]
		Edge(pair<int, int> e_, int cyc) :
			e(e_), cycNum(cyc){}
		bool operator<(const Edge& rhs)const
		{
			return edge_less()(e, rhs.e);
		}
	};

	//Строим запись "головы" двигаясь по огибающей (common) в обе стороны 
	string encodeHead(int firstCycle, vector<int>& commonCh, vector<Edge>& common, int totalCycles)
	{
		//Выбор опорного атома из стартового цикла
		auto& firstChain = chains[firstCycle];
		auto firstIdx = find_if(commonCh.begin(), commonCh.end(), [&firstChain](int v){
			return find(firstChain.begin(), firstChain.end(), v) != firstChain.end();
		});
		auto fIdx = firstIdx - commonCh.begin();
		// обход вперед/назад
		auto fwd = encodeCycle(1, commonCh, fIdx, firstCycle, common, totalCycles);
		auto bwd = encodeCycle(-1, commonCh, fIdx, firstCycle, common, totalCycles);
		return fwd < bwd ? fwd : bwd;
	}

	string encodeCycle(int dir, vector<int>& commonCh, int fIdx, int cycNum, vector<Edge>& common, int totalCycles)
	{
		stringstream cyclic_out;
		int edgeNum = 0; // on the most recent cycle
		int prevEdgeNum = 0; //edge count tracked on previous cycle
		size_t cycCount = 0; //first 2 are output as is, 3rd, 4th and so on need prefixes
		//for (auto& p : common)
		//	cout << p.e.first << ":" << p.e.second << endl;
		auto v1 = chainAt(commonCh, fIdx);
		for (int j = dir; ; j+= dir)
		{
			auto v2 = chainAt(commonCh, fIdx + j);
			auto edge = v1 < v2 ? make_pair(v1, v2) : make_pair(v2, v1);
			v1 = v2;
			auto p = equal_range(common.begin(), common.end(), Edge(edge, 0));
			assert(p.first != p.second);
			//cout << "Go with " << p.first->e.first << ":" << p.first->e.second << endl;
			edgeNum++;
			if (cycNum != p.first->cycNum)
			{
				//cout << "Cycle!" << endl;
				static string tab = "ABCDEFGHK";
				cycCount++;
				if (cycCount > 2)
				{
					assert(prevEdgeNum > 0 && prevEdgeNum < 9); //TODO: error on this
					cyclic_out << setw(1) << tab[prevEdgeNum-1]
						<< setw(1) << cycles[cycNum].size();
				}
				else
					cyclic_out << setw(1) << cycles[cycNum].size();
				cycNum = p.first->cycNum;
				prevEdgeNum = edgeNum;
				edgeNum = 0;
				if (cycCount == totalCycles) // all encoded
					break;
			}
		}
		return cyclic_out.str();
	}

	int pickNonKeyAtom(int dir, int fChain, int firstCycle, vector<int>& commonCh, vector<Edge>& common)
	{
		int res;
		for (int k = 0, j = fChain; k < (int)commonCh.size(); k++, j += dir)
		{
			auto v = chainAt(commonCh, j);
			auto v2 = chainAt(commonCh, j + dir);
			auto p = v < v2 ? make_pair(v, v2) : make_pair(v2, v);
			auto e = Edge(p, 0);
			auto r = equal_range(common.begin(), common.end(), e);
			assert(r.first != r.second);
			if (r.first->cycNum != firstCycle)
			{
				return j - dir;//step back 
			}
		}
		assert(false);
		return 0;
	}

	string encodeHeteroAtoms(int dir, int start, vector<int>& commonCh)
	{
		stringstream s;
		for (int i = 0, j = start; i < commonCh.size(); i++, j += dir)
		{
			auto v = chainAt(commonCh, j);
			string k = keyatom(graph, v);
			if (!k.empty())
			{
				s << k << i + 1;
			}
		}
		return s.str();
	}

	//Строим запись "хвоста" по стартовому циклу, двигаясь по огибающей (common) в обе стороны
	string encodeTail(int firstCycle, vector<int>& commonCh, vector<Edge>& common, int totalCycles)
	{
		auto& firstChain = chains[firstCycle];
		auto firstIdx = find_if(commonCh.begin(), commonCh.end(), [&firstChain](int v){
			return find(firstChain.begin(), firstChain.end(), v) != firstChain.end();
		});
		auto fi = firstIdx - commonCh.begin();
		//Идем в обе стороны и находим последние элементы в нашем цикле
		int firstV, secondV;
		//Правило: ключевые атомы (принадлежащие нескольким циклам) должны иметь наибольший номер
		firstV = pickNonKeyAtom(1, fi, firstCycle, commonCh, common);
		secondV = pickNonKeyAtom(-1, fi, firstCycle, commonCh, common);
		string variant[2];
		//firstV первый неключевой атом в "+" сторону, значит нумеруем в обратную
		variant[0] = encodeHeteroAtoms(-1, firstV, commonCh);
		//тоже для secondV
		variant[1] = encodeHeteroAtoms(1, secondV, commonCh);
		sort(variant, variant + 2);
		return variant[0];
	}

	//копирует ccv, тк будет не однакратно сортирован
	void encodeOnePolyCycle(vector<int> ccv)
	{
		//элементарные циклы кодируются отдельно
		if (ccv.size() == 1)
			return;
		vector<int> intercounts(cycles.size()); //кол-во пересечний с другими циклами
		for (auto n : ccv)
		{
			int c = 0;
			for (auto n2 : ccv)
			{
				if (intermap[n*cycles.size() + n2])
					c++;
			}
			intercounts[n] = c;
		}
		// cout << ccv.size() << endl;
		//Выбираем стартовый цикл - минимальный по числу пересечений (крайний), минимальный по числу элементов
		//ищем сразу 2 стартовых цикла, на случай, где они одинково хорошо подходят
		auto& cys = cycles;
		partial_sort(ccv.begin(), ccv.begin() + 2, ccv.end(), [&intercounts, &cys](int a, int b){
			int ia = intercounts[a];
			int ib = intercounts[b];
			return ia < ib || (ia == ib && cys[a].size() < cys[b].size());
		});
		auto start = ccv.begin();
		//Строим огибающий цикл (ребра)
		//в огибающей записывается номер цикла, которому принадлежит ребро
		vector<Edge> common;
		vector<size_t> idxs(ccv.size()); //idxs[i] --> cycles[ccv[i]][*]
		//алгоритм объединения N множест ребер
		//с исключением по предикату "принадлежит более чем одному множеству"
		for (;;)
		{
			//which edge to pick
			pair<int, int> m_edge(-1, -1);
			size_t smallestCCV = 0; //number of ccv entry having the smalles edge so far
			for (size_t i = 0; i < idxs.size(); i++)
			{
				auto n_edge = idxs[i];
				if (n_edge >= cys[ccv[i]].size())
					continue;
				auto edge = cys[ccv[i]][n_edge];
				if (m_edge.first < 0 || edge_less()(edge, m_edge))
				{
					m_edge = edge;
					smallestCCV = i;
				}
			}
			if (m_edge.first < 0)
				break;
			auto i = idxs[smallestCCV];
			auto& c = cys[ccv[smallestCCV]];
			idxs[smallestCCV]++;
			//проверка - c[i] должен принадлежать только одному эл. циклу
			int matches = 0;
			for (auto x : ccv)
			{
				auto er = equal_range(cys[x].begin(), cys[x].end(), c[i], edge_less());
				if (er.first != er.second)
					matches++;
			}
			if (matches == 1 && (common.empty() || common.back().e != c[i]))
				common.push_back(Edge(c[i], ccv[smallestCCV]));
			assert(is_sorted(common.begin(), common.end()));
		}
		auto commonCh = cycleToChain(common, [](const Edge& e){ return e.e; });
		string head = encodeHead(*start, commonCh, common, ccv.size());
		string head2 = head;
		// еще один стартовый цикл (если одинакового размера)
		if (cycles[*start].size() == cycles[*(start + 1)].size())
		{
			head2 = encodeHead(*(start + 1), commonCh, common, ccv.size());
		}
		head = head < head2 ? head : head2;
		//Суммируем Пи электроны по огибающей
		auto range = vertices(graph);
		auto sz = range.second - range.first;
		vector<bool> visited(sz);
		int piE = 0;
		for (int v : commonCh)
		{
			piE += graph[v].piE;
		}
		bool aromatic = (piE - 2) % 4 == 0;
		//Кодируем "хвост"
		//Выбираем новый цикл - минимальный по числу пересечений (крайний), _максимальный_ по числу элементов
		partial_sort(ccv.begin(), ccv.begin() + 2, ccv.end(), [&intercounts, &cys](int a, int b){
			int ia = intercounts[a];
			int ib = intercounts[b];
			return ia < ib || (ia == ib && cys[a].size() > cys[b].size());
		});
		start = ccv.begin();
		string tail = encodeTail(*start, commonCh, common, ccv.size());
		string tail2 = tail;
		// еще один стартовый цикл (если одинакового размера)
		if (cycles[*start].size() == cycles[*(start + 1)].size())
			tail2 = encodeTail(*(start + 1), commonCh, common, ccv.size());
		if (tail > tail2)
			tail = tail2;
		stringstream buffer;
		vector<int> fragment;

		buffer << head;
		buffer << "," << setfill('0') << setw(2) << (aromatic ? piE : 0);
		buffer << tail;
		for(auto cc : ccv){
			fragment.insert(fragment.end(), chains[cc].begin(), chains[cc].end());
		}
		outputJs(buffer.str(), fragment);
	}

	void recursePoly(vector<int>& poly, map<int, vector<int>>& adj_list, set<vector<int>>& used)
	{
		auto e = poly.back();
		for (auto x : adj_list[e])
		{
			{
				if (find(poly.begin(), poly.end(), x) != poly.end())
					continue;
				assert(intermap[x + e*cycles.size()]);
				poly.push_back(x);
				if (poly.size() > 1)
				{
					auto cp = poly;
					sort(cp.begin(), cp.end());
					if (used.find(cp) == used.end())
					{
						encodeOnePolyCycle(cp);
						using std::move;
						used.insert(move(cp));
					}
				}
				//добавить еще один цикл
				if (poly.size() < adj_list.size())
					recursePoly(poly, adj_list, used);
				poly.pop_back();
			}
		}
	}

	void encodePolyCycles(vector<int>& ccv)
	{
		map<int, vector<int>> adj_list; //ccv[i] --> adj_list[i]
		for (auto x : ccv)
		{
			vector<int> neib;
			for (auto y : ccv)
			{
				if (y != x && intermap[y*cycles.size() + x])
				{
					neib.push_back(y);
				}
			}
			using std::move;
			adj_list.insert(make_pair(x, move(neib)));
		}
		vector<int> poly;
		set<vector<int>> used;
		for (auto x : ccv)
		{
			poly.push_back(x);
			recursePoly(poly, adj_list, used);
			poly.pop_back();
			assert(poly.size() == 0);
		}
	}

	void cyclic(ostream& out)
	{
		// Кодирование простых циклов
		// Оставить только ароматические?
		for (auto c : chains)
		{
			int piE = 0;
			if (graph[c.front()].inAromaCycle)
			{
				for (int v : c)
					piE += graph[v].piE;
			}
			vector<pair<string, int>> hatoms;
			for (int v : c)
			{
				string s = keyatom(graph, v);
				if (!s.empty())
				{
					hatoms.emplace_back(s, v);
				}
			}

			auto pivot = min_element(hatoms.begin(), hatoms.end(), [](const pair<string, int> & lhs, const pair<string, int> &rhs){
				return lhs.first < rhs.first;
			});
			stringstream buffer;
			vector<int> fragment;
			fragment.insert(fragment.end(), c.begin(), c.end());
			buffer << setfill('0') << setw(1) << c.size()
				<< ","
				<< setfill('0') << setw(2) << piE;
			bool heterocycle = pivot != hatoms.end();
			if (heterocycle)
			{
				string left_descr, right_descr;
				auto idx = pivot - hatoms.begin();
				for (int k = 0; k < (int)hatoms.size(); k++)
				{
					auto& nextL = chainAt(hatoms, idx + k);
					auto& nextR = chainAt(hatoms, idx - k);
					auto number = lexical_cast<string>(k + 1);
					left_descr += nextL.first + number;
					right_descr += nextR.first + number;
				}
				buffer << (left_descr < right_descr ? left_descr : right_descr);
			}
			outputJs(buffer.str(), fragment);
		}
		//make a map of intersections
		intermap.resize(cycles.size()*cycles.size());
		intermap.clear();
		vector<pair<int, int>> t;
		for (size_t i = 0; i < cycles.size(); i++)
		for (size_t j = i + 1; j < cycles.size(); j++)
		{
			auto& c1 = cycles[i];
			auto& c2 = cycles[j];
			t.resize(min(c1.size(), c2.size()));
			auto tend = set_intersection(c1.begin(), c1.end(), c2.begin(), c2.end(), t.begin(), edge_less());
			bool intersection = t.begin() != tend;
			intermap[i*cycles.size() + j] = intersection;
			intermap[j*cycles.size() + i] = intersection;
		}
		/*for (size_t i = 0; i < cycles.size(); i++)
		{
			for (size_t j = 0; j < cycles.size(); j++)
			{
				cout << intermap[i*cycles.size() + j];
			}
			cout << endl;
		}*/
		//
		vector<int> ccv; //chain - indices of cycles
		vector<bool> used(cycles.size());
		do
		{
			ccv.clear();
			for (size_t i = 0; i < cycles.size();)
			{
				if (used[i])
				{
					i++;
					continue;
				}
				if (ccv.empty())
				{
					used[i] = true;
					ccv.push_back(i);
					i++; // первый непроверенный, нет нужды начинать с 0
					continue;
				}
				//then must interesect some other cycle in this chain
				bool found = false;
				for (auto j : ccv)
				{
					if (intermap[i*cycles.size() + j])
					{
						found = true;
						break;
					}
				}
				if (found)
				{
					used[i] = true;
					ccv.push_back(i);
					i = 0; //нужно перепроверить циклы
					continue;
				}
				i++;
			}
			if (ccv.size() > 1) //элементарные циклы уже учтены
			{
				encodePolyCycles(ccv);
			}
		} while (!ccv.empty());
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
			ullmann_all(r.piece, graph, [&](ChemGraph::vertex_descriptor a, ChemGraph::vertex_descriptor b){				
				if(!r.piece[a].code.matches(g[b].code.code()))
					return false;
				if(a == r.a1 || a == r.a2){
					if(g[b].inAromaCycle) // none of replacemnt dc are in aroma cycle (e.g. DC 41 is CH3)
						return false;
					return find_if(dcs.begin(), dcs.end(), [&](pair<vd,size_t> dcp){
						return dcp.first == b;
					}) != dcs.end();
				}
				else
					return true;
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
				stringstream buffer;
				vector<int> fragment;
				for(auto& p : m)
					fragment.push_back(p.first);
				addDescriptorAtoms(fragment, fV);
				addDescriptorAtoms(fragment, sV);

				buffer << setfill('0') << setw(2) << fdc
					<< setfill('0') << setw(2) << r.dc
					<< setfill('0') << setw(2) << sdc
					<< r.coupling;
				outputJs(buffer.str(), fragment);
			}
		}
	}

private:
	std::vector<LevelOne> order1;
	std::vector<LevelTwo> order2;
	std::vector<Replacement> repls;
	ChemGraph graph;
	//location of DCs in 'graph' and their numeric value
	vector<pair<vd, int>> dcs;
	map<vd, vector<int>> dcsAtoms; // extra atoms that belong to each DC  
	//
	vector<string> jsons;
	//sorted arrays of edges - basic cycles
	vector<vector<pair<int, int>>> cycles;
	//same basic cycles represented as chains of vertices
	vector<vector<int>> chains;
	//map of intersection between cycles
	vector<bool> intermap;
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
