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
#include <boost/graph/vf2_sub_graph_iso.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/visitors.hpp>

#include "conv.hpp"
#include "ctab.hpp"
#include "descriptors.hpp"
#include "fcsp.hpp"
#include "log.hpp"

enum { NON_PASSABLE = 10000 };
using namespace std;
using namespace boost;

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

bool is_exclusive_dc(int dc)
{
	return dc == 45 || dc == 46; // these are impassable if one is part of chain
}

// Select 1-s from bool matrix so that 
// there is at least one selected per row
// and no more then one per column
// Input:
// matrix - bool matrix
// selections : row idx --> col idx
// cur - number of currently matching row
bool pickChoices(const vector<bool>& matrix, size_t rows, size_t cols, vector<size_t>& selections, size_t cur=0)
{
	if(cur == rows)
		return true;
	for(size_t i=0; i<cols; i++)
	{
		if(matrix[cur*cols + i])
		{
			selections[cur] = i;
			auto end = selections.begin()+cur;
			// not chosen before
			if(find(selections.begin(), selections.begin()+cur, i) == end)
			{
				if(pickChoices(matrix, rows, cols, selections, cur+1))
					return true;
				else
					selections[cur] = -1; // reset and continue
			}
		}
	}
	return false;
}

struct TrackPath: public default_bfs_visitor {
	ChemGraph& g;
	vector<pair<vd, int>>& dcs;
	ChemGraph::vertex_descriptor start, tgt;
	bool pass_4546;
	TrackPath(ChemGraph& graph, ChemGraph::vertex_descriptor s, ChemGraph::vertex_descriptor t, 
		vector<pair<vd, int>>& dcsArr) :
		g(graph), dcs(dcsArr), start(s), tgt(t)
	{
		auto sdc = *find_if(dcs.begin(), dcs.end(), [&](pair<vd,int> p){
			return p.first == start;
		});
		auto edc = *find_if(dcs.begin(), dcs.end(), [&](pair<vd,int> p){
			return p.first == tgt;
		});
		pass_4546 = sdc.second != 45 && sdc.second != 46 && edc.second != 45 && edc.second != 46;
		LOG(TRACE) << (pass_4546 ? "Passing" :"Not passing") <<" DC 45-46 while going "
			<< sdc.second << "-->"<< edc.second << endline; 
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
		//cout << s << "-->" << d << endline;
		if (g[d].inAromaCycle && d != tgt) //all paths  of aromatic cycle  are inpassable
			g[d].path = NON_PASSABLE;
		else if (g[d].inAromaCycle && g[s].inAromaCycle)
			g[d].path = NON_PASSABLE;
		else if (g[d].code != C && d != tgt)
			g[d].path = NON_PASSABLE;
		else if(d != tgt && !pass_4546 && g[d].code == C)
		{
			bool hit4546 = find_if(dcs.begin(), dcs.end(), [d](pair<vd,int> p){
				return p.first == d && (p.second == 45 || p.second == 46);
			}) != dcs.end();

			if(hit4546)
			{
				auto sdc = *find_if(dcs.begin(), dcs.end(), [&](pair<vd,int> p){
					return p.first == start;
				});
				auto edc = *find_if(dcs.begin(), dcs.end(), [&](pair<vd,int> p){
					return p.first == tgt;
				});
				LOG(DEBUG) << "Hit non-passable DC 45/46 while going " 
					<< sdc.second << "-->" << edc.second << endline;
				g[d].path = NON_PASSABLE;
			}
			else
				g[d].path = g[s].path + 1;
		}
		else
			g[d].path = g[s].path + 1;
	}
};

// callback for VF2 algorithm
struct CollectAsVectors{
	ChemGraph& pattern;
	vector<vector<size_t>>& mappings; // mapped atoms in the bigger graph

	
    CollectAsVectors(ChemGraph& pat, vector<vector<size_t>>& _mappings):pattern(pat), mappings(_mappings){}
	template <typename CorrespondenceMap1To2,
          typename CorrespondenceMap2To1>
	bool operator()(CorrespondenceMap1To2 f, CorrespondenceMap2To1 g) const
	{
		auto vtx = vertices(pattern);
		vector<size_t> in_big;
		for(auto p=vtx.first; p != vtx.second; p++)
		{
			in_big.push_back(get(f, *p));
		}
		using std::move;
		mappings.push_back(move(in_big));
		return true;
	}
};

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
	Impl(FCSPOptions opts) :
		order1(std::move(opts.first)), order2(std::move(opts.second)), 
		repls(std::move(opts.replacements)),
		long41(opts.long41), format(opts.format){}

	void load(istream& inp)
	{
		CTab tab = readMol(inp);
		graph = toGraph(tab);
		auto comp = connectedComponents(graph);
	}
	// Очистить все переменные состояния кодировщика
	void clear()
	{
		outPieces.clear();
		outPiecesCycle.clear();
		dcs.clear();
		dcsAtoms.clear();
		reserved_dcs.clear();
		cycles.clear();
	}

	void sortDCs()
	{
		// filter out things that got reserved
		auto before = dcs.size();
		dcs.erase(remove_if(dcs.begin(), dcs.end(), [&](const pair<vd, int>& a){
			auto p = reserved_dcs.find(a.first);
			// LOG(FATAL) << "> " << (p != reserved_dcs.end()) << endline;
			return p != reserved_dcs.end() && p->second != a.second;
		}), dcs.end());
		LOG(INFO) << "Filtered " << before - dcs.size() << " DCs in favor of monolithic patterns."<< endline;
		sort(dcs.begin(), dcs.end(), [](const pair<vd, int>& a, const pair<vd, int>& b){
			return a.second < b.second;
		});
		LOG(INFO) << "Found DCs:" << endline;
		for (auto& dc : dcs)
		{
			LOG(INFO) << dc.first << " -DC-> " << dc.second << endline;
		}
		LOG(INFO) << endline;
	}

	void process(ostream& out, string filename)
	{
		clear(); // clear state
		addHydrogen(graph);
		locatePiElectrons();
		locateCycles(); //adds cyclic DCs
		locateDCs(false);
		locateIrregular();
		sortDCs();
		cyclic(out);
		linear(out);
		locateDCs(true); // REPL-only DCs - 44 so far
		sortDCs();
		replacement(out);
		outputWhole(out, filename);
	}

	void dumpGraph(ostream& out)
	{
		::dumpGraph(graph, out);
	}

	void outputPieceCycle(string code, vector<int>& atoms)
	{
		outputPiece_(code, atoms, true);
	}

	void outputPiece(string code, vector<int>& atoms)
	{
		outputPiece_(code, atoms, false);
	}

	void outputPiece_(string code, vector<int>& atoms, bool cycle)
	{
		if(format == FCSPFMT::JSON)
		{
			stringstream s;
			s << "{" << "\"code\" : \""<<code<<"\",";
			s << "\"place\": ";
			printJsArray(atoms, s); 
			s << "}";
			if(cycle)
				outPiecesCycle.push_back(s.str());
			else
				outPieces.push_back(s.str());
		}
		else if(format == FCSPFMT::CSV || format == FCSPFMT::TXT)
		{
			if(cycle)
				outPiecesCycle.push_back(code);
			else
				outPieces.push_back(code);
		}
	}

	void outputWhole(ostream& out, string filename)
	{
		sort(outPiecesCycle.begin(), outPiecesCycle.end());
		sort(outPieces.begin(), outPieces.end());
		outPieces.insert(outPieces.begin(), outPiecesCycle.begin(), outPiecesCycle.end());
		if(format == FCSPFMT::JSON)
		{
			out << "[";
			for(size_t i=0; i<outPieces.size();i++)
			{
				if(i != 0)
					out << ", ";
				out << outPieces[i];
			}
			out << "]";
			out << endline;
		}
		else if(format == FCSPFMT::CSV || format == FCSPFMT::TXT)
		{
			if(format == FCSPFMT::CSV)
				out << filename << ';';
			for(size_t i=0; i<outPieces.size();i++)
			{
				if(i != 0)
					out << " ";
				out << outPieces[i];
			}
			out << endline;
		}
	}

	void locatePiElectrons()
	{
		auto vrtx = vertices(graph);
		for (auto i = vrtx.first; i != vrtx.second; i++)
		{
			//pre-calculate per-atom properties
			int valency = getValence(graph, *i);
			graph[*i].valence = valency;
			auto dual_tripple = multiCount(graph, *i);
			int piE = countPiElectrons(graph[*i].code, valency, dual_tripple.first, dual_tripple.second);
			graph[*i].piE = piE > 0 ? piE : 0;
		}
	}

	void locateDCs(bool replOnly)
	{
		auto vrtx = vertices(graph);
		for (auto i = vrtx.first; i != vrtx.second; i++)
		{
			int valency = graph[*i].valence;
			auto edges = out_edges(*i, graph);
			LevelOne t{graph[*i].code, valency, 0};
			auto range = equal_range(order1.begin(), order1.end(), t);
			if(!replOnly) // skip level-1 DCs and 45-46 for repl-only DCs
			{
				for (auto j = range.first; j != range.second; ++j)
				{
					dcs.emplace_back(*i, j->dc);
				}
				if(graph[*i].code == C && !graph[*i].inAromaCycle)
				{
					// check for 45 & 46
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
			}
			// select a range of DCs pattern with center matching this atom
			auto range2 = make_pair(order2.begin(), order2.end());
			for (auto j = range2.first; j != range2.second; j++)
			{
				if(j->replOnly != replOnly) // can't use during this stage
					continue;
				if (edges.second - edges.first < (int)j->bonds.size())
					continue;
				if (!j->center.matches(graph[*i].code))
					continue;
				if (j->valence != valency)
					continue;
				LOG(TRACE) << "Candidate DC "<< j->dc <<" CENTER " << j->center.symbol() << " VALENCE "<< valency << endline;

				vector<int> atoms; // atoms in this center
				size_t cand_bnds = edges.second - edges.first;
				size_t smpl_bnds = j->bonds.size();
				// Put ones for combinations that match. A row per edge in a DC pattern (sample).
				// Then we need to pick one in each row, if at least one row is all zeros - no match
				// Same DC may happen twice in the same atom, but we don't accomodate for that (for now)
				vector<bool> mappings(cand_bnds*smpl_bnds); 
				LOG(TRACE) << "MAPPING:" <<endline;
				LOG(TRACE) << "  ";
				for(size_t k=0; k<cand_bnds; k++)
				{
					auto e = (edges.first+k);
					auto t = target(*e, graph);
					LOG(TRACE) << setw(2) << graph[t].code.symbol();
				}
				for(size_t p=0; p<smpl_bnds; p++)
				{
					LOG(TRACE) << endline << setw(2) << j->bonds[p].atom.symbol();
					for(size_t q=0; q<cand_bnds; q++)
					{
						auto e = edges.first + q;
						if(graph[*e].type == j->bonds[p].bondType)
						{
							auto t = target(*e, graph);
							if(j->bonds[p].atom.matches(graph[t].code))
								mappings[p*cand_bnds + q] = true;
						}
						LOG(TRACE) << setw(2) << mappings[p*cand_bnds + q];
					}
				}
				LOG(TRACE) << endline;
				vector<size_t> found_mapping(smpl_bnds); // sample idx --> candidate idx
				if(pickChoices(mappings, smpl_bnds, cand_bnds, found_mapping))
				{
					for(auto idx : found_mapping)
					{
						// get atom by index of edge around this atom
						auto t = target(*(edges.first + idx), graph);
						atoms.push_back(t);
					}
					dcs.emplace_back(*i, j->dc);
					dcsAtoms.insert(make_pair(*i, atoms));
					// only assign DCs to reserve once
					if(j->monolith && reserved_dcs.find(*i) == reserved_dcs.end())
					{ 
						LOG(DEBUG) << "Reserved for DC "<< j->dc << " "<< atoms.size() + 1 <<" atoms"<<endline;
						//reserve atoms that belong to this DC
						reserved_dcs[*i] = j->dc;
						for(auto v : atoms)
						{
							if(graph[v].code != C) //FIXME: should be more sensible
								reserved_dcs[v] = j->dc;
						}
					}
				}
			}
		}
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
			{ 'M', N, 0, 1, 1, 0 },
			{ 'N', N, 0, 2, 0, 0 },
			{ 'Q', O, 2, 0, 0, 0 },
			{ 'R', O, 3, 0, 0, 0 },
			{ 'T', S, 0, 1, 1, 0 },
			{ 'S', S, 2, 2, 0, 0 }
		};
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

	static bool heteroatom(const AtomVertex& atom)
	{
		return !atom.code.matches(C) && !atom.code.matches(H);
	}

	
	void locateCycles()
	{
		cycles = minimalCycleBasis(graph);
		for (auto& ic : cycles)
		{
			auto& vc = ic.chain;
			ic.markAromatic(graph);
			if(ic.aromatic()) // aromatic specific DCs
			{
				for(auto n : vc)
				{
					//any atom in aromatic cycle
					dcs.emplace_back(n, 33);
					//hetero-atom in aromatic cycle - DC 34
					if (!graph[n].code.matches(C))
						dcs.emplace_back(n, 34);
					auto idx = find(vc.begin(), vc.end(), n) - vc.begin();
					assert(idx != (int)vc.size());
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
					//TODO: DC 40 - need to handle charges carefully to test charged atoms
				}
			}
		}
	}

	void locateIrregular() //  TODO: rewrite
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
			//cout << graph[current].path << " -->" << endline;
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
		// DCs are sorted as a[0] < .. < a[n-1]
		for (size_t i = 0; i < dcs.size(); i++)
		for (size_t j = i  + 1; j < dcs.size(); j++)
		{
			vd start = dcs[i].first;
			vd end = dcs[j].first;
			int start_dc = dcs[i].second;
			int end_dc = dcs[j].second;
			breadth_first_search(graph, start, buf,
				TrackPath(graph, start, end, dcs), get(&AtomVertex::color, graph));
			if (graph[end].path && graph[end].path < NON_PASSABLE)
			{
				bool coupled = true; //0-length path is therefore coupled (FIXME: check PI el-s too)
				auto &g = graph;
				vector<int> fragment;
				if(start_dc == 41) // check only the first 
				{
					bool check = true;
					applyPath(start, end, [g, end, &check, &coupled, &fragment](vd v){
						if(check)
						{
							if (v != end && g[v].code.matches(C) && g[v].piE == 0)
								coupled = false;
							check = false;
						}
						fragment.push_back((int)v);
					});
				}
				else
				{
					applyPath(start, end, [g, end, &coupled, &fragment](vd v){
						if (v != end && g[v].code.matches(C) && g[v].piE == 0)
							coupled = false;
						fragment.push_back((int)v);
					});
				}
				int len = graph[end].path - 1;
				//SPECIAL CASE - TODO verifiy correctness
				if(len == 0)
					coupled = true;
				if (len == 0 && is_exclusive_dc(start_dc) && start_dc == end_dc)
				{
					if(g[edge(start, end, g).first].type > 1) // connected with non-single link
					{
						continue; // skip output - this is self-referental
					}
				}
				// SPECIAL CASE - "The long 41" rule
				if(start_dc == 41 && long41)
					len += 1;
				if(end_dc == 41 && long41)
					len += 1;
				stringstream buffer;
				buffer << setfill('0') << setw(2) << start_dc
					<< setfill('0') << setw(2) << len
					<< setfill('0') << setw(2) << end_dc
					<< (coupled ? 1 : 0);
				addDescriptorAtoms(fragment, dcs[i].first);
				addDescriptorAtoms(fragment, dcs[j].first);
				outputPiece(buffer.str(), fragment);
			}
		}
	}

	struct Edge{
		pair<vd, vd> e;
		int cycNum; //e принадлежит циклу cycles[cycNum]
		Edge(pair<vd, vd> e_, int cyc) :
			e(e_), cycNum(cyc){}
		bool operator<(const Edge& rhs)const
		{
			return e < rhs.e;
		}
	};

	//Строим запись "головы" двигаясь по огибающей (common) в обе стороны 
	string encodeHead(int firstCycle, vector<vd>& commonCh, vector<Edge>& common, int totalCycles)
	{
		//Выбор опорного атома из стартового цикла
		auto& firstChain = cycles[firstCycle].chain;
		auto firstIdx = find_if(commonCh.begin(), commonCh.end(), [&firstChain](int v){
			return find(firstChain.begin(), firstChain.end(), v) != firstChain.end();
		});
		auto fIdx = firstIdx - commonCh.begin();
		// обход вперед/назад
		auto fwd = encodeCycle(1, commonCh, fIdx, firstCycle, common, totalCycles);
		auto bwd = encodeCycle(-1, commonCh, fIdx, firstCycle, common, totalCycles);
		return fwd < bwd ? fwd : bwd;
	}

	string encodeCycle(int dir, vector<vd>& commonCh, int fIdx, int cycNum, vector<Edge>& common, int totalCycles)
	{
		stringstream cyclic_out;
		int edgeNum = 0; // on the most recent cycle
		int prevEdgeNum = 0; //edge count tracked on previous cycle
		size_t cycCount = 0; //first 2 are output as is, 3rd, 4th and so on need prefixes
		//for (auto& p : common)
		//	cout << p.e.first << ":" << p.e.second << endline;
		auto v1 = chainAt(commonCh, fIdx);
		for (int j = dir; ; j+= dir)
		{
			auto v2 = chainAt(commonCh, fIdx + j);
			auto edge = v1 < v2 ? make_pair(v1, v2) : make_pair(v2, v1);
			v1 = v2;
			auto p = equal_range(common.begin(), common.end(), Edge(edge, 0));
			assert(p.first != p.second);
			//cout << "Go with " << p.first->e.first << ":" << p.first->e.second << endline;
			edgeNum++;
			if (cycNum != p.first->cycNum)
			{
				//cout << "Cycle!" << endline;
				static string tab = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
				cycCount++;
				if (cycCount > 2)
				{
					assert(prevEdgeNum > 0 && prevEdgeNum < 20); //TODO: error on this
					cyclic_out << setw(1) << tab[prevEdgeNum-1]
						<< setw(1) << cycles[cycNum].edges.size();
				}
				else
					cyclic_out << setw(1) << cycles[cycNum].edges.size();
				cycNum = p.first->cycNum;
				prevEdgeNum = edgeNum;
				edgeNum = 0;
				if (cycCount == totalCycles) // all encoded
					break;
			}
		}
		return cyclic_out.str();
	}

	int pickNonKeyAtom(int dir, int fChain, int firstCycle, vector<vd>& commonCh, vector<Edge>& common)
	{
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

	string encodeHeteroAtoms(int dir, int start, vector<vd>& commonCh)
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
	string encodeTail(int firstCycle, vector<vd>& commonCh, vector<Edge>& common, int totalCycles)
	{
		auto& firstChain = cycles[firstCycle].chain;
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
		//Выбираем стартовый цикл - минимальный по числу пересечений (крайний), минимальный по числу элементов
		//ищем сразу 2 стартовых цикла, на случай, где они одинково хорошо подходят
		auto& cys = cycles;
		partial_sort(ccv.begin(), ccv.begin() + 2, ccv.end(), [&intercounts, &cys](int a, int b){
			int ia = intercounts[a];
			int ib = intercounts[b];
			return ia < ib || (ia == ib && cys[a].edges.size() < cys[b].edges.size());
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
			pair<vd, vd> m_edge(-1, -1);
			size_t smallestCCV = 0; //number of ccv entry having the smalles edge so far
			for (size_t i = 0; i < idxs.size(); i++)
			{
				auto n_edge = idxs[i];
				if (n_edge >= cys[ccv[i]].edges.size())
					continue;
				auto edge = cys[ccv[i]].edges[n_edge];
				if (m_edge.first == (size_t)-1 || edge < m_edge)
				{
					m_edge = edge;
					smallestCCV = i;
				}
			}
			if (m_edge.first == (size_t)-1)
				break;
			auto i = idxs[smallestCCV];
			auto& c = cys[ccv[smallestCCV]];
			idxs[smallestCCV]++;
			//проверка - c[i] должен принадлежать только одному эл. циклу
			int matches = 0;
			for (auto x : ccv)
			{
				auto er = equal_range(cys[x].edges.begin(), cys[x].edges.end(), c.edges[i]);
				if (er.first != er.second)
					matches++;
			}
			if (matches == 1 && (common.empty() || common.back().e != c.edges[i]))
				common.push_back(Edge(c.edges[i], ccv[smallestCCV]));
			assert(is_sorted(common.begin(), common.end()));
		}
		LOG(DEBUG) << "COMMON EDGE SET:";
		for(auto e : common)
			LOG(DEBUG) << e.e << " ";
		LOG(DEBUG) << endline;
		auto commonCh = cycleToChain(common, [](const Edge& e){ return e.e; });
		string head = encodeHead(*start, commonCh, common, ccv.size());
		string head2 = head;
		// еще один стартовый цикл (если одинакового размера)
		if (cycles[*start].edges.size() == cycles[*(start + 1)].edges.size())
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
			return ia < ib || (ia == ib && cys[a].edges.size() > cys[b].edges.size());
		});
		start = ccv.begin();
		string tail = encodeTail(*start, commonCh, common, ccv.size());
		string tail2 = tail;
		// еще один стартовый цикл (если одинакового размера)
		if (cycles[*start].edges.size() == cycles[*(start + 1)].edges.size())
			tail2 = encodeTail(*(start + 1), commonCh, common, ccv.size());
		if (tail > tail2)
			tail = tail2;
		stringstream buffer;
		vector<int> fragment;

		buffer << head;
		buffer << "," << setfill('0') << setw(2) << (aromatic ? piE : 0);
		buffer << tail;
		for(auto cc : ccv){
			fragment.insert(fragment.end(), cycles[cc].chain.begin(), cycles[cc].chain.end());
		}
		outputPieceCycle(buffer.str(), fragment);
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
		for (auto cyc : cycles)
		{
			auto& ch = cyc.chain;
			int piE = 0;
			// NEW RULE - always output piE for singleton cycles and coupling linked systems
			for (vd v : ch)
				piE += graph[v].piE;
			vector<pair<string, vd>> hatoms;
			for (vd v : ch)
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
			fragment.insert(fragment.end(), ch.begin(), ch.end());
			buffer << setfill('0') << setw(1) << cyc.edges.size()
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
					auto number = to<string>(k + 1);
					left_descr += nextL.first + number;
					right_descr += nextR.first + number;
				}
				buffer << (left_descr < right_descr ? left_descr : right_descr);
			}
			outputPieceCycle(buffer.str(), fragment);
		}
		//make a map of intersections
		intermap.resize(cycles.size()*cycles.size());
		intermap.clear();
		for (size_t i = 0; i < cycles.size(); i++)
		for (size_t j = i + 1; j < cycles.size(); j++)
		{
			auto& c1 = cycles[i];
			auto& c2 = cycles[j];
			bool intersection = c1.intersects(c2);
			intermap[i*cycles.size() + j] = intersection;
			intermap[j*cycles.size() + i] = intersection;
		}/*
		for (size_t i = 0; i < cycles.size(); i++)
		{
			for (size_t j = 0; j < cycles.size(); j++)
			{
				cout << intermap[i*cycles.size() + j];
			}
			cout << endline;
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

	// return all descriptors at dcV vertex
	vector<int> mapDCs(size_t dcV)
	{
		vector<int> dc_nums;
		for(auto p : dcs){
			if(p.first == dcV)
				dc_nums.push_back(p.second);
		}
		return dc_nums;
	}

	void replacement(ostream& out)
	{
		//cout << "REPLACEMENTS!" << endline;
		for (auto& r : repls)
		{
			vector<vector<size_t>> mappings;
			auto& g = graph;
			
			vf2_subgraph_iso(r.piece, g, CollectAsVectors(g, mappings), vertex_order_by_mult(r.piece),
				vertices_equivalent([&](ChemGraph::vertex_descriptor a, ChemGraph::vertex_descriptor b){				
					if(!r.piece[a].code.matches(g[b].code))
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
				}).edges_equivalent(
				[&r, &g](ChemGraph::edge_descriptor a, ChemGraph::edge_descriptor b){
					return r.piece[a].type == g[b].type;
				})
			);
			vector<pair<pair<int, int>, pair<int, int>>> used_pairs;
			for (auto& m : mappings)
			{
				int fV = m[r.a1];
				int sV = m[r.a2];				
				// check if this pair of vertex-DC pairs was used before
				auto fdc = mapDCs(fV);
				if (fdc.size() == 0)
					continue;
				auto sdc = mapDCs(sV);
				if (sdc.size() == 0)
					continue;
				for(auto f : fdc)
				for(auto s : sdc)
				{
					auto firstDC = f;
					auto secondDC = s;
					auto firstV = fV;
					auto secondV = sV;
					// the usual rule of smaller DC first
					if(firstDC > secondDC)
					{
						swap(firstDC, secondDC);
						swap(firstV, secondV);
					}
					auto check_val = make_pair(make_pair(firstV, firstDC), make_pair(secondV,secondDC));
					if (find(used_pairs.begin(), used_pairs.end(),check_val) != used_pairs.end())
							continue;
					used_pairs.emplace_back(make_pair(firstV, firstDC), make_pair(secondV,secondDC));
					
					stringstream buffer;
					vector<int> fragment;
					for(auto& p : m)
						fragment.push_back(p);
					addDescriptorAtoms(fragment, fV);
					addDescriptorAtoms(fragment, sV);

					buffer << setfill('0') << setw(2) << firstDC
						<< setfill('0') << setw(2) << r.dc
						<< setfill('0') << setw(2) << secondDC
						<< r.coupling;
					outputPiece(buffer.str(), fragment);
				}
			}
		}
	}

private:
	std::vector<LevelOne> order1;		// patterns for first-order DCs
	std::vector<LevelTwo> order2;		// patterns for second-order DCs
	std::vector<Replacement> repls; // patterns for replacement decsriptors (not DCs)
	bool long41;										// if true - DC #41 adds +1 to the length of chain
	FCSPFMT format;				 					// controls output format
	ChemGraph graph;								// mol graph
	//location of DCs in 'graph' and their numeric value
	vector<pair<vd, int>> dcs;			// sorted by vertex array of vertex->dc mappings
	map<vd, vector<int>> dcsAtoms;  // extra atoms that belong to each DC  
	std::unordered_map<vd, int> reserved_dcs; // atoms reserved by specific monolithic DC
	// unassembled output chunks, assembly depends on format variable
	vector<string> outPiecesCycle; // cycle descriptors go first on assembly
	vector<string> outPieces;
	//sorted arrays of edges - basic cycles
	//vector<vector<pair<int, int>>> cycles;
	//same basic cycles represented as chains of vertices
	//vector<vector<int>> chains;
	vector<Cycle> cycles;
	//map of intersection between cycles
	vector<bool> intermap;
};

FCSP::FCSP(FCSPOptions opts) :
	pimpl(new FCSP::Impl(std::move(opts))){}

void FCSP::load(std::istream& inp)
{
	pimpl->load(inp);
}


void FCSP::dumpGraph(std::ostream& dot)
{
	pimpl->dumpGraph(dot);
}

void FCSP::process(std::ostream& out, string filename)
{
	pimpl->process(out, filename);
}

FCSP::~FCSP(){}
