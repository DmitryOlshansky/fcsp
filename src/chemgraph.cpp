// Std
#include <algorithm>
#include <deque>
#include <unordered_map>
// Boost
#include <boost/graph/graphviz.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/visitors.hpp>
// Our stuff 
#include "chemgraph.hpp"
#include "ctab.hpp"
#include "log.hpp"
#include "gauss.hpp"

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

ChemGraph& addHydrogen(ChemGraph& graph)
{
	std::unordered_map<int, int> valences;
	valences.insert(make_pair(C.code(), 4));
	valences.insert(make_pair(N.code(), 3));
	valences.insert(make_pair(O.code(), 2));
	auto vtx = vertices(graph);
	for(auto p = vtx.first; p != vtx.second; p++)
	{
		//Note: no extra hydrogens for negative ions
		if(valences.find(graph[*p].code.code()) != valences.end()
			&& graph[*p].code.charge() >= 0)
		{
			int normal_valence = valences[graph[*p].code.code()];
			// FIXME: counts 1.5 as 4 but that is "works for me"
			int cur_val = getValence(graph, *p); 
			// fit with hydrogens if not enough valence
			for(int i=cur_val; i<normal_valence; i++) 
			{
				size_t v = add_vertex(AtomVertex(H), graph);
				add_edge(*p, v, Bound(1), graph);
			}
		}
	}
	return graph;
}


class CodeWriter {
public:
	CodeWriter(ChemGraph& graph):g(graph){}
	template <class VertexOrEdge>
	void operator()(ostream& out, const VertexOrEdge& v) const {
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

inline set<pair<size_t, size_t>> chainToEdgeSet(const vector<size_t>& v){
    assert(v.size() > 1);
    set<pair<size_t,size_t>> out;
    for(size_t i=1; i<v.size(); i++){
        if(v[i-1] < v[i])
            out.insert(make_pair(v[i-1], v[i]));
        else
            out.insert(make_pair(v[i], v[i-1]));
    }
    if(v.front() < v.back())
        out.insert(make_pair(v.front(), v.back()));
    else
        out.insert(make_pair(v.back(), v.front()));
    return out;
}


void logCycle(vector<pair<vd,vd>>& c)
{
	for (auto & e : c)
	{
		LOG(DEBUG) << e.first << "--" << e.second << '\n';
	}
	LOG(DEBUG) << endline;
}


Cycle::Cycle(vector<pair<vd, vd>> edges_):edges(edges_)
{
	chain = cycleToChain(edges, [](const pair<vd,vd>&p){ return p; });
}

Cycle::Cycle(vector<vd> chain_):chain(chain_)
{
	auto eset = chainToEdgeSet(chain_);
	edges.resize(eset.size());
	copy(eset.begin(), eset.end(), edges.begin());
}

bool Cycle::intersects(const Cycle& that)const
{
   return intersection(that).size() > 0; // can be optimized futher
}

vector<pair<vd,vd>> Cycle::intersection(const Cycle& c)const
{
    vector<pair<vd,vd>> ret;
    auto& c1 = edges;
    auto& c2 = c.edges;
    set_intersection(c1.begin(), c1.end(), c2.begin(), c2.end(), 
        back_inserter(ret));
    return ret;
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
		for (auto e : edges)
		{
			auto ed = edge(e.first, e.second, g);
			g[ed.first].type = AROMATIC;
		}
	}
	return *this;
}


struct BfsNoop{
    void onVertex(size_t v){}
};

// assuming constant edge weight - no priority queue required
template<class Policy=BfsNoop>
struct ShortestPaths : Policy{
public:
    using G = ChemGraph;
    
    ShortestPaths(G& graph, size_t start, const vector<bool>& m=vector<bool>()):
    g(graph), visited(num_vertices(graph)), edgeTo(num_vertices(graph)), s(start), mask(m){
        queue.push_back(start);
        visited[start] = true;
        edgeTo[start] = 1<<(sizeof(size_t)*8-1);
        bfs();
    }
    // apply functor to each node along the shortest path from v to starting point 
    // except the starting point itself
    template<class Fn>
    void apply(size_t v, Fn&& fn, bool includeStart=false){
        if(!visited[v]) // cut apply for unexplored vertices
            return;
        size_t w = v;
        while(w != s){
            fn(w);
            w = edgeTo[w];
        }
        if(includeStart)
            fn(s);
    }
    // get shortest path to s
    vector<size_t> path(size_t v, bool includeStart=false){
        vector<size_t> vec;
        apply(v, [&vec](size_t w){ vec.push_back(w); }, includeStart);
        return vec;
    }

    // get prior vertex on path to this one
    size_t prev(size_t p){
        return edgeTo[p];
    }
private:
    void bfs(){
        while(!queue.empty()){
            size_t v = queue.front();
            queue.pop_front();
            // push all not visited
            auto adj = adjacent_vertices(v, g);
        	for(auto p = adj.first; p != adj.second; p++ ){
        		auto w = *p;
                if(!visited[w] && (mask.empty() || mask[w])){
                    edgeTo[w] = v;
                    visited[w] = true;
                    queue.push_back(w);
                }
            }
        }
    }
    G& g;
    vector<bool> visited;
    vector<size_t> edgeTo;
    deque<size_t> queue;
    size_t s;
    const vector<bool>& mask;
    
};

ShortestPaths<> shortestPaths(ChemGraph& g, size_t start, const vector<bool>& mask){
    return ShortestPaths<>(g, start, mask);
}


// depth-first search to detect all cycles
// ploicy-based design, final processing is deffered to the inherited policy
template<class Policy>
struct CycleDetect : Policy{
public:
    using G = ChemGraph;
    CycleDetect(G& graph, size_t start):Policy(), g(graph), visited(num_vertices(graph)){
        dfs(start);
    }

private:
    void dfs(size_t v){
        visited[v] = true;
        auto prev = path.size() ? path.back() : -1;
        path.push_back(v);
        auto adj = adjacent_vertices(v, g);
        for(auto p = adj.first; p != adj.second; p++ ){
        	auto w = *p;
            auto e = edge(v, w, g);
            // TODO: turn to predicate or just use Boost DFS
            if(g[e.first].type >= STEREO)
                continue;
            if(!visited[w]){
                dfs(w); //pushes w on path
                path.pop_back();
            }
            else if(w != prev){ // visited and not previous one
                auto i = find(path.begin(), path.end(), w);
                if(i != path.end()) //cycle
                    Policy::onCycle(i, path.end());
            }
        }
    }
    G& g;
    vector<bool> visited;
    vector<int> path;
};

struct FetchCycles{
    vector<vector<size_t>> cycles;
    template<class I>
    void onCycle(I beg, I end){
        cycles.emplace_back(beg, end);
    }
};

//some cycle basis not even Horton's cycle basis
vector<vector<size_t>> cycleBasis(ChemGraph& g, size_t start=0){
    return CycleDetect<FetchCycles>(g, start).cycles;
}

vector<Cycle> minimalCycleBasis(ChemGraph& g){
    vector<vector<size_t>> isolatedCycles; // 
    vector<vector<size_t>> basisCandidates; // that are not isolated
    auto const V = num_vertices(g);
    vector<vd> inCycle(V); // >0 if v is on some cycle
    vector<bool> inCycleSystem(V); // if v is in some cycle system
    LOG(DEBUG)<<"IN CYCLE: "<<inCycle<<endline;
    // find isolated cycles and enumerate vertices belonging to some cycle
    {
        vector<vector<size_t>> someCycles = cycleBasis(g);
        LOG(DEBUG) << "G SIZE:" << V << endline << "CYCLES ARE:" << endline;
        for(auto &c : someCycles)
        	LOG(DEBUG) << c << endline;
        // count number of times a cycles passes through a vertex
        for(auto& c : someCycles){
            for(size_t v : c){
                inCycle[v]++;
                inCycleSystem[v] = 1; // consider every as non-isoalted
            }
        }
        
        // test each cycle
        for(auto& c : someCycles){
            int total = accumulate(c.begin(), c.end(), 0, 
                [&inCycle](int sum, size_t v){
                return sum + inCycle[v];
            });
            // no other cycles passing though isolated cycle's vertices
            // therefore sum == length of cycle
            if(total == (int)c.size()){
                isolatedCycles.push_back(c);
                for(auto v : c) // filter out isolated cycles from cycle systems
                    inCycleSystem[v] = 0;
            }
            
        }
    }
    LOG(DEBUG)<<"IN CYCLE: "<<inCycle<<endline;

    // from now on consider only vertices from some cycles
    // find connection points - vertices with > 2 adjacent
    vector<size_t> connPts; 
    for(size_t v=0; v<num_vertices(g); v++){
    	if(inCycle[v] <= 1)
            continue; //vertex in an isolated cycle
        auto adj = adjacent_vertices(v, g);
        size_t neib = count_if(adj.first, adj.second,
            [&inCycle](vd a){
                return inCycle[a] > 0;
        });
        if(neib > 2){
            connPts.push_back(v);
        }
    }
    LOG(DEBUG) << "Conn pts.:" << connPts <<endline;
    // Modified Horton's algorithm (1987)
    // Each cycle in minimal cycle base
    // has form : Puw + Pvw + {u, v} (where Pxy is shortest path from x -> y)
    // Puw & Pvw = {w} (the paths have no common edges)     
    // Find candidates using these principles

    // Combined with knowledge of connection points, all points being
    // on some cycle & no isolated cycles we can examine only connection points
    // for shortest-path trees
    {
        for(size_t con : connPts){
            auto spf = shortestPaths(g, con, inCycleSystem);
            auto eds = edges(g);
            for(auto ep = eds.first ; ep != eds.second; ep++){
            	auto a = target(*ep, g);
            	auto b = source(*ep, g);
                if(!inCycleSystem[a] || !inCycleSystem[b])
                    continue;
                // drop these along the paths
                if(spf.prev(a) == b || spf.prev(b) == a)
                    continue;
                auto pa = spf.path(a);
                auto pb = spf.path(b);
                LOG(DEBUG) << "Two shortest paths from " << con << ": A "<< a 
                    << "  B " << b << endline;
                LOG(DEBUG) << pa << endline;
                LOG(DEBUG) << pb << endline;
                // both paths w/o starting point 'con'
                size_t missingLink = con;
                if(pa.empty() || pb.empty())
                    continue;
                // if have common suffix - drop
                // they must pass through some connection point
                // that we are going to process anyway
                if(pa.back() == pb.back())
                    continue;
                // add missing link
                pa.push_back(missingLink);
                // and follow 2nd path 
                for_each(pb.rbegin(), pb.rend(), [&pa](size_t v){
                    pa.push_back(v);
                });
                basisCandidates.push_back(pa);
            }
        }
    }
    // now perform Gaussian ellimination of candidate cycles
    sort(basisCandidates.begin(), basisCandidates.end(), 
        [](const vector<size_t>& v1, const vector<size_t>& v2){
            return v1.size() < v2.size();
    });
    vector<set<pair<size_t,size_t>>> basis; // basis accumulated as sets of edges
    vector<size_t> basisIdx; // indices of these candidates that are in final basis
    for(size_t i=0; i<basisCandidates.size(); i++){
        auto s = chainToEdgeSet(basisCandidates[i]);
        if(!elimination(s, basis)){
            basis.push_back(s);
            basisIdx.push_back(i);
        }
    }
    LOG(DEBUG) << "Isolated cycles:\n";
    for(auto & c : isolatedCycles){
        LOG(DEBUG) << c << endline;
    }
    LOG(DEBUG)<<"Cycle base candidates:\n";
    for(auto & c : basisCandidates){
        LOG(DEBUG) << c << endline;
    }
    using std::move;
    vector<Cycle> results;
    copy(isolatedCycles.begin(), isolatedCycles.end(), back_inserter(results));
    for(size_t i : basisIdx){
        results.push_back(Cycle{basisCandidates[i]});
    }
    LOG(DEBUG)<<"Minimal cycle base :\n";
    for(auto& v : results){
        LOG(DEBUG) << v << endline;
    }
    return results;
}