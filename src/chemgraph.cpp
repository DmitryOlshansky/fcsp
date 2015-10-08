#include <boost/graph/graphviz.hpp>
#include "chemgraph.hpp"
#include "ctab.hpp"

using namespace std;

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