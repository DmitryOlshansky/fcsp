#include <algorithm>
#include <string>
#include <sstream>
#include "descriptors.hpp"
#include "conv.hpp"
#include "ctab.hpp"
#include "log.hpp"
#include "parser.hpp"

using namespace std;

Replacement::Replacement(ChemGraph g, int dc_, int coupling_) :
	piece(std::move(g)), dc(dc_), coupling(coupling_)
{
	a1 = a2 = -1;
	auto asym = Code("A1");
	auto bsym = Code("A2");
	auto r = vertices(piece);
	for (auto i = r.first; i != r.second; i++)
	{
		if (piece[*i].code == asym)
		{
			a1 = *i;
			piece[*i].code = Code("R");
		}
		if (piece[*i].code == bsym)
		{
			a2 = *i;
			piece[*i].code = Code("R");
		}
	}
	if(a1 == -1 || a2 == -1)
	{
		throw std::logic_error("Bad replacement loaded");
	}
}

auto read1stOrder(istream& inp) -> vector<LevelOne>
{
	vector<LevelOne> dest;
	string s;
	Parser parser(inp);
	while (!parser.eof())
	{
		string sym = parser.quotedString();
		string valency = parser.quotedString();
		stringstream dcStr(parser.line()); //rest of the line
		int dc;
		dcStr >> dc;
		Code code(sym);
		stringstream str(valency);
		for (;;)
		{
			int val;
			char delim;
			str >> val;
			LevelOne toInsert{code, val, dc};
			auto it = lower_bound(begin(dest), end(dest), toInsert);
			dest.insert(it, toInsert);
			str >> delim;
			if (!str.good() || delim != ',')
				break;
		}
	}
	return dest;
}

auto read2ndOrder(istream& inp) -> vector<LevelTwo>
{
	vector<LevelTwo> dest;
	vector<SDF> sdf = readSdf(inp);
	for (size_t i=0; i<sdf.size(); i++)
	{
		auto& rec = sdf[i];
		auto& b = rec.mol.bounds;
		int center = 1;
		//atom with index 1 as center for 2 atom patterns
		if (rec.mol.atoms.size() > 2)
		{
			int f = b[0].a1;
			auto cnt = count_if(begin(b), end(b), [f](BoundEntry e){
				return e.a1 == f || e.a2 == f;
			});
			if (cnt == b.size()) //present in all bonds
			{
				center = f;
			}
			else
			{
				center = b[0].a2;
			}
		}
		Code c = rec.mol.atoms[center - 1].code;
		int valency = 0;
		vector<Linked> links;
		for_each(begin(b), end(b), [&valency, &rec, &links, center](BoundEntry e){
			valency += e.type;
			if (e.a1 == center)
				links.emplace_back(rec.mol.atoms[e.a2 - 1].code, e.type);
			else
				links.emplace_back(rec.mol.atoms[e.a1 - 1].code, e.type);
		});
		int replOnly = 0; // default to commonly usable DC
		int monolith = 0; // default to allow superpositon of DCs
		int start = 0;
		if(rec.props.find("REPLONLY") != rec.props.end())
		{
			replOnly = to<int>(rec.props["REPLONLY"].front());
		}
		if(rec.props.find("MONOLITH") != rec.props.end())
		{
			monolith = to<int>(rec.props["MONOLITH"].front());
		}
		if(rec.props.find("START") != rec.props.end())
		{
			start = to<int>(rec.props["START"].front());
		}
		if(rec.props.find("DC") == rec.props.end())
			LOG(ERROR) << "No DC found for level-2 pattern #"<<i<<endline;
		else
			for_each(rec.props["DC"].begin(), rec.props["DC"].end(), 
			[&](string s){
				int dc = to<int>(s);
				LevelTwo t{c, start, valency, monolith != 0, replOnly != 0, links, dc};
				auto lb = lower_bound(begin(dest), end(dest), t);
				dest.insert(lb, t);
			});
	}
	
	for (auto &e : dest)
	{
		LOG(DEBUG) << "^^^ " << e.center.symbol() << " " << e.valence << "{ ";
		for (auto lnk : e.bonds)
			LOG(DEBUG) << lnk.atom.symbol() << " ";
		LOG(DEBUG) << "} DC: " << e.dc << endline;
	}
	return dest;
}


auto readReplacements(istream& inp) -> vector<Replacement>
{
	vector<Replacement> repls;
	auto sdfs = readSdf(inp);
	for_each(sdfs.begin(), sdfs.end(), [&repls](SDF& sdf){
		int dc = to<int>(sdf.props["DC"][0]);
		int couple = to<int>(sdf.props["COUPLING"][0]);
		repls.emplace_back(toGraph(sdf.mol), dc, couple);
	});
	return repls;
}
