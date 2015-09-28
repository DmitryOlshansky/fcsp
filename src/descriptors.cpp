#include <algorithm>
#include <string>
#include <sstream>
#include "descriptors.h"
#include "parser.h"
#include "ctab.h"
#include "log.hpp"

using namespace std;

void read1stOrder(istream& inp, vector<LevelOne> &dest)
{
	string s;
	Parser parser(inp);
	while (!parser.eof())
	{
		string sym = parser.quotedString();
		string valency = parser.quotedString();
		stringstream dcStr(parser.line()); //rest of the line
		int dc;
		//cout << valency << endline;
		dcStr >> dc;
		Code code(sym);
		stringstream str(valency);
		for (;;)
		{
			int val;
			char delim;
			str >> val;
			LevelOne toInsert(code, val, dc);
			auto it = lower_bound(begin(dest), end(dest), toInsert);
			dest.insert(it, toInsert);
			str >> delim;
			if (!str.good() || delim != ',')
				break;
		}
	}
	//for (auto a : dest)
	//	cout << "*** " << a.center.symbol() << " " << a.valence << " " << "DC: " << a.dc << endline;
}

void read2ndOrder(istream& inp, vector<LevelTwo> &dest)
{
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
		if(rec.props.find("REPLONLY") != rec.props.end())
		{
			stringstream str2(rec.props["REPLONLY"].front());
			str2 >> replOnly;
		}
		if(rec.props.find("DC") == rec.props.end())
			LOG(ERROR) << "No DC found for level-2 pattern #"<<i<<endline;
		else
			for_each(rec.props["DC"].begin(), rec.props["DC"].end(), 
			[&dest, replOnly, c, valency, &links](string s){
				stringstream str(s);
				int dc;
				str >> dc;
				LevelTwo t(c, valency, replOnly, links, dc);
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
	
}
