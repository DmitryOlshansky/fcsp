#include <algorithm>
#include <string>
#include <sstream>
#include <iostream>
#include "descriptors.h"
#include "parser.h"
#include "ctab.h"

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
		//cout << valency << endl;
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
	for (auto a : dest)
		cout << "*** " << a.center.symbol() << " " << a.valence << " " << "DC: " << a.dc << endl;
}

void read2ndOrder(istream& inp, vector<LevelTwo> &dest)
{
	vector<SDF> sdf = readSdf(inp);
}
