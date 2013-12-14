#include <string>
#include <sstream>
#include <iostream>
#include "descriptors.h"
#include "parser.h"

using namespace std;

std::unordered_map<int, LevelOne> read1Order(std::istream& inp)
{
	unordered_map<int, LevelOne> descr;
	string s;
	Parser parser(inp);
	while (!parser.eof())
	{
		string sym = parser.quotedString();
		string valency = parser.quotedString();
		stringstream dcStr(parser.line()); //rest of the line
		int dc;
		cout << valency << endl;
		dcStr >> dc;
		Code code = atomCode(sym);
		vector<int> v;		
		stringstream str(valency);
		for (;;)
		{
			int k;
			char delim;
			str >> k;			
			v.push_back(k);
			str >> delim;
			if (!str.good() || delim != ',')
				break;
		}
		descr.emplace(code, LevelOne(code, move(v), dc));
	}
	return descr;
}
