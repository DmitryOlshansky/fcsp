#include <algorithm>
#include <stdexcept>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <ctype.h>
#include <cmath>
#include <string.h>
#include "ctab.h"
#include "parser.h"

using namespace std;

void error(const string& msg)
{
	throw logic_error(msg);
}

void warning(const string& msg)
{
	cerr << msg.c_str() << endl;
}

const char* fetchSpec(const char* fmt, FormatSpec& spec)
{
	char marker = *fmt;
	if (!marker)
		error("extra arguments in format string");
	//parse the width-pattern
	if (isalpha(marker) || isdigit(marker))
	{
		// an alpha repeated n times			
		const char* q = fmt;
		while (*++fmt == marker){}
		//fmt is next char
		auto n = fmt - q;
		//folowed by dot ? 
		if (*fmt == '.')
		{
			++fmt;
			if (*fmt != marker)
				error("format error - expected same alpha after '.'");
			const char* q = fmt;
			while (*++fmt == marker){}
			auto m = fmt - q;
			spec.n = n;
			spec.m = m;
		}
		else
		{
			spec.n = n;
			spec.m = 0;
		}
	}
	else
	{
		spec.n = 0;
		spec.m = 0;
	}
	return fmt;
}

CTab readMol(Parser& parser)
{
	CTab tab;
	//MOL Header
	tab.name = parser.line();
	tab.descr = parser.line();
	tab.comment = parser.line();
	//The Counts Line
	//aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv
	//aaa = number of atoms (current max 255)* [Generic]
	//bbb = number of bonds (current max 255)* [Generic]
	//lll = number of atom lists (max 30)* [Query]
	//fff = (obsolete) 
	//ccc = chiral flag: 0=not chiral, 1=chiral [Generic]
	//sss = number of stext entries [ISIS/Desktop]
	//xxx = (obsolete)
	//rrr = (obsolete)
	//ppp = (obsolete)
	//iii = (obsolete)
	//mmm = number of lines of additional properties,  including the M END line.
	// No longer supported, the default is set to 999.
	int aaa, bbb, lll, fff, ccc, sss, xxx, rrr, ppp, iii, mmm;
	char ver[8];
	parser.matchfln("aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv", aaa, bbb, lll, fff,
		ccc, sss, xxx, rrr, ppp, iii, mmm, ver);
	tab.chiral = ccc;
	tab.atomLists = lll;
	if(strcmp(ver, "V2000") != 0)
		warning("counts line has wrong version:"+string(ver));
	tab.atoms = vector<AtomEntry>(aaa);
	for(int i=0; i<aaa; i++)
	{
		//atom line
		//xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee
		double x, y, z;
		char symbol[4]; //aaa
		int ddd, ccc, sss, hhh, bbb, vvv, HHH, rrr, iii, mmm, nnn, eee;
		parser.matchfln("xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee",
			x, y, z, symbol, ddd, ccc, sss, hhh, bbb, vvv, HHH,
			rrr, iii, mmm, nnn, eee);
		auto white = find(symbol, symbol + 4, ' ');
		if (white != symbol + 4)
			*white = 0;
		Code code(symbol);
		tab.atoms[i] = AtomEntry(x, y, z, hhh, code);
	}
	tab.bounds = vector<BoundEntry>(bbb);
	for(int i=0; i<bbb; i++)
	{
		//bound line
		//111222tttsssxxxrrrccc
		int first, second;
		int ttt, sss, xxx, rrr, ccc;
		parser.matchfln("111222tttsssxxxrrrccc", first, second,
			ttt, sss, xxx, rrr, ccc);
		if(first < 0 || second < 0 ||
				first > (int)tab.atoms.size() || second > (int)tab.atoms.size())
			error("bad bounds indices - out of range.");
		tab.bounds[i] = BoundEntry(first, second, ttt);
	}
	// old-style properties M  PROP_NAME ......
	for(;;)
	{
		auto s = parser.line();
		if(s.size() > 6 && s.substr(0, 6) == "M  CHG")
		{
			istringstream iss(s.substr(6)); // continue parsing
			// FIXME: cross fingers and pray that charge and atom count doesn't go up to 100+
			int a, b, c;
			iss >> a >> b >> c;
			b -= 1;
			if(b >= tab.atoms.size() || b < 0)
				error("bad charge record - atom number out of range");
			tab.atoms[b].code.charge(c);
			cerr << "Found charge on "<<tab.atoms[b].code.symbol()<< " = "<< c <<endl;
		}
		else if(s == "M  END" || parser.eof())
			break;
		// cout << "Skipping: " << s << endl;
	}
	return tab;
}

CTab readMol(istream& inp)
{
	Parser parser(inp);
	return readMol(parser);
}

void writefln(ostream& out, const char* fmt)
{
	out << fmt << endl;
}

template<class A>
void write(ostream&out, A a, int n, int m)
{
	error("unsupported type combination");
}

void write(ostream& out, double a, int n, int m)
{
	out << setw(n+m+1)<<setiosflags(ios::right)
		<< setprecision(m) << setiosflags(ios::fixed) << a;
}

template<class F, class ...T>
void writefln(ostream& out, const char* fmt, F front, T... tail)
{
	while(*fmt)
	{
		FormatSpec spec;
		fmt = fetchSpec(fmt, spec);
		if (spec.n)
		{
			if (spec.m)
			{
				write(out, front, spec.n, spec.m);				
			}
			else
			{
				out.width(spec.n);
				out << front;
			}
			return writefln(out, fmt, tail...);
		}
		else
			out << (char)*fmt;
		fmt++;
	}
}

void writeMol(CTab& tab, ostream& out)
{
	out << tab.name << endl << tab.descr << endl << tab.comment << endl;
	//aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv
	writefln(out, "aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv", 
		tab.atoms.size(), tab.bounds.size(), tab.atomLists, 
		0, tab.chiral, 0, 0, 0, 0, 0, 999, "V2000");
	for (auto &e : tab.atoms)
	{
		writefln(out, "xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee",
			e.x, e.y, e.z, e.code.symbol(), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	}
	for (auto &e : tab.bounds)
	{
		writefln(out, "111222tttsssxxxrrrccc", e.a1, e.a2, e.type, 0, 0, 0, 0);
	}
	out << "M  END" << endl;
}

vector<SDF> readSdf(Parser& parser)
{
	vector<SDF> ret;
	while(!parser.eof())
	{
		SDF sdf;
		sdf.mol = readMol(parser);
		// read properties
		for(;;)
		{
			if(parser.eof()) //no properties
				break;
			auto line = parser.line();
			if(line == "$$$$")
				break;
			//> <property_name>
			if(line[0] != '>')
				error("SDF property is expected to start with >:"+line);
			auto left = line.find('<', 1);
			if(left == string::npos)
				error("expected '<' before property name");
			auto right = line.find('>', left+1);
			if(right == string::npos)
				error("expected '>' after property name");
			auto name = line.substr(left+1, right-left-1);
			sdf.props[name] = vector<string>();
			auto& target = sdf.props[name];
			// read property lines
			for(;;)
			{
				line = parser.line();
				//empty line breaks property value list
				if(line.empty())
					break;
				target.push_back(std::move(line));
			}
		}
		ret.push_back(std::move(sdf));
	}
	return ret;
}

vector<SDF> readSdf(std::istream& inp)
{
	Parser parser(inp);
	return readSdf(parser);
}
