#include <assert.h>
#include <vector>
#include <string>
#include <unordered_map>

#include "periodic.h"

using namespace std;

struct Entry{
	int code;
	std::string symbol;
	Entry(int code_, std::string symbol_):
		code(code_), symbol(symbol_){}
};

vector<Entry> atomList;
unordered_map<string, int> symbolToIndex;
//list and mapping for wildcard atoms
vector<Entry> wildList;
unordered_map<string, int> wildIndex;
//
static int create(string sym)
{
	int code = (int)atomList.size();
	atomList.emplace_back(code, sym);
	symbolToIndex.insert(make_pair(sym, code));
	return code;
}

//wildcards are all negative starting with -1
static int createWild(string sym)
{
	int code = (int)wildList.size()+1;
	wildList.emplace_back(-code, sym);
	wildIndex.insert(make_pair(sym, -code));
	return code;
}


struct Module{
	Module(){
		create("H");
		create("C");
		createWild("Z");
		createWild("R");
	}
};

//seed table with predefiend atoms
static Module init_;

Code::Code(const string& sym)
{
	auto it = symbolToIndex.find(sym);
	if (it == symbolToIndex.end())
	{
		it = wildIndex.find(sym);
		if (it == wildIndex.end())
			index = create(sym);
		else
			index = it->second;
	}
	else
		index = it->second;
}

Code::Code(int code)
{
	assert(code < 0 || code < (int)atomList.size());
	index = code;
}

static bool wildMatch(int wildcard, int atom)
{
	switch (wildcard){
	case Z:
		return true;
	case R:
		return atom != H;
	default:
		assert(false);
		return false;
	}
}

bool Code::matches(int code)const
{
	if (index >= 0)
	{
		if (code >= 0)
			return code == index;
		else
			return wildMatch(code, index);
	}
	if (code < 0) //2 wildcards - compare exactly
		return code == index;
	return wildMatch(index, code);
}

const string& Code::symbol()const
{
	return index >= 0 ? atomList[index].symbol : wildList[-index - 1].symbol;
}