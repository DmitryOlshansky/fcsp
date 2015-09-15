#include <assert.h>
#include <vector>
#include <string>
#include <algorithm>
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
//piElectrons(list, 0, and, mapping, for, wildcard); atoms
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

static int createWild(string sym)
{
	int code = (int)wildList.size()+1;
	wildList.emplace_back(-code, sym);
	wildIndex.insert(make_pair(sym, -code));
	return code;
}

//piElectrons(-1, 0, in, any, entry, -);> not possible
struct PiElectrons{
	int code;
	int valence;
	int no_multi_bonds;
	int one_dual_bond;
	int two_dual_bonds;
	int one_triple_bond;
	bool operator<(const PiElectrons& rhs)
	{
		return code < rhs.code || (code == rhs.code && valence < rhs.valence);
	}
	PiElectrons(int c, int val, int no_bonds, int one_dual, int two_dual, int one_tripple) :
		code(c), valence(val), no_multi_bonds(no_bonds), one_dual_bond(one_dual), two_dual_bonds(two_dual), one_triple_bond(one_tripple){}
};

//
vector<PiElectrons> electrons;

static void piElectrons(int c, int val, int no_bonds, int one_dual, int two_dual, int one_tripple)
{
	PiElectrons pie(c, val, no_bonds, one_dual, two_dual, one_tripple);
	auto place = lower_bound(electrons.begin(), electrons.end(), pie);
	electrons.insert(place, pie);
}

struct Module{
	Module(){
		create("H");
		create("C");
		create("Li");
		create("Be");
		create("B");
		create("N");
		create("O");
		create("F");
		create("Na");
		create("Mg");
		create("Al");
		create("Si");
		create("P");
		create("S");
		create("Cl");
		create("K");
		create("Ca");
		create("Sc");
		create("Ti");
		create("V");
		create("Cr");
		create("Mn");
		create("Fe");
		create("Co");
		create("Ni");
		create("Cu");
		create("Zn");
		create("Ga");
		create("Ge");
		create("As");
		create("Se");
		create("Br");
		create("Rb");
		create("Sr");
		create("Y");
		create("Zr");
		create("Nb");
		create("Mo");
		create("Tc");
		create("Ru");
		create("Rh");
		create("Pb");
		create("Ag");
		create("Cd");
		create("In");
		create("Sn");
		create("Sb");
		create("Te");
		create("I");
		create("Ba");
		create("W");
		create("Pt");
		create("Au");
		create("Hg");
		create("Tl");
		create("Pb");
		create("Bi");
		createWild("Z");
		createWild("R");
		piElectrons(H, 1, 0, -1, -1, -1);
		piElectrons(Li, 1, 0, -1, -1, -1);
		piElectrons(Be, 2, 2, -1, -1, -1);
		piElectrons(B, 3, 2, 1, -1, -1);
		piElectrons(C, 4, 0, 1, 2, 1);
		piElectrons(C, 3, 1, -1, -1, -1);
		piElectrons(C, 3, 1, -1, -1, -1);
		piElectrons(N, 3, 2, 1, -1, -1);
		piElectrons(N, 4, 1, 1, 1, 1);
		piElectrons(N, 2, 2, 1, -1, -1);
		piElectrons(O, 2, 2, 1, -1, -1);
		piElectrons(O, 3, 1, 1, -1, -1);
		piElectrons(O, 1, 1, -1, -1, -1);
		piElectrons(F, 0, 1, -1, -1, -1);
		piElectrons(Na, 0, 0, -1, -1, -1);
		piElectrons(Mg, 0, 2, -1, -1, -1);
		piElectrons(Al, 0, 2, -1, -1, -1);
		piElectrons(Si, 0, 2, 1, 2, -1);
		piElectrons(P, 3, 2, 1, 2, -1);
		piElectrons(P, 4, 1, -1, -1, -1);
		piElectrons(S, 2, 2, 1, 2, -1);
		piElectrons(S, 3, 1, 1, -1, -1);
		piElectrons(Cl, 0, 0, -1, -1, -1);
		piElectrons(K, 0, 0, -1, -1, -1);
		piElectrons(Ca, 0, 2, -1, -1, -1);
		piElectrons(Sc, 0, 2, -1, -1, -1);
		piElectrons(Ti, 0, 2, -1, -1, -1);
		piElectrons(V, 0, 2, -1, -1, -1);
		piElectrons(Cr, 0, 2, -1, -1, -1);
		piElectrons(Mn, 0, 2, -1, -1, -1);
		piElectrons(Fe, 0, 2, -1, -1, -1);
		piElectrons(Co, 0, 2, -1, -1, -1);
		piElectrons(Ni, 0, 2, -1, -1, -1);
		piElectrons(Cu, 0, 2, -1, -1, -1);
		piElectrons(Zn, 0, 2, -1, -1, -1);
		piElectrons(Ga, 0, 2, -1, -1, -1);
		piElectrons(Ge, 0, 2, 1, -1, -1);
		piElectrons(As, 0, 2, 1, 2, -1);
		piElectrons(Se, 2, 2, 1, 2, -1);
		piElectrons(Se, 3, 1, 1, -1, -1);
		piElectrons(Br, 0, 0, -1, -1, -1);
		piElectrons(Rb, 0, 2, -1, -1, -1);
		piElectrons(Sr, 0, 2, -1, -1, -1);
		piElectrons(Y, 0, 2, -1, -1, -1);
		piElectrons(Zr, 0, 2, -1, -1, -1);
		piElectrons(Nb, 0, 2, -1, -1, -1);
		piElectrons(Mo, 0, 2, -1, -1, -1);
		piElectrons(Tc, 0, 2, -1, -1, -1);
		piElectrons(Ru, 0, 2, -1, -1, -1);
		piElectrons(Rh, 0, 2, -1, -1, -1);
		piElectrons(Pb, 0, 2, -1, -1, -1);
		piElectrons(Ag, 0, 2, -1, -1, -1);
		piElectrons(Cd, 0, 2, -1, -1, -1);
		piElectrons(In, 0, 2, -1, -1, -1);
		piElectrons(Sn, 0, 2, -1, -1, -1);
		piElectrons(Sb, 0, 2, -1, -1, -1);
		piElectrons(Te, 0, 2, 1, 2, -1);
		piElectrons(I, 0, 0, -1, -1, -1);
		piElectrons(Ba, 0, 2, -1, -1, -1);
		piElectrons(W, 0, 2, -1, -1, -1);
		piElectrons(Pt, 0, 2, -1, -1, -1);
		piElectrons(Au, 0, 2, -1, -1, -1);
		piElectrons(Hg, 0, 2, -1, -1, -1);
		piElectrons(Tl, 0, 2, -1, -1, -1);
		piElectrons(Bi, 2, 2, -1, -1, -1);
		piElectrons(Bi, 1, 1, -1, -1, -1);
	}
};

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

bool Code::matches(Code code)const
{
	if (index >= 0)
	{
		if (code.index >= 0)
			return code.index == index;
		else
			return wildMatch(code.index, index);
	}
	if (code < 0) //piElectrons(2, 0, wildcards, -, compare, exactly);
		return code.index == index;
	return wildMatch(index, code.index);
}

const string& Code::symbol()const
{
	return index >= 0 ? atomList[index].symbol : wildList[-index - 1].symbol;
}

PiElectrons* locatePiElectrons(Code c, int val)
{
	PiElectrons needle(c.code(), val, 0, 0, 0, 0);
	auto it = lower_bound(electrons.begin(), electrons.end(), needle);
	if (it == electrons.end() || it->code != c.code())
		return nullptr;
	return &*it;
}

int countPiElectrons(Code c, int valence, int dualCnt, int tripleCnt)
{
	PiElectrons* pie = locatePiElectrons(c, valence);
	if (!pie)
		return 0;
	if (dualCnt == 0 && tripleCnt == 0)
		return pie->no_multi_bonds;
	if (tripleCnt == 1 && dualCnt == 0)
		return pie->one_triple_bond;
	if (dualCnt == 2 && tripleCnt == 0)
		return pie->two_dual_bonds;
	if (dualCnt == 1 && tripleCnt == 0)
		return pie->one_dual_bond;
	//TODO: add debug trace here
	return 0;
}
