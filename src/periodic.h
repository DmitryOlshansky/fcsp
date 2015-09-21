#pragma once
#include <string>

enum {
//predefined
	H = 0,
	C = 1,
	Li,
	Be,
	B,
	N,
	O,
	F,
	Na,
	Mg,
	Al,
	Si,
	P,
	S,
	Cl,
	K,
	Ca,
	Sc,
	Ti,
	V,
	Cr,
	Mn,
	Fe,
	Co,
	Ni,
	Cu,
	Zn,
	Ga,
	Ge,
	As,
	Se,
	Br,
	Rb,
	Sr,
	Y,
	Zr,
	Nb,
	Mo,
	Tc,
	Ru,
	Rh,
	Pb,
	Ag,
	Cd,
	In,
	Sn,
	Sb,
	Te,
	I,
	Ba,
	W,
	Pt,
	Au,
	Hg,
	Tl,
	Bi,
//wild cards
	Z = -1, //anything
	R = -2, //anything but H
};

//an atom or a template
struct Code {
public:
	Code() :index(-99999), _charge(0){}
	//create/get from symbol code
	Code(int c);
	//create/get from symbol name
	Code(const std::string& symbol);
	//matches given index
	bool matches(Code code)const;
	//ref to a symbol in a global table
	const std::string& symbol()const;
	int code()const{ return index; }
	int charge()const{ return _charge; }
	Code& charge(int chargeMod){ _charge = chargeMod; return *this; }
	bool operator<(Code rhs)const{ return index < rhs.index; }
	bool operator<=(Code rhs)const{ return index <= rhs.index; }
	bool operator>(Code rhs)const{ return index > rhs.index; }
	bool operator>=(Code rhs)const{ return index >= rhs.index; }
	bool operator==(Code rhs)const{ return index == rhs.index; }
	bool operator!=(Code rhs)const{ return index != rhs.index; }
	bool operator==(int rhs)const{ return index == rhs; }
	bool operator!=(int rhs)const{ return index != rhs; }
private:
	int index; //index or wild-card
	int _charge; // charge modifier
};

int countPiElectrons(Code c, int valence, int dualCnt, int tripleCnt);

namespace std {
	template <> struct hash<Code>
	{
		size_t operator()(Code x) const
		{
			return hash<int>()(x.code());
		}
	};
}