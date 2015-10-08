#pragma once
#include <string>

enum {
//predefined
	H_code = 0,
	C_code = 1,
	Li_code,
	Be_code,
	B_code,
	N_code,
	O_code,
	F_code,
	Na_code,
	Mg_code,
	Al_code,
	Si_code,
	P_code,
	S_code,
	Cl_code,
	K_code,
	Ca_code,
	Sc_code,
	Ti_code,
	V_code,
	Cr_code,
	Mn_code,
	Fe_code,
	Co_code,
	Ni_code,
	Cu_code,
	Zn_code,
	Ga_code,
	Ge_code,
	As_code,
	Se_code,
	Br_code,
	Rb_code,
	Sr_code,
	Y_code,
	Zr_code,
	Nb_code,
	Mo_code,
	Tc_code,
	Ru_code,
	Rh_code,
	Pb_code,
	Ag_code,
	Cd_code,
	In_code,
	Sn_code,
	Sb_code,
	Te_code,
	I_code,
	Ba_code,
	W_code,
	Pt_code,
	Au_code,
	Hg_code,
	Tl_code,
	Bi_code,
//wild cards
	Z_code = -1, // anything
	R_code = -2, // anything but H
	X_code = -3,  // anything but O
	Y1_code = -4 // heteroatom - anything but C or H
};

// An atom or a template a-la ^
// may hold a charge modifier
struct Code {
public:
	Code() :index(-99999), _charge(0){}
	//create/get from symbol code
	explicit Code(int c);
	// by wildcard code
	static int WildCard(int c);
	//create/get from symbol name
	Code(const std::string& symbol);
	//matches given index
	bool matches(Code code)const;
	//ref to a symbol in a global table
	const std::string& symbol()const;
	bool isWild()const;
	int code()const{ return index; }
	int charge()const{ return _charge; }
	Code& charge(int chargeMod){ _charge = chargeMod; return *this; }
	bool operator<(Code rhs)const{ return index < rhs.index; }
	bool operator<=(Code rhs)const{ return index <= rhs.index; }
	bool operator>(Code rhs)const{ return index > rhs.index; }
	bool operator>=(Code rhs)const{ return index >= rhs.index; }
	bool operator==(Code rhs)const{ return index == rhs.index; }
	bool operator!=(Code rhs)const{ return index != rhs.index; }
private:
	int index; //index or wild-card
	int _charge; // charge modifier
};

extern Code H;
extern Code C;
extern Code Li;
extern Code Be;
extern Code B;
extern Code N;
extern Code O;
extern Code F;
extern Code Na;
extern Code Mg;
extern Code Al;
extern Code Si;
extern Code P;
extern Code S;
extern Code Cl;
extern Code K;
extern Code Ca;
extern Code Sc;
extern Code Ti;
extern Code V;
extern Code Cr;
extern Code Mn;
extern Code Fe;
extern Code Co;
extern Code Ni;
extern Code Cu;
extern Code Zn;
extern Code Ga;
extern Code Ge;
extern Code As;
extern Code Se;
extern Code Br;
extern Code Rb;
extern Code Sr;
extern Code Y;
extern Code Zr;
extern Code Nb;
extern Code Mo;
extern Code Tc;
extern Code Ru;
extern Code Rh;
extern Code Pb;
extern Code Ag;
extern Code Cd;
extern Code In;
extern Code Sn;
extern Code Sb;
extern Code Te;
extern Code I;
extern Code Ba;
extern Code W;
extern Code Pt;
extern Code Au;
extern Code Hg;
extern Code Tl;
extern Code Bi;

extern Code Z;
extern Code R;


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