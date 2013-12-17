#include <algorithm>
#include <string.h>
#include "periodic.h"

using namespace std;

struct AtomData{
	const char *symbol, *name;
};

//Atom number --> data
AtomData data[] = {
	{ "H", "Hydrogen" },
	{ "He", "Helium" },
	{ "Li", "Lithium" },
	{ "Be", "Beryllium" },
	{ "B", "Boron" },
	{ "C", "Carbon" },
	{ "N", "Nitrogen" },
	{ "O", "Oxygen" },
	{ "F", "Fluorine" },
	{ "Ne", "Neon" },
	{ "Na", "Sodium" },
	{ "Mg", "Magnesium" },
	{ "Al", "Aluminium" },
	{ "Si", "Silicon" },
	{ "P", "Phosphorus" },
	{ "S", "Sulfur" },
	{ "Cl", "Chlorine" },
	{ "Ar", "Argon" },
	{ "K", "Potassium" },
	{ "Ca", "Calcium" },
	{ "Sc", "Scandium" },
	{ "Ti", "Titanium" },
	{ "V", "Vanadium" },
	{ "Cr", "Chromium" },
	{ "Mn", "Manganese" },
	{ "Fe", "Iron" },
	{ "Co", "Cobalt" },
	{ "Ni", "Nickel" },
	{ "Cu", "Copper" },
	{ "Zn", "Zinc" },
	{ "Ga", "Gallium" },
	{ "Ge", "Germanium" },
	{ "As", "Arsenic" },
	{ "Se", "Selenium" },
	{ "Br", "Bromine" },
	{ "Kr", "Krypton" },
	{ "Rb", "Rubidium" },
	{ "Sr", "Strontium" },
	{ "Y", "Yttrium" },
	{ "Zr", "Zirconium" },
	{ "Nb", "Niobium" },
	{ "Mo", "Molybdenum" },
	{ "Tc", "Technetium" },
	{ "Ru", "Ruthenium" },
	{ "Rh", "Rhodium" },
	{ "Pd", "Palladium" },
	{ "Ag", "Silver" },
	{ "Cd", "Cadmium" },
	{ "In", "Indium" },
	{ "Sn", "Tin" },
	{ "Sb", "Antimony" },
	{ "Te", "Tellurium" },
	{ "I", "Iodine" },
	{ "Xe", "Xenon" },
	{ "Cs", "Caesium" },
	{ "Ba", "Barium" },
	{ "La", "Lanthanum" },
	{ "Ce", "Cerium" },
	{ "Pr", "Praseodymium" },
	{ "Nd", "Neodymium" },
	{ "Pm", "Promethium" },
	{ "Sm", "Samarium" },
	{ "Eu", "Europium" },
	{ "Gd", "Gadolinium" },
	{ "Tb", "Terbium" },
	{ "Dy", "Dysprosium" },
	{ "Ho", "Holmium" },
	{ "Er", "Erbium" },
	{ "Tm", "Thulium" },
	{ "Yb", "Ytterbium" },
	{ "Lu", "Lutetium" },
	{ "Hf", "Hafnium" },
	{ "Ta", "Tantalum" },
	{ "W", "Tungsten" },
	{ "Re", "Rhenium" },
	{ "Os", "Osmium" },
	{ "Ir", "Iridium" },
	{ "Pt", "Platinum" },
	{ "Au", "Gold" },
	{ "Hg", "Mercury" },
	{ "Tl", "Thallium" },
	{ "Pb", "Lead" },
	{ "Bi", "Bismuth" },
	{ "Po", "Polonium" },
	{ "At", "Astatine" },
	{ "Rn", "Radon" },
	{ "Fr", "Francium" },
	{ "Ra", "Radium" },
	{ "Ac", "Actinium" },
	{ "Th", "Thorium" },
	{ "Pa", "Protactinium" },
	{ "U", "Uranium" },
	{ "Np", "Neptunium" },
	{ "Pu", "Plutonium" },
	{ "Am", "Americium" },
	{ "Cm", "Curium" },
	{ "Bk", "Berkelium" },
	{ "Cf", "Californium" },
	{ "Es", "Einsteinium" },
	{ "Fm", "Fermium" },
	{ "Md", "Mendelevium" },
	{ "No", "Nobelium" },
	{ "Lr", "Lawrencium" },
	{ "Rf", "Rutherfordium" },
	{ "Db", "Dubnium" },
	{ "Sg", "Seaborgium" },
	{ "Bh", "Bohrium" },
	{ "Hs", "Hassium" },
	{ "Mt", "Meitnerium" },
	{ "Ds", "Darmstadtium" },
	{ "Rg", "Roentgenium" },
	{ "Cn", "Copernicium" },
	{ "Uut", "Ununtrium" },
	{ "Fl", "Flerovium" },
	{ "Uup", "Ununpentium" },
	{ "Lv", "Livermorium" }
};

std::string atomSymbol(Code atom)
{
	auto idx = (atom & ELEMENT_MASK) - 1;
	if (idx >= sizeof(data) / sizeof(data[0]))
		return "?";
	string str = data[idx].symbol;
	if (atom & POSITIVE_ION)
		str += "+";
	else if (atom & NEGATIVE_ION)
		str += "-";
	return str;
}


struct AtomStruct{
	Code code;
	const char* str;
};

//sorted by symbol
AtomStruct symbols[] = {
		{ Ac, "Ac"},
		{ Ag, "Ag"},
		{ Al, "Al"},
		{ Am, "Am"},
		{ Ar, "Ar"},
		{ As, "As"},
		{ At, "At"},
		{ Au, "Au"},
		{ B, "B"},
		{ Ba, "Ba"},
		{ Be, "Be"},
		{ Bh, "Bh"},
		{ Bi, "Bi"},
		{ Bk, "Bk"},
		{ Br, "Br"},
		{ C, "C"},
		{ Ca, "Ca"},
		{ Cd, "Cd"},
		{ Ce, "Ce"},
		{ Cf, "Cf"},
		{ Cl, "Cl"},
		{ Cm, "Cm"},
		{ Cn, "Cn"},
		{ Co, "Co"},
		{ Cr, "Cr"},
		{ Cs, "Cs"},
		{ Cu, "Cu"},
		{ Db, "Db"},
		{ Ds, "Ds"},
		{ Dy, "Dy"},
		{ Er, "Er"},
		{ Es, "Es"},
		{ Eu, "Eu"},
		{ F, "F"},
		{ Fe, "Fe"},
		{ Fl, "Fl"},
		{ Fm, "Fm"},
		{ Fr, "Fr"},
		{ Ga, "Ga"},
		{ Gd, "Gd"},
		{ Ge, "Ge"},
		{ H, "H"},
		{ He, "He"},
		{ Hf, "Hf"},
		{ Hg, "Hg"},
		{ Ho, "Ho"},
		{ Hs, "Hs"},
		{ I, "I"},
		{ In, "In"},
		{ Ir, "Ir"},
		{ K, "K"},
		{ Kr, "Kr"},
		{ La, "La"},
		{ Li, "Li"},
		{ Lr, "Lr"},
		{ Lu, "Lu"},
		{ Lv, "Lv"},
		{ Md, "Md"},
		{ Mg, "Mg"},
		{ Mn, "Mn"},
		{ Mo, "Mo"},
		{ Mt, "Mt"},
		{ N, "N"},
		{ Na, "Na"},
		{ Nb, "Nb"},
		{ Nd, "Nd"},
		{ Ne, "Ne"},
		{ Ni, "Ni"},
		{ No, "No"},
		{ Np, "Np"},
		{ O, "O"},
		{ Os, "Os"},
		{ P, "P"},
		{ Pa, "Pa"},
		{ Pb, "Pb"},
		{ Pd, "Pd"},
		{ Pm, "Pm"},
		{ Po, "Po"},
		{ Pr, "Pr"},
		{ Pt, "Pt"},
		{ Pu, "Pu"},
		{ Ra, "Ra"},
		{ Rb, "Rb"},
		{ Re, "Re"},
		{ Rf, "Rf"},
		{ Rg, "Rg"},
		{ Rh, "Rh"},
		{ Rn, "Rn"},
		{ Ru, "Ru"},
		{ S, "S"},
		{ Sb, "Sb"},
		{ Sc, "Sc"},
		{ Se, "Se"},
		{ Sg, "Sg"},
		{ Si, "Si"},
		{ Sm, "Sm"},
		{ Sn, "Sn"},
		{ Sr, "Sr"},
		{ Ta, "Ta"},
		{ Tb, "Tb"},
		{ Tc, "Tc"},
		{ Te, "Te"},
		{ Th, "Th"},
		{ Ti, "Ti"},
		{ Tl, "Tl"},
		{ Tm, "Tm"},
		{ U, "U"},
		{ Uup, "Uup"},
		{ Uut, "Uut"},
		{ V, "V"},
		{ W, "W"},
		{ Xe, "Xe"},
		{ Y, "Y"},
		{ Yb, "Yb"},
		{ Zn, "Zn"},
		{ Zr, "Zr"}
};
#include <iostream>
Code atomCode(const string& symbol)
{
	if (symbol.empty())
		return Code(0);
	int start = 0;
	size_t len = symbol.size();
	if (symbol.back() == '-'){
		len -= 1;
		start += NEGATIVE_ION;
	}
	else if (symbol.back() == '+'){
		len -= 1;
		start += POSITIVE_ION;
	}
	auto it = std::lower_bound(symbols, symbols + sizeof(symbols) / sizeof(symbols[0]), symbol,
		[len](const AtomStruct& lhs, const string& rhs){
		return strncmp(lhs.str, rhs.c_str(), len) < 0;
	});
	if(it == symbols+sizeof(symbols))
		return Code(0);
	//cout << symbol.c_str() << "=" << it->code << endl;
	return (Code)((int)it->code+start);
}
