#pragma once
#include <string>

enum {
//predefined
	H = 0,
	C = 1,
//wild cards
	Z = -1, //anything
	R = -2, //anything but H
};

//an atom or a template
struct Code {
public:
	Code() :index(-99999){}
	//create/get from symbol code
	explicit Code(int c);
	//create/get from symbol name
	Code(const std::string& symbol);
	//matches given index
	bool matches(int code)const;
	//ref to a symbol in a global table
	const std::string& symbol()const;
	int code()const{ return index; }
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
};

namespace std {
	template <> struct hash<Code>
	{
		size_t operator()(Code x) const
		{
			return hash<int>()(x.code());
		}
	};
}