/**
	Small module for string<->some type<->string conversions.

	API inspired by D's std.conv.
*/
#pragma once

#include <sstream>

template<class T, class U> class Converter;

template<class T>
class Converter<T, std::string>{
public:
	static T convert(const std::string& s)
	{
		T val;
		std::stringstream str(s);
		str >> val;
		return val;
	}
};

template<class T>
class Converter<std::string,T>{
public:
	static std::string convert(const T& val)
	{
		std::stringstream str;
		str << val;
		return str.str(); // does make a copy
	}
};

template<class T, class U>
T to(U from)
{
	return Converter<T,U>::convert(from);
}
