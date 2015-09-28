// An output sink to be used in place of stderr
// with various verbosity levels
#pragma once
#include <iostream>
#include <iomanip>

#define LOG(level) LogSink(level)

extern int logLevel;

 // C++'s endl is boatload of template magic, so define our own
#define endline EndLine()

struct EndLine{};

struct LogSink{
	int level;
	LogSink(int lvl): level(lvl){}
};

template<class T>
LogSink operator<<(LogSink sink, const T& arg){
	if(sink.level <= logLevel)
		std::cerr << arg;
	return sink;
}

// do conditional "magic" for our LogSink
inline LogSink operator<<(LogSink sink, EndLine line){
	if(sink.level <= logLevel)
		std::cerr << std::endl;
	return sink;
}

// pass-through for other streams
template<class Stream>
Stream& operator<<(Stream& sink, EndLine line){
	sink << std::endl;
	return sink;
}

enum {
	FATAL = 1,
	ERROR = 2,
	WARN = 3,
	INFO = 4,
	DEBUG = 5,
	TRACE = 6
};
