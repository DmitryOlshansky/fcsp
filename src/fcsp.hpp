/*
 * fcsp.hpp
 *
 *  Created on: Dec 16, 2013
 *      Author: dmitry
 */
#pragma once

#include <stdexcept>
#include <vector>
#include <string>
#include <iostream>
#include <memory>
#include "periodic.hpp"
#include "ctab.hpp"
#include "descriptors.hpp"

enum FCSPFMT {
	JSON, // array of JSON arrays with pairs : (code,bindings)
	CSV, // CSV - 2 columns: file-name,codes
	TXT // TXT - line per file, whitespace separated codes
};

struct FCSPOptions{
	std::vector<LevelOne> first;
	std::vector<LevelTwo> second;
	std::vector<Replacement> replacements;
	bool long41;
	FCSPFMT format;
};

struct FCSP {
	FCSP(FCSPOptions opts);
	void load(std::istream& inp);
	void dumpGraph(std::ostream& dot);
	void process(std::ostream& out, std::string filename="");
	~FCSP();
private:
	struct Impl;
	std::unique_ptr<Impl> pimpl;
};
