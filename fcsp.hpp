/*
 * fcsp.hpp
 *
 *  Created on: Dec 16, 2013
 *      Author: dmitry
 */
#include <vector>
#include <unordered_map>
#include <iostream>
#include <memory>
#include "descriptors.h"

#pragma once

struct FCSP {
	FCSP(std::unordered_map<int, LevelOne> first, std::unordered_map<int, LevelTwo> second);
	void load(std::istream& inp);
	void dumpGraph(std::ostream& dot);
	void linear(std::ostream& out);
	~FCSP();
private:
	struct Impl;
	std::unique_ptr<Impl> pimpl;
};
