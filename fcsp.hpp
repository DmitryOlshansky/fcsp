/*
 * fcsp.hpp
 *
 *  Created on: Dec 16, 2013
 *      Author: dmitry
 */
#pragma once

#include <vector>
#include <iostream>
#include <memory>
#include "descriptors.h"


struct FCSP {
	FCSP(std::vector<LevelOne> first, std::vector<LevelTwo> second);
	void load(std::istream& inp);
	void dumpGraph(std::ostream& dot);
	void linear(std::ostream& out);
	~FCSP();
private:
	struct Impl;
	std::unique_ptr<Impl> pimpl;
};
