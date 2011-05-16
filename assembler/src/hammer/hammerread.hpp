/*
 * hammerread.hpp
 *
 *  Created on: 11.05.2011
 *      Author: snikolenko
 */
 
#ifndef HAMMER_READ_HPP
#define HAMMER_READ_HPP

#include<cassert>
#include<cmath>
#include<string>
#include<iostream>
#include<sstream>
#include<fstream>
#include<cstdlib>
#include<vector>
#include<map>
#include<list>
#include<queue>
#include<cstdarg>
#include<algorithm>

using namespace std;

class HammerRead {
public:
	string id;
	string seq;
	int count;
	float freq;
	bool operator<(HammerRead x) const {
		return x.seq < seq;
	}
	bool operator==(string x) const {
		return x == this->seq;
	}
};

#endif 
