#include <iostream>
#include <algorithm>
#include <map>
#include <string>
#include <vector>

#include "utils.hpp"

std::string reverseComplement(const std::string& read) {
	std::map<char, char> reverse;
	reverse['C'] = 'G';
	reverse['G'] = 'C';
	reverse['T'] = 'A';
	reverse['A'] = 'T';
	reverse['N'] = 'N';

	std::vector<char> res;
	for(int i = 0; i < (int) read.length(); ++i) {
		res.push_back(reverse[read[i]]);
	}

	std::reverse(res.begin(), res.end());
	return std::string(res.begin(), res.end());
}
