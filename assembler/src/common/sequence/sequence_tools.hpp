#ifndef SEQUENCE_TOOLS_HPP_
#define SEQUENCE_TOOLS_HPP_

#include <sstream>
#include <string>
#include <vector>

#include "nucl.hpp"

inline std::string Reverse(const std::string &s) {
	return std::string(s.rbegin(), s.rend());
}

inline std::string Complement(const std::string &s) {
	std::string res = s;
	for (size_t i = 0; i < s.size(); i++) {
		if (res[i] != 'N' && res[i] != 'n') {
			res[i] = nucl_complement(res[i]);
		}
	}
	return res;
}

inline std::string ReverseComplement(const std::string &s) {
	return Complement(Reverse(s));
}


#endif /* SEQUENCE_TOOLS_HPP_ */
