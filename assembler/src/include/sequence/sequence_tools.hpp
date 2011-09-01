#ifndef SEQUENCE_TOOLS_HPP_
#define SEQUENCE_TOOLS_HPP_

#include <sstream>
#include <string>
#include <vector>

#include "sequence/nucl.hpp"

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

class UniformPositionAligner {
private:
	size_t upper_length_;
	size_t lower_length_;
public:
	UniformPositionAligner(size_t upper_length, size_t lower_length) :
		upper_length_(upper_length), lower_length_(lower_length) {
	}

	size_t GetPosition(size_t upper_position) {
		return (2 * upper_position + 1) * lower_length_ / (2 * upper_length_);
	}
};

#endif /* SEQUENCE_TOOLS_HPP_ */
