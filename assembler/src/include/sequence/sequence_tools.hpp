#ifndef SEQUENCE_TOOLS_HPP_
#define SEQUENCE_TOOLS_HPP_

#include <sstream>
#include <string>
#include <vector>

#include "sequence/nucl.hpp"
#include "sequence/sequence.hpp"
#include "levenshtein.hpp"

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

inline const Sequence MergeOverlappingSequences(const vector<const Sequence*>& ss, size_t overlap) {
	SequenceBuilder sb;
	if (ss.empty())
		return Sequence(); 
	sb.append(ss.front()->Subseq(0, overlap));
	for (auto it = ss.begin(); it != ss.end(); ++it) {
		sb.append((*it)->Subseq(overlap));
	}
	return sb.BuildSequence();
}

inline size_t EditDistance(const Sequence& s1, const Sequence& s2) {
	return edit_distance(s1.str(), s2.str());
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

class EnsureEndsPositionAligner {
private:
	size_t upper_length_;
	size_t lower_length_;
public:
	EnsureEndsPositionAligner(size_t upper_length, size_t lower_length) :
		upper_length_(upper_length), lower_length_(lower_length) {
	}

	size_t GetPosition(size_t upper_position) {
		if (lower_length_ == 1)
			return 0;
		return (2 * upper_position * lower_length_ + upper_length_) / (2
				* upper_length_);
	}
};

#endif /* SEQUENCE_TOOLS_HPP_ */
