//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef SEQUENCE_TOOLS_HPP_
#define SEQUENCE_TOOLS_HPP_

#include <sstream>
#include <string>
#include <vector>

#include "nucl.hpp"
#include "sequence.hpp"
#include "utils/levenshtein.hpp"

inline const std::string Reverse(const std::string &s) {
    return std::string(s.rbegin(), s.rend());
}

inline const std::string Complement(const std::string &s) {
    std::string res(s.size(), 0);
    transform(s.begin(), s.end(), res.begin(), nucl_complement);
    return res;
}

inline const Sequence MergeOverlappingSequences(std::vector<Sequence>& ss,
        size_t overlap, bool safe_merging = true) {
    if (ss.empty()) {
        return Sequence();
    }
    SequenceBuilder sb;
    Sequence prev_end = ss.front().Subseq(0, overlap);
    sb.append(prev_end);
    for (auto it = ss.begin(); it != ss.end(); ++it) {
        if(safe_merging)
            VERIFY(prev_end == it->Subseq(0, overlap));
        sb.append(it->Subseq(overlap));
        prev_end = it->Subseq(it->size() - overlap);
    }
    return sb.BuildSequence();
}

inline size_t EditDistance(const Sequence& s1, const Sequence& s2) {
    return edit_distance(s1.str(), s2.str());
}

inline bool Relax(int& val, int new_val) {
    if (new_val > val) {
        val = new_val;
        return true;
    }
    return false;
}

inline std::pair<size_t, size_t> LocalSimilarity(const Sequence& s1, const Sequence& s2) {
    size_t m = s1.size();
    size_t n = s2.size();
  std::vector<std::vector<int>> a(m + 1);
    for (size_t i = 0; i <= m; ++i) {
        a[i].resize(n + 1);
    }
    for (size_t i = 0; i <= m; ++i) {
        for (size_t j = 0; j <= n; ++j) {
            a[i][j] = 0;
        }
    }
    for (size_t i = 1; i <= m; ++i) {
        for (size_t j = 1; j <= n; ++j) {
            Relax(a[i][j], a[i - 1][j] - 1);
            Relax(a[i][j], a[i][j - 1] - 1);
            if (s1[i - 1] == s2[j - 1]) {
                Relax(a[i][j], a[i - 1][j - 1] + 1);
            } else {
                Relax(a[i][j], a[i - 1][j - 1] - 1);
            }
        }
    }

    //finding local alignment
    int answer = 0;
    size_t i_m = 0;
    size_t j_m = 0;
    for (size_t i = 0; i <= m; ++i) {
        for (size_t j = 0; j <= n; ++j) {
            if (Relax(answer, a[i][j])) {
                i_m = i;
                j_m = j;
            }
        }
    }

    //finding alignment lengths
    size_t i = i_m;
    size_t j = j_m;
    while (a[i][j] > 0) {
        if (a[i][j] == a[i][j - 1] - 1) {
            j--;
        } else if (a[i][j] == a[i-1][j] - 1) {
            i--;
        } else if (a[i][j] == a[i-1][j-1] + 1) {
            VERIFY(s1[i-1] == s2[j-1]);
            i--;
            j--;
        } else {
            VERIFY(a[i-1][j-1] - 1 == a[i][j] && s1[i-1] != s2[j-1]);
            i--;
            j--;
        }
    }
    return std::make_pair(size_t(answer), std::min(i_m - i, j_m - j));
}

inline const std::string ReverseComplement(const std::string &s) {
    std::string res(s.size(), 0);
    transform(s.begin(), s.end(), res.rbegin(), nucl_complement); // only difference with reverse is rbegin() instead of begin()
    return res;
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
        if (upper_position * 2 + 1 >= upper_length_)
            return (2 * upper_position + 1) * lower_length_
                    / (2 * upper_length_);
        else
            return lower_length_ - 1
                    - GetPosition(upper_length_ - 1 - upper_position);
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
        VERIFY(upper_position > 0);
        if (lower_length_ == 1)
            return 1;
        return (2 * upper_position * lower_length_ + upper_length_)
            / (2 * upper_length_);
    }
};

#endif /* SEQUENCE_TOOLS_HPP_ */
