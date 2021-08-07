#pragma once

#include "common_enums.hpp"

#include "utils/verify.hpp"

#include <vector>

namespace error_analyzer {

struct SeqFragment {
    unsigned long long start_pos;
    RangeType type;

    SeqFragment(unsigned long long from, RangeType type)
        : start_pos(from)
        , type(type)
    {}

    bool operator < (SeqFragment const & other) const noexcept {
        return start_pos < other.start_pos;
    }
};

struct SeqFragments : std::vector<SeqFragment> {
    size_t total_len = 0;

    size_t FindFragmentIndex(unsigned long long pos) const noexcept;

    std::pair<size_t, size_t> GetAllIntersectedFragments(unsigned long long start_pos, size_t len) const noexcept;
};

} // namespace error_analyzer
