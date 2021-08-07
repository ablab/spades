#include "sequence_fragments.hpp"

namespace error_analyzer {

size_t SeqFragments::FindFragmentIndex(unsigned long long pos) const noexcept {
    size_t l = 0, r = size();
    while (l + 1 < r) {
        auto mid = (l + r) / 2;
        if (pos < (*this)[mid].start_pos)
            r = mid;
        else
            l = mid;
    }
    return l;
}

std::pair<size_t, size_t> SeqFragments::GetAllIntersectedFragments(unsigned long long start_pos, size_t len) const noexcept {
    auto first_fragment = FindFragmentIndex(start_pos);
    auto end_pos = start_pos + len;
    VERIFY(end_pos <= total_len);
    auto end_fragment = first_fragment + 1;
    while (end_fragment < size() && (*this)[end_fragment].start_pos < end_pos)
        ++end_fragment;
    return {first_fragment, end_fragment - 1};
}

} // namespace error_analyzer
