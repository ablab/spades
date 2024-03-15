#include "id_distributor.hpp"

using namespace omnigraph;

uint64_t ReclaimingIdDistributor::next_free(uint64_t n) const {
    for (size_t i = n; i < free_map_.size(); ++i) {
        if (free_map_[i])
            return i;
    }

    return free_map_.size();
}

void ReclaimingIdDistributor::resize(size_t sz) {
    //fprintf(stderr, "!!!RESIZE!!!! %llu\n", sz);
    free_map_.resize(sz, true);
}

uint64_t ReclaimingIdDistributor::allocate(uint64_t offset) {
    // First hint: see if we could find any spot after last allocated
    uint64_t hint = last_allocated_ + offset;
    uint64_t n = next_free(hint);
    if (n == free_map_.size()) {
        // No luck, start from the beginning
        n = next_free();
    }

    // Still no luck, resize
    if (n == free_map_.size())
        resize(free_map_.size() * 2);

    last_allocated_ = n;
    free_map_[n] = false;
    return n + bias_;
}

size_t ReclaimingIdDistributor::free() const {
    size_t res = 0;
    for (bool flag : free_map_)
        res += flag;
    return res;
}
