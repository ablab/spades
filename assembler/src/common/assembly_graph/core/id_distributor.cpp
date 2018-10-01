#include "id_distributor.hpp"

using namespace omnigraph;

uint64_t ReclaimingIdDistributor::round_up_stride(uint64_t n) const {
    return uint64_t(n + stride_- 1 ) / stride_d_ * stride_;
}

    
uint64_t ReclaimingIdDistributor::next_free(uint64_t n, uint64_t stride) const {
    for (size_t i = n; i < free_map_.size(); i+= stride) {
        if (free_map_[i])
            return i;
    }

    return free_map_.size();
}

void ReclaimingIdDistributor::resize(size_t sz) {
    fprintf(stderr, "!!!RESIZE!!!!\n");
    sz = round_up_stride(sz);
    free_map_.resize(sz, true);
}

uint64_t ReclaimingIdDistributor::allocate(uint64_t offset) {
    // First hint: see if we could find any spot after last allocated
    uint64_t hint = round_up_stride(last_allocated_) + offset;
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

void ReclaimingIdDistributor::release(uint64_t id) {
    free_map_[id - bias_] = true;
}

