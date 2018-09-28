//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <vector>
#include <cstddef>

namespace omnigraph {

class ReclaimingIdDistributor {
  public:
    ReclaimingIdDistributor(uint64_t bias = 0, size_t initial_size = 1, uint64_t stride = 1)
            : bias_(bias), stride_(stride), stride_d_(stride) {
        resize(initial_size);
    }

    void resize(size_t sz) {
        sz = round_up_stride(sz);
        free_map_.resize(sz, true);
    }

    uint64_t allocate(uint64_t offset = 0);
    void release(uint64_t id);

    size_t size() const {
        return free_map_.size();
    }

    size_t free() const {
        size_t res = 0;
        for (bool flag : free_map_)
            res += flag;
        return res;
    }


  private:
    uint64_t round_up_stride(uint64_t n) const;    
    uint64_t next_free(uint64_t n, uint64_t stride) const;

    uint64_t next_free(uint64_t n = 0) const {
        return next_free(n, stride_);
    }

    uint64_t last_allocated_;

    uint64_t bias_;
    uint64_t stride_;
    uint64_t stride_d_;
    std::vector<bool> free_map_;
};

}
