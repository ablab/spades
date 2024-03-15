//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "concurrent_dsu.hpp"

#include "io/kmers/mmapped_writer.hpp"
#include "utils/logger/logger.hpp"

#include <parallel_hashmap/phmap.h>
#include <fstream>

namespace dsu {

size_t ConcurrentDSU::extract_to_file(const std::string &Prefix) {
    INFO("Connecting to root");
    // First, touch all the sets to make them directly connect to the root
#   pragma omp parallel for
    for (size_t x = 0; x < data_.size(); ++x)
        (void) find_set(x);

    phmap::flat_hash_map<size_t, size_t> sizes;
#if 0
    for (size_t x = 0; x < size; ++x) {
        if (data_[x].parent != x) {
            size_t t = data_[x].parent;
            VERIFY(data_[t].parent == t)
        }
    }
#endif

    INFO("Calculating counts");
    // Insert all the root elements into the map
    sizes.reserve(num_sets());
    for (size_t x = 0; x < data_.size(); ++x) {
        if (is_root(x))
            sizes[x] = 0;
    }

    // Now, calculate the counts. We can do this in parallel, because we know no
    // insertion can occur.
#   pragma omp parallel for
    for (size_t x = 0; x < data_.size(); ++x) {
        size_t &entry = sizes[parent(x)];
#       pragma omp atomic
        entry += 1;
    }

    // Now we know the sizes of each cluster. Go over again and calculate the
    // file-relative (cumulative) offsets.
    size_t off = 0;
    for (size_t x = 0; x < data_.size(); ++x) {
        if (is_root(x)) {
            size_t &entry = sizes[x];
            size_t noff = off + entry;
            entry = off;
            off = noff;
        }
    }

    INFO("Writing down entries");
    // Write down the entries
    std::vector<size_t> out(off);
    for (size_t x = 0; x < data_.size(); ++x) {
        size_t &entry = sizes[parent(x)];
        out[entry++] = x;
    }
    std::ofstream os(Prefix, std::ios::binary | std::ios::out);
    os.write((char *) &out[0], out.size() * sizeof(out[0]));
    os.close();

    // Write down the sizes
    MMappedRecordWriter <size_t> index(Prefix + ".idx");
    index.reserve(sizes.size());
    size_t *idx = index.data();
    for (size_t x = 0, i = 0, sz = 0; x < data_.size(); ++x) {
        if (is_root(x)) {
            idx[i++] = sizes[x] - sz;
            sz = sizes[x];
        }
    }

    return sizes.size();
}

}
