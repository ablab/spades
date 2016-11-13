//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef CONCURRENTDSU_HPP_
#define CONCURRENTDSU_HPP_

#include "io/kmers/mmapped_writer.hpp"

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstdarg>
#include <cstdint>

#include <algorithm>
#include <vector>
#include <unordered_map>
#include <atomic>
#include <fstream>

// Silence bogus gcc warnings
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"

class ConcurrentDSU {
    struct atomic_set_t {
        uint64_t data  : 61;
        uint64_t aux   : 2;
        bool root  : 1;
    } __attribute__ ((packed));

    static_assert(sizeof(atomic_set_t) == 8, "Unexpected size of atomic_set_t");

public:
    ConcurrentDSU(size_t size)
            : data_(size) {

        for (size_t i = 0; i < size; i++)
            data_[i] = {.data = 1, .aux = 0, .root = true};
    }

    ~ConcurrentDSU() { }

    void unite(size_t x, size_t y) {
        uint64_t x_size, y_size;
        uint64_t x_aux, y_aux;

        // Step one: update the links
        while (true) {
            x = find_set(x);
            y = find_set(y);
            if (x == y)
                return;

            atomic_set_t x_entry = data_[x], y_entry = data_[y];
            // If someone already changed roots => retry
            if (!x_entry.root || !y_entry.root)
                continue;

            // We need to link the smallest subtree to the largest
            x_size = x_entry.data, y_size = y_entry.data;
            x_aux = x_entry.aux, y_aux = y_entry.aux;
            if (x_size > y_size || (x_size == y_size && x > y)) {
                std::swap(x, y);
                std::swap(x_size, y_size);
                std::swap(x_aux, y_aux);
                std::swap(x_entry, y_entry);
            }

            // Link 'x' to 'y'. If someone already changed 'x' => try again.
            atomic_set_t new_x_entry = {.data = y, .aux = x_aux, .root = false};
            if (!data_[x].compare_exchange_strong(x_entry, new_x_entry))
                continue;

            break;
        }

        // Step two: update the size.  We already linked 'x' to 'y'. Therefore we
        // need to add 'x_size' to whichever value is currently inside 'y'.
        while (true) {
            y = find_set(y);
            atomic_set_t y_entry = data_[y];
            // If someone already changed the roots => retry
            if (!y_entry.root)
                continue;

            // Update the size. If someone already changed 'y' => try again.
            atomic_set_t new_y_entry = {.data = x_size + y_entry.data, .aux = y_aux, .root = true};
            if (!data_[y].compare_exchange_strong(y_entry, new_y_entry))
                continue;

            break;
        }
    }

    size_t set_size(size_t i) const {
        while (true) {
            size_t el = find_set(i);
            atomic_set_t entry = data_[el];
            if (!entry.root)
                continue;

            return entry.data;
        }
    }

    size_t find_set(size_t x) const {
        // Step one: find the root
        size_t r = x;
        atomic_set_t r_entry = data_[r];
        while (!r_entry.root) {
            r = r_entry.data;
            r_entry = data_[r];
        }

        // Step two: traverse the path from 'x' to root trying to update the links
        // Note that the links might change, therefore we stop as soon as we'll
        // end at 'some' root.
        while (x != r) {
            atomic_set_t x_entry = data_[x];
            if (x_entry.root)
                break;

            // Try to update parent (may fail, it's ok)
            atomic_set_t new_x_entry = {.data = r, .aux = x_entry.aux, .root = false};
            data_[x].compare_exchange_weak(x_entry, new_x_entry);
            x = x_entry.data;
        }

        return x;
    }

    bool same(size_t x, size_t y) const {
        while (true) {
            x = find_set(x);
            y = find_set(y);
            if (x == y)
                return true;
            if (data_[x].load().root)
                return false;
        }
    }

    size_t num_sets() const {
        size_t count = 0;
        for (const auto &entry : data_) {
            count += entry.load(std::memory_order_relaxed).root;
        }

        return count;
    }

    bool is_root(size_t x) const {
        return data_[x].load(std::memory_order_relaxed).root;
    }

    uint64_t aux(size_t x) const {
        return data_[x].load(std::memory_order_relaxed).aux;
    }

    uint64_t root_aux(size_t x) const {
        while (true) {
            x = find_set(x);
            atomic_set_t entry = data_[x];

            if (!entry.root)
                continue;

            return entry.aux;
        }
    }

    void set_aux(size_t x, uint64_t data) {
        while (true) {
            atomic_set_t x_entry = data_[x];
            atomic_set_t new_x_entry = {.data = x_entry.data, .aux = data, .root = x_entry.root};
            if (!data_[x].compare_exchange_strong(x_entry, new_x_entry))
                continue;

            break;
        }
    }

    void set_root_aux(size_t x, uint64_t data) {
        while (true) {
            x = find_set(x);
            atomic_set_t x_entry = data_[x];
            if (!x_entry.root)
                continue;

            atomic_set_t new_x_entry = {.data = x_entry.data, .aux = data, .root = true};
            if (!data_[x].compare_exchange_strong(x_entry, new_x_entry))
                continue;

            break;
        }
    }

    size_t extract_to_file(const std::string &Prefix) {
        // First, touch all the sets to make them directly connect to the root
#   pragma omp parallel for
        for (size_t x = 0; x < data_.size(); ++x)
            (void) find_set(x);

        std::unordered_map<size_t, size_t> sizes;

#if 0
    for (size_t x = 0; x < size; ++x) {
        if (data_[x].parent != x) {
            size_t t = data_[x].parent;
            VERIFY(data_[t].parent == t)
        }
    }
#endif

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
        MMappedRecordWriter<size_t> index(Prefix + ".idx");
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

    void get_sets(std::vector<std::vector<size_t> > &otherWay) {
        otherWay.resize(data_.size());
        for (size_t i = 0; i < data_.size(); i++) {
            size_t set = find_set(i);
            otherWay[set].push_back(i);
        }
        otherWay.erase(remove_if(otherWay.begin(), otherWay.end(), zero_size),
                       otherWay.end());
    }

private:
    size_t parent(size_t x) const {
        atomic_set_t val = data_[x];
        return (val.root ? x : val.data);
    }

    static bool zero_size(const std::vector<size_t> &v) {
        return v.size() == 0;
    }

    mutable std::vector<std::atomic<atomic_set_t> > data_;
};

#pragma GCC diagnostic pop

#endif /* CONCURRENTDSU_HPP_ */
