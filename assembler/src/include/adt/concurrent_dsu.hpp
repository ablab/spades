//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef CONCURRENTDSU_HPP_
#define CONCURRENTDSU_HPP_

#include "io/mmapped_writer.hpp"

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstdarg>
#include <cstdint>

#include <algorithm>
#include <vector>
#include <unordered_map>

class ConcurrentDSU {
  union atomic_set_t {
    uint64_t raw;
    struct {
      uint32_t next  : 32;
      uint32_t size  : 31;
      uint32_t dirty : 1;
    };
  } __attribute__ ((packed));

 public:
  ConcurrentDSU(size_t size) :
      size(size) {
    if (size > 0xFFFFFFFFULL) {
      std::cerr << "Error, size greater than 2^32 -1";
      exit(-1);
    }

    data = new atomic_set_t[size];
    for (size_t i = 0; i < size; i++) {
      data[i].next = (uint32_t) i;
      data[i].size = 1;
      data[i].dirty = 0;
    }
  }

  void unite(unsigned x, unsigned y) {
    while (true) {
      x = find_set(x);
      y = find_set(y);
      if (x == y) {
        return;
      }

      unsigned x_size = data[x].size;
      unsigned y_size = data[y].size;
      if (x_size > y_size || (x_size == y_size && x > y)) {
        std::swap(x, y);
        std::swap(x_size, y_size);
      }

      if (lock(y)) {
        if (update_root(x, x_size, y, x_size)) {
          while (data[x].dirty) {
            // SPIN LOCK
          }
          while (true) {
            atomic_set_t old = data[y];
            atomic_set_t nnew = old;
            nnew.size = (uint32_t) ((nnew.size + data[x].size) & 0x7fffffff);
            if (__sync_bool_compare_and_swap(&data[y].raw, old.raw, nnew.raw)) {
              break;
            }
          }
        }
        unlock(y);
      }
    }
  }

  size_t set_size(unsigned i) const {
    unsigned el = find_set(i);
    return data[el].size;
  }

  unsigned find_set(unsigned x) const {
    unsigned ans = x;
    while (ans != data[ans].next) {
      unsigned t = data[ans].next;
      atomic_set_t old = data[ans];
      atomic_set_t nnew = old;
      nnew.next = data[t].next;
      __sync_bool_compare_and_swap(&data[ans].raw, old.raw, nnew.raw);
      ans = data[t].next;
    }

    return ans;
  }

  size_t num_sets() const {
    size_t count = 0;
    for (size_t i = 0; i < size; i++) {
      if (data[i].next == i)
        count++;
    }
    return count;
  }

    void get_sets2(const std::string& FileName) {
    // First, touch all the sets to make them directly connect to the root
#   pragma omp parallel for
    for (size_t x = 0; x < size; ++x)
        (void)find_set(x);

    std::unordered_map<size_t, size_t> sizes;

    // Insert all the root elements into the map
    for (size_t x = 0; x < size; ++x) {
        if (data[x].next == x)
            sizes[x] = 0;
    }

    // Now, calculate the counts. We can do this in parallel, because we know no
    // insertion can occur.
#   pragma omp parallel for
    for (size_t x = 0; x < size; ++x) {
        size_t& entry = sizes[data[x].next];
#       pragma omp atomic
        entry += 1;
    }

    // Now we know the sizes of each cluster. Go over again and calculate the
    // file-relative (cumulative) offsets.
    size_t off = 0;
    for (size_t x = 0; x < size; ++x) {
        if (data[x].next == x) {
            size_t& entry = sizes[x];
            size_t noff = off + entry;
            entry = off;
            off = noff;
        }
    }

    // Write down the entries
    MMappedRecordWriter<size_t> os(FileName);
    os.reserve(off);
    size_t *out = os.data();
    for (size_t x = 0; x < size; ++x) {
        size_t& entry = sizes[data[x].next];
        out[entry++] = x;
    }

    // Write down the sizes
    MMappedRecordWriter<size_t> index(FileName + ".idx");
    index.reserve(sizes.size());
    size_t *idx = index.data();
    for (size_t x = 0, i = 0, sz = 0; x < size; ++x) {
        if (data[x].next == x) {
            idx[i++] = sizes[x] - sz;
            sz = sizes[x];
        }
    }

  }

  void get_sets(std::vector<std::vector<size_t> > &otherWay) {
    otherWay.resize(size);
    for (size_t i = 0; i < size; i++) {
      unsigned set = find_set(i);
      otherWay[set].push_back(i);
    }
    otherWay.erase(remove_if(otherWay.begin(), otherWay.end(), zero_size),
                   otherWay.end());
  }

private:
  bool lock(unsigned y) {
    while (true) {
      atomic_set_t old = data[y];
      if (old.next != y) {
        return false;
      }
      old.dirty = 0;
      atomic_set_t nnew = old;
      nnew.dirty = 1;
      if (__sync_bool_compare_and_swap(&data[y].raw, old.raw, nnew.raw)) {
        return true;
      }
    }
    return false;
  }

  void unlock(int y) {
    data[y].dirty = 0;
  }

  static bool zero_size(const std::vector<size_t> & v) {
    return v.size() == 0;
  }

  bool update_root(uint32_t x, uint32_t oldrank, uint32_t y,
                   uint32_t newrank) {
    atomic_set_t old = data[x];
    if (old.next != x || old.size != oldrank) {
      return false;
    }
    atomic_set_t nnew = old;
    nnew.next = y;
    nnew.size = newrank & 0x7fffffff;
    return __sync_bool_compare_and_swap(&data[x].raw, old.raw, nnew.raw);
  }

  mutable atomic_set_t *data;
  size_t size;
};

#endif /* CONCURRENTDSU_HPP_ */
