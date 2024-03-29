
//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "lemiere_mod_reduce.hpp"

#include <functional>
#include <vector>
#include <atomic>

#include <cassert>

namespace bf {

namespace {
inline constexpr uint64_t cell_num(uint64_t x, uint64_t num) {
    return mod_reduce::multiply_high_u64(x, num);
}

inline constexpr uint64_t pairhash(uint64_t x, uint64_t y, uint64_t i) {
    return x + i * y + i * i;
}
}


/// The ordinary Bloom filter.
template<class T>
class bloom_filter {
    bloom_filter(const bloom_filter &) = delete;
    bloom_filter &operator=(const bloom_filter &) = delete;

protected:
    static constexpr size_t cells_per_entry_ = 8 * sizeof(uint64_t);

public:
    /// The hash digest type.F
    typedef uint64_t digest;

    /// The hash function type.
    typedef std::function<digest(const T &, uint64_t seed)> hasher;

    // FIXME disable default constructor
    bloom_filter() = default;

    ~bloom_filter() = default;

    /// Constructs a Bloom filter.
    /// @param h The hasher.
    /// @param cells The number of cells.
    /// @param num_hashes The number of hash functions to use
    /// The memory consumption will be cells bits
    bloom_filter(hasher h,
                 size_t cells, size_t num_hashes = 3)
            : hasher_(std::move(h)),
              num_hashes_(num_hashes),
              cells_(cells),
              data_((cells + cells_per_entry_ - 1) / cells_per_entry_) {}

    /// Move-constructs a Bloom filter.
    bloom_filter(bloom_filter &&) = default;

    /// Adds an element to the Bloom filter.
    /// @tparam T The type of the element to insert.
    /// @param x An instance of type `T`.
    bool add(const T &o) {
        bool dup = true;
        digest d1 = hasher_(o, 0xDEAD), d2 = hasher_(o, 0xBEEF);

        for (size_t i = 0; i < num_hashes_; ++i) {
            digest d = pairhash(d1, d2, i);
            size_t cell_id = cell_num(d, cells_);
            size_t pos = cell_id / cells_per_entry_;
            size_t epos = cell_id - pos * cells_per_entry_;

            uint64_t oldval = data_[pos].fetch_or(uint64_t(1) << epos);
            dup &= (oldval >> epos) & 1;
        }

        return !dup;
    }

    /// Retrieves the count of an element.
    /// @tparam T The type of the element to query.
    /// @param x An instance of type `T`.
    /// @return A frequency estimate for *x*.
    size_t lookup(const T &o) const {
        digest d1 = hasher_(o, 0xDEAD), d2 = hasher_(o, 0xBEEF);

        for (size_t i = 0; i < num_hashes_; ++i) {
            digest d = pairhash(d1, d2, i);
            size_t cell_id = cell_num(d, cells_);
            size_t pos = cell_id / cells_per_entry_;
            size_t epos = cell_id - pos * cells_per_entry_;
            if (((data_[pos] >> epos) & 1) == 0)
                return 0;
        }

        return 1;
    }

    /// Removes all items from the Bloom filter.
    void clear() {
        std::fill(data_.begin(), data_.end(), 0);
    }

    void merge(const bloom_filter<T> &other) {
        VERIFY(data_.size() == other.data_.size());
        VERIFY(num_hashes_ == other.num_hashes_);
        VERIFY(cells_ = other.cells_);

        for (size_t cell_id = 0; cell_id < cells_; ++cell_id) {
            size_t pos = cell_id / cells_per_entry_;
            data_[pos] |= other.data_[pos];
        }
    }

    template <typename Archive>
    void BinArchiveSave(Archive &ar) const {
        ar(num_hashes_, cells_, data_.size());
        ar.raw_array(data_.data(), data_.size());
    }

    template <typename Archive>
    void BinArchiveLoad(Archive &ar) {
        size_t size;
        ar(num_hashes_, cells_, size);
        if (data_.size() != size) {
            // data_.resize(size); // vector of atomics could not be resized
            data_ = std::vector<std::atomic<uint64_t>>(size);
        }
        ar.raw_array(data_.data(), data_.size());
    }

protected:
    hasher hasher_;
    size_t num_hashes_;
    size_t cells_;
    std::vector<std::atomic<uint64_t>> data_;
};

template<class T, size_t depth_>
class cascading_bloom_filter {
    cascading_bloom_filter(const cascading_bloom_filter &) = delete;
    cascading_bloom_filter &operator=(const cascading_bloom_filter &) = delete;

public:
    using hasher = typename bloom_filter<T>::hasher;
    using digets = typename bloom_filter<T>::digest;

    cascading_bloom_filter(hasher h,
                           size_t cells, size_t num_hashes = 3, double damp_factor = 0.1) {
        for (size_t i = 0; i < depth_; ++i) {
            filters_.emplace_back(h, cells, num_hashes);
            cells = size_t(double(cells) * damp_factor);
            if (cells < 1000)
                cells = 1000;
        }
    }

    size_t add(const T &o) {
        for (size_t i = 0; i < depth_; ++i) {
            if (filters_[i].lookup(o))
                continue;

            if (filters_[i].add(o))
                return i;
        }

        return depth_;
    }

    size_t lookup(const T &o) {
        for (size_t i = 0; i < depth_; ++i) {
            if (filters_[i].lookup(o))
                continue;

            return i;
        }

        return depth_;
    }

private:
    std::vector<bloom_filter<T>> filters_;
};

/// The counting Bloom filter.
template<class T, unsigned width_ = 4>
class counting_bloom_filter {
    counting_bloom_filter(const counting_bloom_filter &) = delete;
    counting_bloom_filter &operator=(const counting_bloom_filter &) = delete;

protected:
    static constexpr uint64_t cell_mask_ = (1ull << width_) - 1;
    static constexpr size_t cells_per_entry_ = 8 * sizeof(uint64_t) / width_;

public:
    /// The hash digest type.F
    typedef size_t digest;

    /// The hash function type.
    typedef std::function<digest(const T &, uint64_t seed)> hasher;

    // FIXME disable default constructor
    counting_bloom_filter() = default;

    ~counting_bloom_filter() = default;

    /// Constructs a counting Bloom filter.
    /// @param h The hasher.
    /// @param cells The number of cells.
    /// @param num_hashes The number of hash functions to use
    /// The memory consumption will be cells * width bits
    counting_bloom_filter(hasher h,
                          size_t cells, size_t num_hashes = 3)
            : hasher_(std::move(h)),
              num_hashes_(num_hashes),
              cells_(cells),
              data_((cells * width_ + 8 * sizeof(uint64_t) - 1) / 8 / sizeof(uint64_t)) {
        static_assert((width_ & (width_ - 1)) == 0, "Width must be power of two");
    }

    /// Move-constructs a counting Bloom filter.
    counting_bloom_filter(counting_bloom_filter &&) = default;

    /// Adds an element to the Bloom filter.
    /// @tparam T The type of the element to insert.
    /// @param x An instance of type `T`.
    void add(const T &o) {
        digest d1 = hasher_(o, 0xDEAD), d2 = hasher_(o, 0xBEEF);
        for (size_t i = 0; i < num_hashes_; ++i) {
            digest d = pairhash(d1, d2, i);
            size_t cell_id = cell_num(d, cells_);
            size_t pos = cell_id / cells_per_entry_;
            size_t epos = cell_id - pos * cells_per_entry_;
            auto &entry = data_[pos];
            uint64_t mask = cell_mask_ << (width_ * epos);

            // Add counter
            while (true) {
                uint64_t val = entry.load();

                // Overflow, do nothing
                if ((val & mask) == mask)
                    break;

                uint64_t newval = val + (1ull << (width_ * epos));
                if (!entry.compare_exchange_strong(val, newval))
                    continue;

                break;
            }

        }
    }

    /// Retrieves the count of an element.
    /// @tparam T The type of the element to query.
    /// @param x An instance of type `T`.
    /// @return A frequency estimate for *x*.
    size_t lookup(const T &o) const {
        size_t val = (1ull << width_) - 1;
        digest d1 = hasher_(o, 0xDEAD), d2 = hasher_(o, 0xBEEF);
        for (size_t i = 0; i < num_hashes_; ++i) {
            digest d = pairhash(d1, d2, i);
            size_t cell_id = cell_num(d, cells_);
            size_t pos = cell_id / cells_per_entry_;
            size_t epos = cell_id - pos * cells_per_entry_;
            size_t cval = (data_[pos] >> (width_ * epos)) & cell_mask_;
            if (val > cval)
                val = cval;
            if (val == 0)
                break;
        }

        return val;
    }

    /// Removes all items from the Bloom filter.
    void clear() {
        std::fill(data_.begin(), data_.end(), 0);
    }

    double load_factor() const {
      size_t loaded = 0;
      // FIXME: inefficient
      for (size_t cell_id = 0; cell_id < cells_; ++cell_id) {
          size_t pos = cell_id / cells_per_entry_;
          size_t epos = cell_id - pos * cells_per_entry_;
          const auto &entry = data_[pos];
          uint64_t mask = cell_mask_ << (width_ * epos);
          uint64_t val = entry & mask;
          loaded += val != 0;
      }
      return double(loaded) / double(cells_);
    }
    
    void merge(const counting_bloom_filter<T, width_> &other) {
        VERIFY(data_.size() == other.data_.size());
        VERIFY(num_hashes_ == other.num_hashes_);
        VERIFY(cells_ = other.cells_);

        for (size_t cell_id = 0; cell_id < cells_; ++cell_id) {
            size_t pos = cell_id / cells_per_entry_;
            size_t epos = cell_id - pos * cells_per_entry_;
            auto &entry = data_[pos];
            const auto &other_entry = other.data_[pos];
            uint64_t mask = cell_mask_ << (width_ * epos);
            uint64_t val = entry & mask;
            uint64_t other_val = other_entry & mask;
            uint64_t newval = (entry + other_val) & mask;
            if (newval < val) {
                newval = mask;
            }
            entry &= ~mask;
            entry |= newval;
        }
    }

    template <typename Archive>
    void BinArchiveSave(Archive &ar) const {
        ar(num_hashes_, cells_, data_.size());
        ar.raw_array(data_.data(), data_.size());
    }

    template <typename Archive>
    void BinArchiveLoad(Archive &ar) {
        size_t size;
        ar(num_hashes_, cells_, size);
        if (data_.size() != size) {
            // data_.resize(size); // vector of atomics could not be resized
            data_ = std::vector<std::atomic<uint64_t>>(size);
        }
        ar.raw_array(data_.data(), data_.size());
    }

protected:
    hasher hasher_;
    size_t num_hashes_;
    size_t cells_;
    std::vector<std::atomic<uint64_t>> data_;
};

/// The counting Bloom filter.
template<class T, unsigned width_ = 4>
class bitcounting_bloom_filter : public counting_bloom_filter<T, width_> {
    using typename counting_bloom_filter<T, width_>::digest;
    using typename counting_bloom_filter<T, width_>::hasher;

public:
    bitcounting_bloom_filter(hasher h,
                             size_t cells, size_t num_hashes = 3)
            : counting_bloom_filter<T, width_>(h, cells, num_hashes) { }

    /// Adds an element to the Bloom filter.
    /// @tparam T The type of the element to insert.
    /// @param x An instance of type `T`.
    void add(const T &o) {
        digest d1 = this->hasher_(o, 0xDEAD), d2 = this->hasher_(o, 0xBEEF);
        for (size_t i = 0; i < this->num_hashes_; ++i) {
            digest d = pairhash(d1, d2, i);
            size_t cell_id = cell_num(d, this->cells_);
            size_t pos = cell_id / this->cells_per_entry_;
            size_t epos = cell_id - pos * this->cells_per_entry_;
            auto &entry = this->data_[pos];
            uint64_t mask = this->cell_mask_ << (width_ * epos);

            // Add counter
            while (true) {
                uint64_t val = entry.load() & mask;

                // Overflow, do nothing
                if (val == mask)
                    break;

                uint64_t cellval = val >> width_ * epos;
                size_t cnt = (cellval == 0 ? 0 : 64 - __builtin_clzll(cellval)) + width_ * epos;

                if ((std::atomic_fetch_or(&entry, uint64_t(1) << cnt) & mask) != val)
                    continue;

                break;
            }
        }
    }

    /// Retrieves the count of an element.
    /// @tparam T The type of the element to query.
    /// @param x An instance of type `T`.
    /// @return A frequency estimate for *x*.
    size_t lookup(const T &o) const {
        size_t val = (1ull << width_) - 1;
        digest d1 = this->hasher_(o, 0xDEAD), d2 = this->hasher_(o, 0xBEEF);
        for (size_t i = 0; i < this->num_hashes_; ++i) {
            digest d = pairhash(d1, d2, i);
            size_t cell_id = cell_num(d, this->cells_);
            size_t pos = cell_id / this->cells_per_entry_;
            size_t epos = cell_id - pos * this->cells_per_entry_;
            uint64_t entry = (this->data_[pos] >> (width_ * epos)) & this->cell_mask_;
            size_t cval = (entry == 0 ? 0 : 64 - __builtin_clzll(entry));

            if (val > cval)
                val = cval;
        }

        return val;
    }

    void merge(...) {
        VERIFY_MSG(false, "Not implemented");
    }
};


} // namespace bf
