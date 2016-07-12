#pragma once

#include <functional>
#include <vector>
#include <atomic>

#include <cassert>

namespace bf {

/// The counting Bloom filter.
template<class T, unsigned width_ = 4>
class counting_bloom_filter {
  counting_bloom_filter(const counting_bloom_filter&) = delete;
  counting_bloom_filter& operator=(const counting_bloom_filter&) = delete;

  static constexpr uint64_t cell_mask_ = (1ull << width_) - 1;
  static constexpr size_t cells_per_entry_ = 8 * sizeof(uint64_t) / width_;
    
public:
  /// The hash digest type.
  typedef size_t digest;

  /// The hash function type.
  typedef std::function<digest(const T&, uint64_t seed)> hasher;

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
      data_((cells * width_ + 8 * sizeof(uint64_t) - 1)/ 8 / sizeof(uint64_t)) {
    static_assert((width_ & (width_ - 1)) == 0, "Width must be power of two");
  }

  /// Move-constructs a counting Bloom filter.
  counting_bloom_filter(counting_bloom_filter&&) = default;

  /// Adds an element to the Bloom filter.
  /// @tparam T The type of the element to insert.
  /// @param x An instance of type `T`.
  void add(const T &o) {
    for (size_t i = 0; i < num_hashes_; ++i) {
      digest d = hasher_(o, i);
      size_t cell_id = d - cells_ * (d / cells_);
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
    for (size_t i = 0; i < num_hashes_; ++i) {
      digest d = hasher_(o, i);
      size_t cell_id = d - cells_ * (d / cells_);
      size_t pos = cell_id / cells_per_entry_;
      size_t epos = cell_id - pos * cells_per_entry_;
      size_t cval = (data_[pos] >> (width_ * epos)) & cell_mask_;
      if (val > cval)
        val = cval;
    }

    return val;    
  }
  
  /// Removes all items from the Bloom filter.
  void clear() {
    std::fill(data_.begin(), data_.end(), 0);
  }
  

private:
  hasher hasher_;
  size_t num_hashes_;
  size_t cells_;
  std::vector<std::atomic<uint64_t>> data_;
};
  
} // namespace bf
