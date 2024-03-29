
//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <vector>
#include <functional>
#include <numeric>
#include <cmath>

namespace hll {

template<unsigned precision = 24>
class hll {
    static constexpr uint64_t m_ = 1ull << precision;
    static constexpr uint64_t mask_ = (m_ - 1) << (64 - precision);

    constexpr double alpha(unsigned p) const {
      // constexpr switches are C++14 only :(
      return (p > 6 ?
              0.7213 / (1.0 + 1.079 / double(1ull << p)) :
              p == 6 ? 0.709 : p == 5 ? 0.697 : 0.673);
    }

public:
    /// The hash digest type.
    typedef uint64_t digest;

    hll()
      : data_(1ull << precision, 0) { }

    void add(digest d) {
      // Split digest into parts
      size_t id = (d & mask_) >> (64 - precision);
      uint8_t rho = uint8_t(((d & ~mask_) == 0 ? 64 : __builtin_clzll(d & ~mask_)) - precision + 1);
      if (data_[id] < rho)
        data_[id] = rho;
    }

    void merge(const hll &other) {
      VERIFY(data_.size() == other.data_.size());
      for (size_t i = 0; i < data_.size(); ++i)
        data_[i] = std::max(data_[i], other.data_[i]);
    }

    double cardinality() const {
      // FIXME: Precision loss?
      // FIXME: Bias correction!
      double res = alpha(precision) * m_ * m_;
      double E = std::accumulate(data_.begin(), data_.end(),
                                 0.0, [](double a, uint8_t b) { return a + exp2(-(double) b); });
      uint64_t zeros_bucket_cnt = std::accumulate(data_.begin(), data_.end(), 0,
          [](uint64_t a, uint8_t b) {return a + (b == 0);});
      res /= E;
      if (res <= 5.0 * m_/2 && zeros_bucket_cnt > 0) {
          return m_ * (std::log((double)m_) - std::log((double)zeros_bucket_cnt));
      } else {
          return res;
      }
    }

    double upper_bound_cardinality() const {
        return 1.1 * cardinality();
    }

    void clear() {
      std::fill(data_.begin(), data_.end(), 0);
    }

    template <typename Archive>
    void BinArchive(Archive &ar) {
        ar(data_);
    }

private:
    std::vector<uint8_t> data_;
};

template<class T, unsigned precision = 24>
class hll_with_hasher : public hll<precision> {
public:
    using typename hll<precision>::digest;
    using hll<precision>::add;

    /// The hash function type.
    typedef std::function<digest(const T)> hasher;

    hll_with_hasher(hasher h = nullptr)
            : hasher_(std::move(h)) { }

    /// @tparam T The type of the element to insert.
    /// @param o An instance of type `T`.
    void add(const T &o) {
      digest d = hasher_(o);
      add(d);
    }

  private:
    hasher hasher_;
};

} // hll
