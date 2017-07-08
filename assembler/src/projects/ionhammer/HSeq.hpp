//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef __HAMMER_HSEQ_HPP__
#define __HAMMER_HSEQ_HPP__

#include <city/city.h>
#include "sequence/nucl.hpp"

#include <array>
#include <deque>
#include <string>
#include <vector>

#include <cstdint>

namespace hammer {

union HomopolymerRun {
  uint8_t raw;
  struct {
    uint8_t len : 6;
    uint8_t nucl : 2;
  };

  HomopolymerRun() : raw(0) {}
  HomopolymerRun(uint8_t nucl, uint8_t len = 1)
      : len(len & 63), nucl(nucl & 3) {}

  bool operator==(const HomopolymerRun &that) const { return raw == that.raw; }

  bool operator!=(const HomopolymerRun &that) const { return raw != that.raw; }

  bool operator<(const HomopolymerRun &that) const { return raw < that.raw; }

  inline char Nucl() const { return nucl; }

  inline char Len() const { return len; }

  std::string str() const { return std::string(len, ::nucl(nucl)); }
};

namespace iontorrent {
// Container shall have push_back method
template <typename Container>
void toHomopolymerRuns(const std::string &seq, Container &runs) {
  if (seq.empty()) return;

  char nucl = seq[0];
  uint8_t len = 1;
  for (size_t i = 1; i < seq.size(); ++i) {
    if (seq[i] != nucl) {
      runs.push_back(HomopolymerRun(dignucl(nucl), len));
      len = 1;
      nucl = seq[i];
    } else {
      ++len;
    }
  }
  if (len > 0) {
    runs.push_back(HomopolymerRun(dignucl(nucl), len));
  }
}

};  // namespace iontorrent

template <size_t N = 16>
class HSeq {
 public:
  typedef std::array<HomopolymerRun, N> StorageType;

 private:
  StorageType data_;

 public:
  HSeq() {}

  template <class Iterator>
  HSeq(Iterator start, Iterator end) {
    std::copy(start, end, data_.begin());
  }

  typedef HomopolymerRun DataType;
  const static size_t DataSize = N;
  const static size_t TotalBytes = sizeof(DataType) * DataSize;

  static size_t GetDataSize(size_t size) {
    VERIFY(size == N);
    return N * sizeof(HomopolymerRun);
  }

  typename StorageType::const_iterator begin() const { return data_.begin(); }

  typename StorageType::const_iterator end() const { return data_.end(); }

  typename StorageType::const_reverse_iterator rbegin() const {
    return data_.rbegin();
  }

  typename StorageType::const_reverse_iterator rend() const {
    return data_.rend();
  }

  const HomopolymerRun *data() const { return data_.data(); }

  size_t data_size() const { return DataSize; }

  HomopolymerRun &operator[](size_t idx) { return data_[idx]; }

  const HomopolymerRun &operator[](size_t idx) const { return data_[idx]; }

  HSeq<N> operator!() const {
    HSeq<N> res(*this);

    for (size_t i = 0; i < N / 2; ++i) {
      HomopolymerRun front = res[i], back = res[N - i - 1];
      front.nucl = complement(front.nucl) & 3;
      back.nucl = complement(back.nucl) & 3;
      res[i] = back;
      res[N - i - 1] = front;
    }

    if (N & 1) res[N / 2].nucl = complement(res[N / 2].nucl) & 3;

    return res;
  }

  HSeq<N> operator<<(char nucl) const {
    if (is_nucl(nucl)) nucl = dignucl(nucl);

    HSeq<N> res(*this);
    // Easy case - just add to run
    HomopolymerRun &last = res[N - 1];
    if (last.nucl == nucl) {
      last.len += 1;
      return res;
    }

    // Hard case - have to shift the stuff
    for (size_t i = 0; i < N - 1; ++i) res[i] = res[i + 1];
    res[N - 1].nucl = nucl;
    res[N - 1].len = 1;

    return res;
  }

  HSeq<N> operator<<(HomopolymerRun run) const {
    HSeq<N> res(*this);
    // Easy case - just add to run
    HomopolymerRun &last = res[N - 1];
    if (last.nucl == run.nucl) {
      last.len += run.len;
      return res;
    }

    // Hard case - have to shift the stuff
    for (size_t i = 0; i < N - 1; ++i) res[i] = res[i + 1];
    res[N - 1] = run;

    return res;
  }

  HSeq<N> &operator<<=(char nucl) {
    if (is_nucl(nucl)) nucl = dignucl(nucl);

    // Easy case - just add to run
    HomopolymerRun &last = data_[N - 1];
    if (last.nucl == nucl) {
      last.len = (last.len + 1) & 63;
      return *this;
    }

    // Hard case - have to shift the stuff
    for (size_t i = 0; i < N - 1; ++i) data_[i] = data_[i + 1];
    data_[N - 1].nucl = nucl & 3;
    data_[N - 1].len = 1;

    return *this;
  }

  HSeq<N> &operator<<=(HomopolymerRun run) {
    // Easy case - just add to run
    HomopolymerRun &last = data_[N - 1];
    if (last.nucl == run.nucl) {
      last.len = (last.len + run.len) & 63;
      return *this;
    }

    // Hard case - have to shift the stuff
    for (size_t i = 0; i < N - 1; ++i) data_[i] = data_[i + 1];
    data_[N - 1] = run;
    return *this;
  }

  HSeq<N> operator>>(HomopolymerRun run) const {
    HSeq<N> res(*this);
    // Easy case - just add to run
    HomopolymerRun &first = res[0];
    if (first.nucl == run.nucl) {
      first.len += run.len;
      return res;
    }

    // Hard case - have to shift the stuff
    for (size_t i = 0; i < N - 1; ++i) res[i + 1] = res[i];
    res[0].nucl = run.nucl;
    res[0].len = run.len;

    return res;
  }

  HSeq<N> operator>>(char nucl) const {
    if (is_nucl(nucl)) nucl = dignucl(nucl);

    HSeq<N> res(*this);
    // Easy case - just add to run
    HomopolymerRun &first = res[0];
    if (first.nucl == nucl) {
      first.len += 1;
      return res;
    }

    // Hard case - have to shift the stuff
    for (size_t i = 0; i < N - 1; ++i) res[i + 1] = res[i];
    res[0].nucl = nucl;
    res[0].len = 1;

    return res;
  }

  bool operator==(const HSeq<N> &that) const { return (data_ == that.data_); }
  bool operator!=(const HSeq<N> &that) const { return (data_ != that.data_); }

  size_t size() const {
    size_t res = 0;
    for (size_t i = 0; i < N; ++i) res += data_[i].len;

    return res;
  }
  size_t max_run_length() const {
    size_t res = 0;
    for (size_t i = 0; i < N; ++i) res = std::max((size_t)(data_[i].len), res);

    return res;
  }

  std::string str() const {
    std::string res;
    for (size_t i = 0; i < N; ++i) res += data_[i].str();

    return res;
  }

  static size_t GetHash(const DataType *data, size_t sz = DataSize,
                        uint32_t seed = 0) {
    return CityHash64WithSeed((const char *)data, sz * sizeof(DataType),
                              0x9E3779B9 ^ seed);
  }

  size_t GetHash(uint32_t seed = 0) const {
    return GetHash(data_.data(), DataSize, seed);
  }

  struct hash {
    size_t operator()(const HSeq<N> &seq, uint32_t seed = 0) const {
      return seq.GetHash(seed);
    }

    size_t operator()(const DataType *data, size_t sz = DataSize,
                      uint32_t seed = 0) const {
      return GetHash(data, sz, seed);
    }
  };

  struct less2_fast {
    bool operator()(const HSeq<N> &l, const HSeq<N> &r) const {
      for (size_t i = 0; i < N; ++i) {
        const uint8_t lr = l[i].raw, rr = r[i].raw;
        if (lr != rr) return lr < rr;
      }

      return false;
    }
  };
};

template <size_t N>
std::ostream &operator<<(std::ostream &os, const HSeq<N> &seq) {
  os << seq.str();
  return os;
}

namespace internal {
template <size_t N>
inline size_t getSize(const hammer::HSeq<N> &) {
  return N;
}

template <typename T>
inline size_t getSize(const T &a) {
  return a.size();
}
}  // namespace internal

};  // namespace hammer

#endif  // __HAMMER_HSEQ_HPP__
