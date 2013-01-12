#ifndef __HAMMER_HSEQ_HPP__
#define __HAMMER_HSEQ_HPP__

#include "sequence/nucl.hpp"

#include <array>

#include <cstdint>

namespace hammer {

union HomopolymerRun {
  uint8_t raw;
  struct {
    uint8_t len  : 6;
    uint8_t nucl : 2;
  };

  HomopolymerRun()
      : raw(0) {}
  HomopolymerRun(uint8_t nucl, uint8_t len)
      : len(len), nucl(nucl) {}

  bool operator==(const HomopolymerRun &that) const {
    return raw == that.raw;
  }

  std::string str() const {
    return std::string(len, ::nucl(nucl));
  }
};

template <size_t N = 16>
class HSeq {
  typedef std::array<HomopolymerRun, N> StorageType;
  StorageType data_;

  const static size_t PrimeNum = 239;

public:
  HSeq() {}

  HSeq(typename StorageType::const_iterator Start,
       typename StorageType::const_iterator End) {
    std::copy(Start, End, data_.begin());
  }

  typedef HomopolymerRun DataType;
  const static size_t DataSize = N;
  const static size_t TotalBytes = sizeof(DataType) * DataSize;

  static size_t GetDataSize(size_t size) {
    VERIFY(size == N);
    return N * sizeof(HomopolymerRun);
  }

  const HomopolymerRun *data() const {
    return data_.data();
  }

  size_t data_size() const {
    return DataSize;
  }

  HomopolymerRun &operator[](size_t idx) {
    return data_[idx];
  }

  const HomopolymerRun &operator[](size_t idx) const {
    return data_[idx];
  }

  HSeq<N> operator!() const {
    HSeq<N> res(*this);

    for (size_t i = 0; i < N / 2; ++i) {
      HomopolymerRun front = res[i], back = res[N - i - 1];
      front.nucl = complement(front.nucl);
      back.nucl = complement(back.nucl);
      res[i] = back;
      res[N - i + 1] = front;
    }

    if (N & 1)
      res[N/2].nucl = complement(res[N/2].nucl);

    return res;
  }

  HSeq<N> operator<<(char nucl) const {
    if (is_nucl(nucl))
      nucl = dignucl(nucl);

    HSeq<N> res(*this);
    // Easy case - just add to run
    HomopolymerRun &last = res[N-1];
    if (last.nucl == nucl) {
      last.len += 1;
      return res;
    }

    // Hard case - have to shift the stuff
    for (size_t i = 0; i < N - 1; ++i)
      res[i] = res[i + 1];
    res[N - 1].nucl = nucl;
    res[N - 1].len = 1;

    return res;
  }

  HSeq<N>& operator<<=(char nucl) {
    if (is_nucl(nucl))
      nucl = dignucl(nucl);

    // Easy case - just add to run
    HomopolymerRun &last = data_[N-1];
    if (last.nucl == nucl) {
      last.len += 1;
      return *this;
    }

    // Hard case - have to shift the stuff
    for (size_t i = 0; i < N - 1; ++i)
      data_[i] = data_[i + 1];
    data_[N - 1].nucl = nucl;
    data_[N - 1].len = 1;

    return *this;
  }

  HSeq<N> operator>>(char nucl) const {
    if (is_nucl(nucl))
      nucl = dignucl(nucl);

    HSeq<N> res(*this);
    // Easy case - just add to run
    HomopolymerRun &first = res[0];
    if (first.nucl == nucl) {
      first.len += 1;
      return res;
    }

    // Hard case - have to shift the stuff
    for (size_t i = 0; i < N - 1; ++i)
      res[i + 1] = res[i];
    res[0].nucl = nucl;
    res[0].len = 1;

    return res;
  }

  bool operator==(const HSeq<N> &that) const {
    return (data_ == that.data_);
  }
  bool operator!=(const HSeq<N> &that) const {
    return (data_ != that.data_);
  }

  size_t size() const {
    size_t res = 0;
    for (size_t i = 0; i < N; ++i)
      res += data_[i].len;

    return res;
  }

  std::string str() const {
    std::string res;
    for (size_t i = 0; i < N; ++i)
      res += data_[i].str();

    return res;
  }

  static size_t GetHash(const DataType *data, size_t sz = DataSize) {
    size_t hash = PrimeNum;

    for (size_t i = 0; i < sz; i++) {
      hash = ((hash << 5) - hash) + data[i].raw;
    }
    return hash;
  }

  size_t GetHash() const {
    return GetHash(data_.data());
  }

  struct hash {
    size_t operator()(const HSeq<N> &seq) const {
      return seq.GetHash();
    }

    size_t operator()(const DataType *data, size_t sz = DataSize) {
      return GetHash(data, sz);
    }
  };

  struct less2_fast {
    bool operator()(const HSeq<N> &l, const HSeq<N> &r) const {
      for (size_t i = 0; i < N; ++i) {
        const uint8_t lr = l[i].raw, rr = r[i].raw;
        if (lr != rr)
          return lr < rr;
      }

      return false;
    }
  };
};

template<size_t N>
std::ostream& operator<<(std::ostream& os, const HSeq<N> &seq) {
  os << seq.str();
  return os;
}

};

#endif // __HAMMER_HSEQ_HPP__
