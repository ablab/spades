//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <cstdlib>
#include <cstdint>
#include <algorithm>

namespace cap {

namespace precalc {

template <uint64_t base, uint64_t degree, class Enable = void>
struct Power {
};
template <uint64_t base, uint64_t degree>
struct Power<base, degree, typename boost::enable_if_c<degree % 2 == 1>::type > {
  enum { value = Power<base * base, degree / 2>::value * base };
};
template <uint64_t base, uint64_t degree>
struct Power<base, degree, typename boost::enable_if_c<degree % 2 == 0>::type > {
  enum { value = Power<base * base, degree / 2>::value };
};

template <uint64_t base>
struct Power<base, 0, void> {
  enum { value = 1 };
};

template <uint64_t base>
struct Power<base, 1, void> {
  enum { value = base };
};

}

namespace utils {

template <class T>
inline T FastPow(T base, size_t pow) {
  if (pow == 0) {
    return 1;
  }
  if (pow & 1) {
    return FastPow(base * base, pow / 2) * base;
  }
  return FastPow(base * base, pow / 2);
}

template <size_t ind>
struct PrimeHolder {
};
template <>
struct PrimeHolder<0> {
  enum { val = 239 };
};
template <>
struct PrimeHolder<1> {
  enum { val = 241 };
};
template <>
struct PrimeHolder<2> {
  enum { val = 251 };
};
template <>
struct PrimeHolder<3> {
  enum { val = 257 };
};
template <>
struct PrimeHolder<4> {
  enum { val = 263 };
};
template <>
struct PrimeHolder<5> {
  enum { val = 269 };
};
/*
template <size_t ind>
class PrimeHolder {
  enum { val = ({239, 241, 251, 257, 263, 269, 271, 277, 281, 283})[ind] };
};
*/

}

template <size_t prime = 239, class HashT = uint64_t>
class PolynomialHash : private boost::noncopyable {
  static const size_t kHashP = prime;
  HashT kHashPDeg;
  HashT hash_;
//  uint8_t last_chars_; // assume every char <= 2bits (preferably)

  inline static HashT GenPolyDeg(unsigned polydeg) {
    return utils::FastPow<HashT>(prime, polydeg - 1);
  }

  explicit PolynomialHash()
    : kHashPDeg(0),
      hash_(0) {
  }
  
 public:
  typedef HashT DataType;

  PolynomialHash(unsigned polydeg)
      : kHashPDeg(GenPolyDeg(polydeg)),
        hash_(0)/*,
        last_chars_(0)*/ {
  }

  PolynomialHash(unsigned polydeg, const PolynomialHash<prime, HashT> &other)
      : kHashPDeg(GenPolyDeg(polydeg)),
        hash_(other.hash_)/*,
        last_chars_(0)*/ {
  }

  inline void CopyFrom(unsigned polydeg, const PolynomialHash<prime, HashT> &other) {
    kHashPDeg = GenPolyDeg(polydeg);
    hash_ = other.hash_;
    //last_chars_ = other.last_chars_;
  }


  inline void Update(HashT front) {
    hash_ = hash_ * kHashP + front;
    //last_chars_ = (last_chars_ << 2) ^ (front & 8);
  }
  inline void Update(HashT front, HashT back) {
    hash_ = (hash_ - back * kHashPDeg) * kHashP + front;
    //last_chars_ = (last_chars_ << 2) ^ (front & 8);
  }
  inline void UpdateBack(HashT back) {
    hash_ = hash_ + kHashPDeg * back;
    //last_chars_ >>= 2;
  }
  inline void UpdateBack(HashT back, HashT front) {
    hash_ = (hash_ - front) * precalc::Power<kHashP, (1ull << 63) - 1>::value + back * kHashPDeg;
    //last_chars_ >>= 2;
  }
  inline HashT GetHash() const {
    return hash_;
  }
  inline HashT GetHash(HashT /*seed*/) const {
    return hash_;// + last_chars_ * seed * kHashPDeg * kHashP;
  }

  inline bool operator ==(const PolynomialHash<prime, HashT> &other) const {
    return hash_ == other.hash_;
  }
  inline bool operator !=(const PolynomialHash<prime, HashT> &other) const {
    return !(operator==(other));
  }
};

template <size_t size, class StorageT>
class HashTuple {
  typedef HashTuple<size - 1, StorageT> ChildClass;

  ChildClass child_data_;
  StorageT data_;

 public:
  HashTuple() : child_data_(), data_() {
  }
  HashTuple(ChildClass child_data, StorageT data)
      : child_data_(child_data),
        data_(data) {
  }
  template <size_t pos>
  inline void set(StorageT data, typename boost::enable_if_c<pos == size - 1>::type* = 0) {
    data_ = data;
  }
  template <size_t pos>
  inline void set(StorageT data, typename boost::enable_if_c<pos != size - 1>::type* = 0) {
    child_data_.template set<pos>(data);
  }

  template <size_t pos>
  inline StorageT get(typename boost::enable_if_c<pos == size - 1>::type* = 0) {
    return data_;
  }
  template <size_t pos>
  inline StorageT get(typename boost::enable_if_c<pos != size - 1>::type* = 0) {
    return child_data_.template get<pos>();
  }
};

template <class StorageT>
class HashTuple<1, StorageT> {
  StorageT data_;

 public:
  HashTuple() : data_() {
  }
  HashTuple(StorageT data) : data_(data) {
  }

  template <size_t pos>
  inline void set(StorageT data) {
    data_ = data;
  }
  template <size_t pos>
  inline StorageT get() const {
    return data_;
  }
};

// nobody knows how to do it in a good way
template <size_t size, class HashT = uint64_t>
class MultiPolynomialHash {
  typedef MultiPolynomialHash<size - 1, HashT> ChildClass;

  ChildClass child_hash_;
  PolynomialHash<utils::PrimeHolder<size - 1>::val, HashT> hash_;

 public:
  typedef HashTuple<size, HashT> DataType;
  typedef HashT AtomType;

  MultiPolynomialHash(unsigned polydeg) : child_hash_(polydeg), hash_(polydeg) {
  }

  MultiPolynomialHash(unsigned polydeg, const MultiPolynomialHash<size, HashT> &other)
      : child_hash_(polydeg, other.child_hash_),
        hash_(polydeg, other.hash_) {
  }
  inline void CopyFrom(unsigned polydeg, const MultiPolynomialHash<size, HashT> &other) {
    child_hash_.CopyFrom(polydeg, other.child_hash_);
    hash_.CopyFrom(polydeg, other.hash_);
  };

  inline void Update(HashT front) {
    child_hash_.Update(front);
    hash_.Update(front);
  }
  inline void Update(HashT front, HashT back) {
    child_hash_.Update(front, back);
    hash_.Update(front, back);
  }
  inline void UpdateBack(HashT back) {
    child_hash_.UpdateBack(back);
    hash_.UpdateBack(back);
  }
  inline void UpdateBack(HashT back, HashT front) {
    child_hash_.UpdateBack(back, front);
    hash_.UpdateBack(back, front);
  }

  inline HashTuple<size, HashT> GetHash() const {
    return HashTuple<size, HashT>(child_hash_.GetHash(), hash_.GetHash());
  }
  inline HashTuple<size, HashT> GetHash(HashT seed) const {
    return HashTuple<size, HashT>(child_hash_.GetHash(seed), hash_.GetHash(seed));
  }

  inline size_t GetHashInt() const {
    return child_hash_.GetHashInt() ^ hash_.GetHash();
  }


  inline bool operator ==(const MultiPolynomialHash<size, HashT> &other) const {
    return child_hash_ == other.child_hash_ && hash_ == other.hash_;
  }
  inline bool operator !=(const MultiPolynomialHash<size, HashT> &other) const {
    return !(operator==(other));
  }
};

template <class HashT>
class MultiPolynomialHash<1, HashT> {
  PolynomialHash<utils::PrimeHolder<0>::val, HashT> hash_;

 public:
  typedef HashTuple<1, HashT> DataType;

  MultiPolynomialHash(unsigned polydeg) : hash_(polydeg) {
  }

  MultiPolynomialHash(unsigned polydeg, const MultiPolynomialHash<1, HashT> &other)
      : hash_(polydeg, other.hash_) {
  }
  inline void CopyFrom(unsigned polydeg, const MultiPolynomialHash<1, HashT> &other) {
    hash_.CopyFrom(polydeg, other.hash_);
  };

  inline void Update(HashT front) {
    hash_.Update(front);
  }
  inline void Update(HashT front, HashT back) {
    hash_.Update(front, back);
  }
  inline void UpdateBack(HashT back) {
    hash_.UpdateBack(back);
  }
  inline void UpdateBack(HashT back, HashT front) {
    hash_.UpdateBack(back, front);
  }

  inline HashTuple<1, HashT> GetHash() const {
    return HashTuple<1, HashT>(hash_.GetHash());
  }
  
  inline HashTuple<1, HashT> GetHash(HashT seed) const {
    return HashTuple<1, HashT>(hash_.GetHash(seed));
  }

  inline size_t GetHashInt() const {
    return hash_.GetHash();
  }


  inline bool operator ==(const MultiPolynomialHash<1, HashT> &other) const {
    return hash_ == other.hash_;
  }
  inline bool operator !=(const MultiPolynomialHash<1, HashT> &other) const {
    return !(operator==(other));
  }
};

/*
template <class HashT>
class MultiPolynomialHash<2, HashT> {
  PolynomialHash<239, HashT> h1;
  PolynomialHash<269, HashT> h2;

 public:
  typedef std::pair<HashT, HashT> DataType;

  MultiPolynomialHash(unsigned polydeg) : h1(polydeg), h2(polydeg) {
  }

  MultiPolynomialHash(unsigned polydeg, const MultiPolynomialHash<2, HashT> &other)
      : h1(polydeg, other.h1),
        h2(polydeg, other.h2) {
  }
  inline void CopyFrom(unsigned polydeg, const MultiPolynomialHash<2, HashT> &other) {
    h1.CopyFrom(polydeg, other.h1);
    h2.CopyFrom(polydeg, other.h2);
  };

  inline void Update(HashT front) {
    h1.Update(front);
    h2.Update(front);
  }
  inline void Update(HashT front, HashT back) {
    h1.Update(front, back);
    h2.Update(front, back);
  }
  inline void UpdateBack(HashT back) {
    h1.UpdateBack(back);
    h2.UpdateBack(back);
  }
  inline void UpdateBack(HashT back, HashT front) {
    h1.UpdateBack(back, front);
    h2.UpdateBack(back, front);
  }

  inline std::pair<HashT, HashT> GetHash() const {
    return std::pair<HashT, HashT>(h1.GetHash(), h2.GetHash());
  }

  inline size_t GetHashInt() const {
    return h1.GetHash() ^ h2.GetHash();
  }


  inline bool operator ==(const MultiPolynomialHash<2, HashT> &other) const {
    return GetHash() == other.GetHash();
  }
  inline bool operator !=(const MultiPolynomialHash<2, HashT> &other) const {
    return !(operator==(other));
  }
};

template <class HashT = uint64_t>
class DoublePolynomialHash : public MultiPolynomialHash<2, HashT> {
  typedef MultiPolynomialHash<2, HashT> base;

 public:
  DoublePolynomialHash(unsigned polydeg) : base(polydeg) {
  }
  DoublePolynomialHash(unsigned polydeg, const DoublePolynomialHash<HashT> &other)
      : base(polydeg, other) {
  }

};
*/

}

namespace std {

template<size_t size, class HashT>
std::ostream& operator<<(std::ostream& os, const cap::MultiPolynomialHash<size, HashT> &hash) {
  os << hash.GetHashInt();
  return os;
}

}
