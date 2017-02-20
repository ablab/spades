//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <cstdlib>
#include <cstdint>
#include "polynomial_hash.hpp"
#include "math/log.hpp"
#include "sequence/sequence.hpp"
#include "utils/parallel/openmp_wrapper.h"

namespace cap {

typedef unsigned char uchar;

template <class HashT = PolynomialHash<> >
class LongSeq {
  typedef LongSeq<HashT> ThisType;
 public:
  /**
   * @variable Number of bits in type T (e.g. 8 for char)
   * @example 8: 2^8 = 256 or 16
   */
  const static size_t TBits = 8;

  /**
   * @variable Number of nucleotides that can be stored in one type T (e.g. 4 for char)
   * TNucl MUST be a power of two
   * @example 4: 8/2 = 4 or 16/2 = 8
   */
  const static size_t TNucl = 1;

  /**
   * @variable Number of bits in TNucl (e.g. 2 for char). Useful for shifts instead of divisions.
   */
  const static size_t TNuclBits = log_<TNucl, 2>::value;

  /**
   * @variable Number of Ts which required to store all sequence.
   */
  const static size_t DataSize = (16 + TNucl - 1) >> TNuclBits;

  /**
   * @variable Number of meaningful bytes in whick seq is stored
   */
  const static size_t TotalBytes = sizeof(char) * DataSize;

  typedef char DataType;
  typedef typename HashT::DataType HashType;

  static size_t GetDataSize(size_t size) {
    return (size + TNucl - 1) >> TNuclBits;
  }

  static const char kNoNextNucl = char(69);

 private:
  // Maybe SequenceData?
  Sequence sequence_;

  unsigned size_;
  size_t pos_;
  /*
   * last_char_ contains information about extra first/last chars
   * storage: xxxyxxxy
   *         frnt--back
   * xxx codes character itself
   * y codes existance of extra char
   */
  uchar last_char_;
  HashT hash_;

  static const uchar kNoFrontChar = uchar(8);
  static const uchar kNoBackChar = uchar(128);
  static const uchar kNoChar = uchar(255);

  inline void InitHash() {
    for (size_t i = 0; i < size_; ++i) {
      hash_.Update(sequence_[pos_ + i]);
    }
  }

  // for fast copy
  LongSeq(unsigned size, const Sequence &sequence, size_t pos, const HashT &hash)
      : sequence_(sequence),
        size_(size),
        pos_(pos),
        last_char_(kNoChar),
        hash_(size, hash) {
  }

  // convenient methods
  inline bool HasExtraLastChar() const {
    return !(last_char_ & kNoBackChar);
  }
  inline bool HasExtraFrontChar() const {
    return !(last_char_ & kNoFrontChar);
  }
  inline char LastChar() const {
    if (!HasExtraLastChar()) {
      return operator[](size_ - 1);
    }
    return (last_char_ >> 4) & 7;
  }
  inline char FirstChar() const {
    if (!HasExtraFrontChar()) {
      return operator[](0);
    }
    return last_char_ & 7;
  }
  inline void SetFirstChar(uchar c) {
    last_char_ = uchar((last_char_ & 0xF0) | c);
  }
  inline void SetLastChar(uchar c) {
    last_char_ = uchar((last_char_ & 0x0F) | (c << 4));
  }

 public:
  LongSeq()
      : sequence_(),
        size_(0),
        pos_(0),
        last_char_(kNoChar),
        hash_(0) {
  }
  explicit LongSeq(unsigned size)
      : sequence_(std::string(size, 'A')), // optimize by setting bad first char
        size_(size),
        pos_(0),
        last_char_(kNoChar),
        hash_(size_) {
  }
  LongSeq(unsigned size, const Sequence &sequence, size_t pos = 0)
      : sequence_(sequence),
        size_(size),
        pos_(pos),
        last_char_(kNoChar),
        hash_(size_) {
    InitHash();
    VERIFY(pos_ + size_ <= sequence_.size());
  }

  // Weird constructor for constructing from `data' origined from another LongSeq
  LongSeq(unsigned /* size */, const LongSeq<HashT> &other)
      : sequence_(other.sequence_),
        size_(other.size_),
        pos_(other.pos_),
        last_char_(other.last_char_),
        hash_(size_, other.hash_) {
  }

  LongSeq(const LongSeq<HashT> &other)
      : sequence_(other.sequence_),
        size_(other.size_),
        pos_(other.pos_),
        last_char_(other.last_char_),
        hash_(size_, other.hash_) {
  }

  LongSeq<HashT> operator =(const LongSeq<HashT> &other) {
    sequence_ = other.sequence_;
    size_ = other.size_;
    pos_ = other.pos_;
    last_char_ = other.last_char_;
    hash_.CopyFrom(size_, other.hash_);

    return *this;
  }

  unsigned size() const {
    return size_;
  }

  // TODO check consistensy
  const LongSeq<HashT> &data() const {
    return *this;
  }

  size_t data_size() const {
    return DataSize;
  }

  char operator [](size_t pos) const {
    return sequence_[pos_ + pos];
  }

  LongSeq<HashT> operator !() const {
    return LongSeq<HashT>(size_, !sequence_, sequence_.size() - size_ - pos_);
  }

  void Shift() {
    if (pos_ + size_ < sequence_.size()) {
      hash_.Update(sequence_[pos_ + size_], sequence_[pos_]);
    }
    pos_++;
  }

  void operator <<=(char c) {
    if (is_nucl(c)) {
      c = dignucl(c);
    }

    if (pos_ + size_ < sequence_.size() && c == sequence_[pos_ + size_]) {
      // do nothing
    } else {
      // actually this can be only once during the life
      // of Kmer.
      VERIFY(!HasExtraLastChar());
      SetLastChar(c);
    }

    char front_char;
    if (HasExtraFrontChar()) {
      front_char = FirstChar();
      SetFirstChar(kNoFrontChar);
    } else {
      front_char = sequence_[pos_];
    }
    // Updating hash
    hash_.Update(c, front_char);

    pos_++;

  }

  LongSeq<HashT> operator <<(char c) const {
    LongSeq<HashT> other_seq = *this;
    other_seq <<= c;
    return other_seq;
  }

  LongSeq<HashT> pushBack(char c) const {
    if (is_nucl(c)) {
      c = dignucl(c);
    }

    LongSeq<HashT> result(size_ + 1, sequence_, pos_, hash_);

    if (pos_ + size_ < sequence_.size() && c == sequence_[pos_ + size_]) {
      // do nothing
    } else {
      VERIFY(!HasExtraLastChar());
      result.SetLastChar(c);
    }

    result.hash_.Update(c);

    return result;
  }

  void operator >>=(char c) {
    if (is_nucl(c)) {
      c = dignucl(c);
    }

    if (pos_ > 0 && c == sequence_[pos_ - 1]) {
      // do nothing
    } else {
      // actually this can be only once during the life
      // of Kmer.
      VERIFY(!HasExtraFrontChar());
      SetFirstChar(c);
    }

    pos_--;

    char back_char;
    if (HasExtraLastChar()) {
      back_char = LastChar();
      SetLastChar(kNoFrontChar);
    } else {
      back_char = sequence_[pos_ + size_];
    }
    // Updating hash
    hash_.UpdateBack(c, back_char);

  }

  LongSeq<HashT> operator >>(char c) const {
    LongSeq<HashT> other_seq = *this;
    other_seq >>= c;
    return other_seq;
  }
  LongSeq<HashT> pushFront(char c) const {
    if (is_nucl(c)) {
      c = dignucl(c);
    }

    LongSeq<HashT> result(size_ + 1, sequence_, pos_ - 1, hash_);

    if (pos_ > 0 && c == sequence_[pos_ - 1]) {
      // do nothing
    } else {
      VERIFY(!HasExtraFrontChar());
      result.SetFirstChar(c);
    }

    result.hash_.UpdateBack(c);

    return result;
  }

  bool IsValid() const {
    return pos_ + size_ <= sequence_.size();
  }

  bool operator ==(const LongSeq<HashT> &other) const {

    VERIFY(size_ == other.size_);
    if (hash_ != other.hash_) {
//      (was a kind of joke, yea)
//      VERIFY(str() != other.str());
      return false;
    }
    for (size_t i = 1; i + 1 < size_; ++i) {
      if (operator[](i) != other[i]) {
//        VERIFY(str() != other.str());
        return false;
      }
    }
//    VERIFY((str() == other.str()) == (
//          FirstChar() == other.FirstChar() &&
//             LastChar() == other.LastChar()
//          ));
    return FirstChar() == other.FirstChar() &&
             LastChar() == other.LastChar();
  }

  bool operator !=(const LongSeq<HashT> &other) const {
    return !(operator==(other));
  }


  std::string str() const {
    if (size_ > 1) {
      std::string res(size_, '-');
      res[0] = nucl(FirstChar());
      res[size_ - 1] = nucl(LastChar());
      for (size_t i = 1; i + 1 < size_; ++i) {
        res[i] = nucl(operator[](i));
      }
      return res;
    }

    if (size_ == 0) {
      return "";
    }
    //if (size_ == 1)
    if (HasExtraFrontChar()) {
      return std::string("") + nucl(FirstChar());
    } else if (HasExtraLastChar()) {
      return std::string("") + nucl(LastChar());
    } else {
      return std::string("") + nucl(operator[](0));
    }
  }
  std::string err() const {
    std::ostringstream oss;
    oss << "{ sequence_=" << sequence_.err() <<
      "[" << sequence_.str() << "]" <<
      ", size_=" << size_ <<
      ", pos_=" << pos_ <<
      ", last_char_=" << size_t(last_char_) <<
      ", hash_=" << hash_.GetHash() << " }";
    return oss.str();
  }

  typename HashT::DataType GetHash() const {
    return hash_.GetHash();
  }
  typename HashT::DataType GetHash(typename HashT::AtomType seed) const {
    return hash_.GetHash(seed);
  }

  char GetNextNucl() const {
    if (pos_ + size_ >= sequence_.size()) {
      return kNoNextNucl;
    }
    return sequence_[pos_ + size_];
  }

  /*
  void UpdateTransition(char symbol, LongSeq<HashT> *link) {
    char my_symbol = GetNextNucl();
    if (symbol == my_symbol) {
      if (transition_storage_ == NULL) {
        transition_storage_ = std::shared_ptr<KmerJumper<ThisType> >(new SingleKmerJumper<ThisType>);
      }
      if (!transition_storage_->HasTransition()) {
        transition_storage_->SetTransition(symbol, link);
      }
      return;
    }

    if (transition_storage_ == NULL) {
      transition_storage_ = std::shared_ptr<KmerJumper<ThisType> >(new MultiKmerJumper<ThisType>);
    } else if (transition_storage_->Arity() == 1) {
      LongSeq<HashT> *old_link = transition_storage_->GetTransitionLink(my_symbol);
      transition_storage_ = std::shared_ptr<KmerJumper<ThisType> >(new MultiKmerJumper<ThisType>);
      if (my_symbol != kNoNextNucl) {
        transition_storage_->SetTransition(my_symbol, old_link);
      }
    }
    transition_storage_->SetTransition(symbol, link);
  }
  */

  struct hash {
    inline size_t operator()(const LongSeq<HashT>& seq) const {
      return seq.hash_.GetHashInt();
    }

    size_t operator()(const DataType * /*data*/, size_t /*sz*/ = DataSize) {
      VERIFY(false);
      return 0;
    }
  };

  struct equal_to {
    inline bool operator()(const LongSeq<HashT>& l, const LongSeq<HashT>& r) const {
      return l == r;
    }
  };

  struct fast_equal_to {
    inline bool operator()(const LongSeq<HashT>& l, const LongSeq<HashT>& r) const {
      return l.hash_ == r.hash_;
    }
  };

  struct less2 {
    bool operator()(const LongSeq<HashT> &l, const LongSeq<HashT> &r) const {
      VERIFY(l.size() == r.size());
      size_t size = l.size();
      for (size_t i = 1; i + 1 < size; ++i) {
        if (l[i] != r[i]) {
          return (l[i] < r[i]);
        }
      }
      return l.FirstChar() == r.FirstChar() && l.LastChar() < r.LastChar();
    }
  };

  bool operator<(const LongSeq<HashT>& that) const{
      static less2 comp;
      return comp(*this, that);
  }

  /**
   * Denotes some (weird) order on k-mers. Works fast.
  struct less2_fast {
    bool operator()(const LongSeq<HashT> &l, const LongSeq<HashT> &r) const {
      return l.GetHash() < r.GetHash();
    }
  };
   */

  bool IsMinimal() const {
      return true;
  }
};

typedef LongSeq<MultiPolynomialHash<3, uint64_t> > LSeq;

}

namespace std {

template<typename HashT>
std::ostream& operator<<(std::ostream& os, const cap::LongSeq<HashT> &seq) {
  os << seq.str();
  return os;
}

}
