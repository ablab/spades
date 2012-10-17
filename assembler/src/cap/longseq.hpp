#pragma once

namespace cap {

template <size_t size_, typename HashT = u_int64_t>
class LongSeq {
  friend class LongSeq<size_ - 1, HashT>;

  // Maybe SequenceData?
  Sequence sequence_;

  size_t pos_;
  char last_char_;
  HashT hash_;

  static const char kNoChar = char(-1);
  static const HashT kHashP = 239;
  static const HashT kHashPDegSize = InitHashPDegSize();

  inline static HashT InitHashPDegSize() const {
    HashT result = 1;
    for (size_t i = 1; i < size_; ++i) {
      result *= kHashP;
    }
    return result;
  }

  inline void LongSeq::InitHash() {
    hash_ = 0;
    for (size_t i = 0; i < size_; ++i) {
      hash_ = hash_ * kHashP + sequence_[i];
    }
  }

  // for fast copy
  LongSeq(const Sequence &sequence, size_t pos, HashT hash)
      : sequence_(sequence),
        pos_(pos),
        hash_(hash),
        last_char_(kNoChar) {
  }

  // convenient method
  inline char LastChar() const {
    if (last_char_ == kNoChar) {
      return operator[](size_ - 1);
    }
    return last_char_;
  }

 public:
  LongSeq(const Sequence &sequence, size_t pos = 0)
      : sequence_(sequence),
        pos_(pos),
        last_char_(kNoChar) {
    InitHash();
  }

  size_t size() const {
    return size_;
  }

  char operator [](size_t pos) const {
    return data_[pos_ + pos];
  }

  LongSeq<size_, HashT> operator !() const {
    return LongSeq<size_, HashT>(!sequence_, sequence_.size() - size_ - pos_);
  }

  void operator <<=(char c) {
    if (is_nucl(c)) {
      c = dignucl(c);
    }

    if (c == sequence_[pos_ + size_]) {
      // do nothing
    } else {
      // actually this can be only once during the life
      // of Kmer.
      VERIFY(last_char_ == kNoChar);
      last_char_ = c;
    }

    // Updating hash
    hash_ -= kHashPDegSize * sequence_[pos_];
    hash_ *= kHashP;
    hash_ += c;

    pos_++;
  }

  LongSeq<size_, HashT> operator <<(char c) const {
    LongSeq<size_, HashT> other_seq = *this;
    other_seq <<= c;
    return other_seq;
  }

  LongSeq<size_ + 1, HashT> pushBack(char c) const {
    if (is_nucl(c)) {
      c = dignucl(c);
    }

    LongSeq<size_+ 1, HashT> result(sequence_, pos_, hash_);

    if (c == sequence_[pos_ + size_ + 1]) {
      // do nothing
    } else {
      VERIFY(last_char_ == kNoChar);
      result.last_char_ = c;
    }

    result.hash_ *= kHashP;
    result.hash_ += c;

    return result;
  }

  bool operator ==(const LongSeq<size_, HashT> &other) const {
    if (hash_ != other.hash_) {
      return false;
    }
    for (size_t i = 0; i + 1 < size_; ++i) {
      if (operator[](i) != other[i]) {
        return false;
      }
    }
    return LastChar() == other.LastChar();
  }

  bool operator !=(const LongSeq<size_, HashT> &other) const {
    return !(operator==(other));
  }


  std::string str() const {
    std::string res(size_, '-');
    for (size_t i = 0; i < size_; ++i) {
      res[i] = nucl(operator[](i));
    }
    return res;
  }

  HashT GetHash() const {
    return hash_;
  }

};

}
