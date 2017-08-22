//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "storing_traits.hpp"

namespace utils {

template<typename Key, class HashFunction>
class SimpleKeyWithHash {
public:
    typedef Key KeyType;
private:
    typedef typename HashFunction::IdxType IdxType;
    const HashFunction &hash_;
    Key key_;
    mutable IdxType idx_; //lazy computation
    mutable bool ready_;

    void CountIdx() const {
        ready_ = true;
        idx_ = hash_.seq_idx(key_);
    }

    void SetKey(const Key &key) {
        ready_ = false;
        key_ = key;
    }
public:

    SimpleKeyWithHash(Key key, const HashFunction &hash)
            : hash_(hash), key_(key), idx_(0), ready_(false) {}

    Key key() const {
        return key_;
    }

    IdxType idx() const {
        if (!ready_)
            CountIdx();

        return idx_;
    }

    SimpleKeyWithHash &operator=(const SimpleKeyWithHash &that) {
        VERIFY(&this->hash_ == &that.hash_);
        this->key_= that.key_;
        this->idx_ = that.idx_;
        this->ready_ = that.ready_;
        return *this;
    }

    bool operator==(const SimpleKeyWithHash &that) const {
        VERIFY(&this->hash_ == &that.hash_);
        if (this->ready_ && that.ready_)
            return this->idx_ == that.idx_ && this->is_minimal_ == that.is_minimal_;
        return this->key_ == that.key_;
    }

    bool operator!=(const SimpleKeyWithHash &that) const {
        VERIFY(&this->hash_ == &that.hash_);
        return this->key_ != that.key_;
    }

    SimpleKeyWithHash operator!() const {
        return SimpleKeyWithHash(!key_, hash_);
    }

    SimpleKeyWithHash operator<<(char nucl) const {
        return SimpleKeyWithHash(key_ << nucl, hash_);
    }

    SimpleKeyWithHash operator>>(char nucl) const {
        return SimpleKeyWithHash(key_ >> nucl, hash_);
    }

    void operator<<=(char nucl) {
        SetKey(key_ << nucl);
    }

    void operator>>=(char nucl) {
        SetKey(key_ >> nucl);
    }

    char operator[](size_t i) const {
        return key_[i];
    }

    bool is_minimal() const {
        return true;
    }
};

template<class stream, class Key, class Index>
stream &operator<<(stream &s, const SimpleKeyWithHash<Key, Index> &kwh) {
    return s << "SKWH[" << kwh.key() << ", " << kwh.idx() << "]";
}

//Would it make sense to also store inverted kmer for not minimal kwh?
template<typename Key, class HashFunction>
class InvertableKeyWithHash {
private:
    typedef typename HashFunction::IdxType IdxType;

    const HashFunction &hash_;
    Key key_;
    mutable IdxType idx_; //lazy computation
    mutable bool is_minimal_;
    mutable bool ready_;

    void CountIdx() const {
        ready_ = true;
        is_minimal_ = key_.IsMinimal();
        if (is_minimal_)
            idx_ = hash_.seq_idx(key_);
        else{
            idx_ = hash_.seq_idx(!key_);
        }
    }

    InvertableKeyWithHash(Key key, const HashFunction &hash, bool is_minimal,
                          size_t idx, bool ready)
            : hash_(hash), key_(key), idx_(idx),
              is_minimal_(is_minimal), ready_(ready) {
    }
  public:

    InvertableKeyWithHash(Key key, const HashFunction &hash)
            : hash_(hash), key_(key), idx_(0), is_minimal_(false), ready_(false) {}

    const Key &key() const {
        return key_;
    }

    IdxType idx() const {
        if (!ready_)
            CountIdx();

        return idx_;
    }

    bool is_minimal() const {
        if(!ready_) {
            return key_.IsMinimal();
        }
        return is_minimal_;
    }

    bool ready() const {
        return ready_;
    }

    InvertableKeyWithHash &operator=(const InvertableKeyWithHash &that) {
        VERIFY(&this->hash_ == &that.hash_);
        this->key_= that.key_;
        this->idx_ = that.idx_;
        this->ready_ = that.ready_;
        this->is_minimal_ = that.is_minimal_;
        return *this;
    }

    bool operator==(const InvertableKeyWithHash &that) const {
        VERIFY(&this->hash_ == &that.hash_);
        return this->key_ == that.key_;
    }

    bool operator!=(const InvertableKeyWithHash &that) const {
        VERIFY(&this->hash_ == &that.hash_);
        return this->key_ != that.key_;
    }

    InvertableKeyWithHash operator!() const {
        if (!ready_)
            return InvertableKeyWithHash(!key_, hash_);
        return InvertableKeyWithHash(!key_, hash_, !is_minimal_, idx_, ready_);
    }

    InvertableKeyWithHash operator<<(char nucl) const {
        return InvertableKeyWithHash(key_ << nucl, hash_);
    }

    InvertableKeyWithHash operator>>(char nucl) const {
        return InvertableKeyWithHash(key_ >> nucl, hash_);
    }

    void operator<<=(char nucl) {
        key_ <<= nucl;
        ready_ = false;
    }

    void operator>>=(char nucl) {
        key_ >>= nucl;
        ready_ = false;
    }

    char operator[](size_t i) const {
        return key_[i];
    }
};

template<class stream, class Key, class Index>
stream &operator<<(stream &s, const InvertableKeyWithHash<Key, Index> &kwh) {
    s << "IKWH[" << kwh.key();
    if(kwh.ready()) {
        return s << ", " << kwh.is_minimal() << ", " << kwh.idx() << "]";
    } else {
        return s << ", not ready]";
    }
}

template<class K, class Index, class StoringType>
struct StoringTraits;

template<class K, class Index>
struct StoringTraits<K, Index, SimpleStoring> {
    typedef SimpleKeyWithHash<K, Index> KeyWithHash;
};

template<class K, class Index>
struct StoringTraits<K, Index, InvertableStoring> {
    typedef InvertableKeyWithHash<K, Index> KeyWithHash;
};

}
