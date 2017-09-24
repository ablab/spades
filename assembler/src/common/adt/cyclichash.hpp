#pragma once

#include <cstdint>
#include <cstdlib>

//FIXME add namespace
template<typename chartype = uint8_t>
class DNASeqHash {
    /// The hash digest type.
    typedef uint64_t digest;

    static constexpr digest dA = 0x3c8bfbb395c60474;
    static constexpr digest dC = 0x3193c18562a02b4c;
    static constexpr digest dG = 0x20323ed082572324;
    static constexpr digest dT = 0x295549f54be24456;
    
    digest dtab_[256] = {
        [0] = dA,   [1] = dC, [2] = dG, [3] = dT,
        ['a'] = dA, ['c'] = dC, ['g'] = dG, ['t'] = dT,
        ['A'] = dA, ['C'] = dC, ['G'] = dG, ['T'] = dT
    };

  public:
    DNASeqHash(digest maxval, uint32_t seed = 42)
        : maxval_(maxval), seed_(seed) 
    {}

    digest operator()(chartype val) const {
        digest hval = dtab_[(unsigned)val];
        return hval ^ ((hval * seed_) >> 23);
    }

  private:
    digest maxval_;
    uint32_t seed_;
};

template<typename chartype = uint8_t>
class NDNASeqHash {
    /// The hash digest type.
    typedef uint64_t digest;

    static constexpr digest dA = 0x3c8bfbb395c60474;
    static constexpr digest dC = 0x3193c18562a02b4c;
    static constexpr digest dG = 0x20323ed082572324;
    static constexpr digest dT = 0x295549f54be24456;
    
  public:
    NDNASeqHash(digest maxval, uint32_t seed = 42)
        : maxval_(maxval), seed_(seed) 
    {}

    digest operator()(chartype val) const {
        digest hval = 0;
        switch (val) {
          case 0:
            hval = dA;
            break;
          case 1:
            hval = dC;
            break;
          case 2:
            hval = dG;
            break;
          case 3:
            hval = dT;
            break;
        }
        
        return hval ^ ((hval * seed_) >> 23);
    }

  private:
    digest maxval_;
    uint32_t seed_;
};

//FIXME bring to standard hash interface
template<unsigned precision = 64, typename chartype = uint8_t,
         typename hasher = DNASeqHash<chartype>>
class CyclicHash {
  public:
    /// The hash digest type.
    typedef uint64_t digest;
    typedef chartype char_t;

  private:
     constexpr digest mask(unsigned bits) const {
        return (digest(1) << (bits - 1)) ^
               ((digest(1) << (bits - 1)) - 1);
    }

    constexpr digest roln(digest x) const {
       return ((x & maskn_) << r_) | (x >> (precision-r_));
    }

    constexpr digest rol1(digest x) const {
        return ((x & mask1_) << 1) | (x >> (precision-1));
    }

    constexpr digest ror1(digest x) const {
        return (x >> 1) | ((x & 1) << (precision-1));
    }

  public:
    // @param n is the length of the sequences, e.g., 3 means that you want to hash
    // sequences of 3 characters
    CyclicHash(unsigned n)
        : n_(n), r_(n % precision),
          hasher_(mask(precision)),
          mask1_(mask(precision - 1)),
          maskn_(mask(precision - r_)) {
        static_assert(precision <= 8 * sizeof(digest), "Precision is too much");
    }

    // This is a convenience function, use eat, update and .hashvalue to use as
    // a rolling hash function
    template<class Seq>
    digest hash(const Seq &s) const {
        digest answer(0);
        for (size_t k = 0; k < s.size(); ++k)
            answer = rol1(answer) ^ hasher_(s[k]);

        return answer;
    }

    digest hashz(chartype outchar, unsigned n) const {
        digest answer = hasher_.hashvalues[static_cast<unsigned int>(outchar)];
        for (unsigned k = 0; k < n; ++k)
            answer = rol1(answer);
        return answer;
    }

//    // Add inchar as an input and remove outchar, the hashvalue is updated. This
//    // function can be used to update the hash value from the hash value of
//    // [outchar]ABC to the hash value of ABC[inchar]
//    digest update(chartype outchar, chartype inchar) const {
//        hashvalue = rol1(hashvalue) ^ roln(hasher_(outchar)) ^ hasher_(inchar);
//        return hashvalue;
//    }

//    // This is the reverse of the update function.  This function can be used to
//    // update the hash value from the hash value of ABC[inchar] to the hash
//    // value of [outchar]ABC
//    digest reverse_update(chartype outchar, chartype inchar) const {
//        hashvalue ^= roln(hasher_(outchar)) ^ hasher_(inchar);
//        hashvalue = ror1(hashvalue);
//        return hashvalue;
//    }

//    // Add inchar as an input, this is used typically only at the start the hash
//    // value is updated to that of a longer string (one where inchar was
//    // appended)
//    digest eat(chartype inchar) const {
//        hashvalue = rol1(hashvalue) ^ hasher_(inchar);
//        return hashvalue;
//    }

    // Add inchar as an input and remove outchar, the hashvalue is updated. This
    // function can be used to update the hash value from the hash value of
    // [outchar]ABC to the hash value of ABC[inchar]
    digest hash_update(digest hashvalue, chartype outchar, chartype inchar) const {
        return rol1(hashvalue) ^ roln(hasher_(outchar)) ^ hasher_(inchar);
    }

//    // For an n-gram X it returns hash value of (n + 1)-gram XY without changing
//    // the object X. For example, if X = "ABC", then X.hash_extend("D") returns
//    // value of "ABCD" without changing the state of X
//    digest hash_extend(chartype Y) const {
//        return rol1(hashvalue) ^ hasher_(Y);
//    }

//    //  Same as hash_extend, but with prepending the n-gram with character Y. If
//    //  X = "ABC", then X.hash_prepend("D") returns value of "DABC" without
//    //  changing the state of X
//    digest hash_prepend(chartype Y) const {
//        return roln(hasher_(Y)) ^ hashvalue;
//    }

//    digest hashvalue;
  private:
    const unsigned n_;
    const unsigned r_;

    hasher hasher_;
    const digest mask1_;
    const digest maskn_;
};

//Ask Anton about maxval and seed
template<typename chartype = uint8_t, typename digest = uint64_t>
class TrivialDNASeqHash {
    static constexpr digest dA = 0x3c8bfbb395c60474;
    static constexpr digest dC = 0x3193c18562a02b4c;
    static constexpr digest dG = 0x20323ed082572324;
    static constexpr digest dT = 0x295549f54be24456;
    
    //FIXME make static?
    digest dtab_[256] = {
        [0] = dA,   [1] = dC, [2] = dG, [3] = dT,
        ['a'] = dA, ['c'] = dC, ['g'] = dG, ['t'] = dT,
        ['A'] = dA, ['C'] = dC, ['G'] = dG, ['T'] = dT
    };

  public:
    TrivialDNASeqHash() {}

    digest operator()(chartype val) const {
        VERIFY(universal_is_nucl((char)val));
        return (digest) dtab_[(unsigned)val];
    }
};

//template<typename chartype = uint8_t,
//         typename digesttype = uint64_t,
//         typename hasher = TrivialDNASeqHash<chartype, digesttype>>
//class SymmetricCyclicHash {
//  public:
//    //todo rename
//    typedef digesttype digest;
//    typedef chartype char_t;
//  private:
//    static const unsigned precision = std::numeric_limits<digest>::digits;
//    hasher hasher_;
//    unsigned k_;
//    digest fwd_;
//    digest rvs_;
//    chartype first_;
//
//    digest rol(digest x, unsigned s = 1) {
//        return x << s | x >> (precision - s);
//    }
//
//    digest ror(digest x, unsigned s = 1) {
//        return x >> s | x << (precision - s);
//    }
//
//  public:
//    
//    /*
//     * Reimplementation of ntHash
//     * Adapted from https://bioinformatics.stackexchange.com/questions/19/are-there-any-rolling-hash-functions-that-can-hash-a-dna-sequence-and-its-revers
//     */
//    SymmetricCyclicHash(unsigned k): k_(k), fwd_(0), rvs_(0) {
//    }
//
//    template<class Seq>
//    digest operator()(const Seq &s) {
//        VERIFY(s.size() == k_);
//        fwd_ = 0;
//        rvs_ = 0;
//        for (size_t i = 0; i < s.size(); ++i)
//            fwd_ = rol(fwd_) ^ hasher_(s[i]);
//        for (size_t i = 0; i < s.size(); ++i)
//            rvs_ = rol(rvs_) ^ hasher_(universal_complement(s[k_-1-i]));
//        return operator()();   
//    }
//    
//    digest update(chartype outchar, chartype inchar) {
//        fwd_ = rol(fwd_) ^ rol(hasher_(outchar), k_) ^ hasher_(inchar);
//        rvs_ = ror(rvs_) ^ ror(hasher_(universal_complement(outchar))) ^ rol(hasher_(universal_complement(inchar)), k_ - 1);
//        return operator()();   
//    }
//
//    digest operator()() const {
//        return std::min(fwd_, rvs_);
//    }
//
//    //backwards compatibility
//    template<class Seq>
//    digest hash(const Seq &s) {
//        return operator()(s);
//    }
//
//    digest hash_update(digest hash, chartype outchar, chartype inchar) {
//        VERIFY(hash == operator()());
//        return update(outchar, inchar);
//    }
//
//};

template<typename chartype = uint8_t,
         typename digesttype = uint64_t,
         typename hasher = TrivialDNASeqHash<chartype, digesttype>>
class SymmetricCyclicHash {
  public:
    struct CyclicDigest {
        typedef digesttype digest;
        digesttype fwd;
        digesttype rvs;

        CyclicDigest() : fwd(0), rvs(0) {}

        digesttype value() const {
            return std::min(fwd, rvs);
        }

        explicit operator digesttype() {
            return value();
        }

        std::string str() const {
            std::stringstream os;
            os << "fwd " << fwd << "; rvs " << rvs << "; val " << value();
            return os.str();
        }
        
    };

    //todo rename
    typedef CyclicDigest digest;
    typedef chartype char_t;
  private:
    static const unsigned precision = std::numeric_limits<digesttype>::digits;
    hasher hasher_;
    unsigned k_;

    static digesttype rol(digesttype x, unsigned s = 1) {
        return x << s | x >> (precision - s);
    }

    static digesttype ror(digesttype x, unsigned s = 1) {
        return x >> s | x << (precision - s);
    }

  public:
    
    /*
     * Reimplementation of ntHash
     * Adapted from https://bioinformatics.stackexchange.com/questions/19/are-there-any-rolling-hash-functions-that-can-hash-a-dna-sequence-and-its-revers
     */
    SymmetricCyclicHash(unsigned k): k_(k) {
    }

    template<class Seq>
    digest operator()(const Seq &s) const {
        //std::cout << "Hashing " << s << std::endl;
        VERIFY(k_ <= s.size());
        digest answer;
        for (size_t i = 0; i < k_; ++i)
            answer.fwd = rol(answer.fwd) ^ hasher_(s[i]);
        for (size_t i = 0; i < k_; ++i) {
            //std::cout << "nucl " << nucl(s[k_-1-i]) << " ; compl " << nucl(universal_complement(s[k_-1-i])) << std::endl;
            //std::cout << "hasher val " << hasher_(universal_complement(s[k_-1-i])) << std::endl;
            answer.rvs = rol(answer.rvs) ^ hasher_(universal_complement(s[k_-1-i]));
        }
        return answer;   
    }

    //backwards compatibility
    template<class Seq>
    digest hash(const Seq &s) const {
        return operator()(s);
    }

    digest hash_update(digest hash, chartype outchar, chartype inchar) const {
        digest answer;
        answer.fwd = rol(hash.fwd) ^ rol(hasher_(outchar), k_) ^ hasher_(inchar);
        answer.rvs = ror(hash.rvs) ^ ror(hasher_(universal_complement(outchar))) ^ rol(hasher_(universal_complement(inchar)), k_ - 1);
        return answer;
    }

};

