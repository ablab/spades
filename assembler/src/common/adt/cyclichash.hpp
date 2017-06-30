#ifndef CYCLICHASH_HPP
#define CYCLICHASH_HPP

#include <cstdint>
#include <cstdlib>

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
        : hashvalue(0),
          n_(n), r_(n % precision),
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

    // Add inchar as an input and remove outchar, the hashvalue is updated. This
    // function can be used to update the hash value from the hash value of
    // [outchar]ABC to the hash value of ABC[inchar]
    digest update(chartype outchar, chartype inchar) const {
        hashvalue = rol1(hashvalue) ^ roln(hasher_(outchar)) ^ hasher_(inchar);
        return hashvalue;
    }

    // This is the reverse of the update function.  This function can be used to
    // update the hash value from the hash value of ABC[inchar] to the hash
    // value of [outchar]ABC
    digest reverse_update(chartype outchar, chartype inchar) const {
        hashvalue ^= roln(hasher_(outchar)) ^ hasher_(inchar);
        hashvalue = ror1(hashvalue);
        return hashvalue;
    }

    // Add inchar as an input, this is used typically only at the start the hash
    // value is updated to that of a longer string (one where inchar was
    // appended)
    digest eat(chartype inchar) const {
        hashvalue = rol1(hashvalue) ^ hasher_(inchar);
        return hashvalue;
    }

    // Add inchar as an input and remove outchar, the hashvalue is updated. This
    // function can be used to update the hash value from the hash value of
    // [outchar]ABC to the hash value of ABC[inchar]
    digest hash_update(digest hashvalue, chartype outchar, chartype inchar) const {
        return rol1(hashvalue) ^ roln(hasher_(outchar)) ^ hasher_(inchar);
    }

    // For an n-gram X it returns hash value of (n + 1)-gram XY without changing
    // the object X. For example, if X = "ABC", then X.hash_extend("D") returns
    // value of "ABCD" without changing the state of X
    digest hash_extend(chartype Y) const {
        return rol1(hashvalue) ^ hasher_(Y);
    }

    //  Same as hash_extend, but with prepending the n-gram with character Y. If
    //  X = "ABC", then X.hash_prepend("D") returns value of "DABC" without
    //  changing the state of X
    digest hash_prepend(chartype Y) const {
        return roln(hasher_(Y)) ^ hashvalue;
    }

    digest hashvalue;
  private:
    const unsigned n_;
    const unsigned r_;

    hasher hasher_;
    const digest mask1_;
    const digest maskn_;
};

#endif 
