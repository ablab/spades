#pragma once

#include <cstdint>
#include <cstdlib>

namespace rolling_hash {

/// The hash digest type.
typedef uint64_t digest;
typedef uint8_t chartype;

class NDNASeqHash {
    static constexpr digest dA = 0x3c8bfbb395c60474;
    static constexpr digest dC = 0x3193c18562a02b4c;
    static constexpr digest dG = 0x20323ed082572324;
    static constexpr digest dT = 0x295549f54be24456;

public:
    NDNASeqHash(uint32_t seed = 0)
            : seed_(seed) {}

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
            default:
                VERIFY_DEV(false);
                break;
        }

        return hval ^ ((hval * seed_) >> 23);
    }

private:
    uint32_t seed_;
};

class DNASeqHash {
    NDNASeqHash inner_hash_;

public:
    DNASeqHash(uint32_t seed = 0)
            : inner_hash_(seed) {}

    digest operator()(chartype val) const {
        return inner_hash_(dignucl(val));
    }
};


//FIXME bring to standard hash interface
//FIXME extract commented code to separate class
//FIXME does it make sense to have a different precision?
template<typename hasher = NDNASeqHash, unsigned precision = 64>
class CyclicHash {
    constexpr digest mask(unsigned bits) const {
        return (digest(1) << (bits - 1)) ^
               ((digest(1) << (bits - 1)) - 1);
    }

    constexpr digest roln(digest x) const {
        return ((x & maskn_) << r_) | (x >> (precision - r_));
    }

    constexpr digest rol1(digest x) const {
        return ((x & mask1_) << 1) | (x >> (precision - 1));
    }

    constexpr digest ror1(digest x) const {
        return (x >> 1) | ((x & 1) << (precision - 1));
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
    digest operator()(const Seq &s) const {
        VERIFY(n_ <= s.size());
        digest answer(0);
        for (size_t k = 0; k < n_; ++k)
            answer = rol1(answer) ^ hasher_(s[k]);

        return answer;
    }

    template<class Seq>
    digest hash(const Seq &s) const {
        return operator()(s);
    }

//    digest hashz(chartype outchar, unsigned n) const {
//        digest answer = hasher_.hashvalues[static_cast<unsigned int>(outchar)];
//        for (unsigned k = 0; k < n; ++k)
//            answer = rol1(answer);
//        return answer;
//    }

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

template<typename hasher = NDNASeqHash>
class SymmetricCyclicHash {
public:
    struct CyclicDigest {
        digest fwd;
        digest rvs;

        CyclicDigest() : fwd(0), rvs(0) {}

        digest value() const {
            return fwd + rvs;
        }

        explicit operator digest() {
            return value();
        }

        std::string str() const {
            std::stringstream os;
            os << "fwd " << fwd << "; rvs " << rvs << "; val " << value();
            return os.str();
        }
    };

private:
    static constexpr uint8_t precision = uint8_t(std::numeric_limits<digest>::digits);
    hasher hasher_;
    unsigned n_;

    static digest rol(digest x, unsigned s = 1) {
        return x << s | x >> (precision - s);
    }

    static digest ror(digest x, unsigned s = 1) {
        return x >> s | x << (precision - s);
    }

public:

    /*
     * Reimplementation of ntHash
     * Adapted from https://bioinformatics.stackexchange.com/questions/19/are-there-any-rolling-hash-functions-that-can-hash-a-dna-sequence-and-its-revers
     */
    SymmetricCyclicHash(unsigned n)
            : n_(n) {}

    template<class Seq>
    CyclicDigest operator()(const Seq &s) const {
        VERIFY(n_ <= s.size());
        CyclicDigest answer;
        for (size_t i = 0; i < n_; ++i)
            answer.fwd = rol(answer.fwd) ^ hasher_(s[i]);

        for (size_t i = 0; i < n_; ++i)
            answer.rvs = rol(answer.rvs) ^ hasher_(nucl_complement(s[n_ - 1 - i]));

        return answer;
    }

    template<class Seq>
    CyclicDigest hash(const Seq &s) const {
        return operator()(s);
    }

    CyclicDigest hash_update(CyclicDigest hash, chartype outchar, chartype inchar) const {
        CyclicDigest answer;
        answer.fwd = rol(hash.fwd) ^ rol(hasher_(outchar), n_) ^ hasher_(inchar);
        answer.rvs = ror(hash.rvs) ^ ror(hasher_(nucl_complement(outchar))) ^
                     rol(hasher_(nucl_complement(inchar)), n_ - 1);
        return answer;
    }

};
}
