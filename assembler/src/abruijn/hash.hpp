#ifndef HASH_HPP_
#define HASH_HPP_

#include "logging.hpp"
#include "parameters.hpp"

namespace hashing {

LOGGER("a.hashing");

// Type of hash values.
typedef unsigned int hash_t;

/**
 * Maximum value of type hash_t
 */
const hash_t kMax = -1;

const hash_t kXor = 1845724623;

/**
 * Multiplies by x, fast.
 */
inline hash_t mult(hash_t v) {
	return (v << 5) - v;
}

/**
 * Calculates x^n.
 */
hash_t power(size_t k);

/**
 * x^K
 */
const hash_t kXK = power(K);

template< typename T >
struct Hash {
public:
	hash_t operator() (const T &seq) const {
		hash_t h = 0;
		for (size_t i = 0; i < seq.size(); i++) {
			h = mult(h) + seq[i];
		}
		return h ^ kXor;
	}
};

template< typename T >
struct HashSym {
	Hash<T> h;
public:
	hash_t operator() (const T &s) const {
		return h(s) ^ h(!s) ^ kXor;
	};

	/**
	 * Counts polynomial hashes of all k-mers in s, and puts it
	 * into array/container ha
	 */
	template<typename S>
	void kmers(const T &s, S &ha) {
		TRACE("hashing k-mers of " << s);
		size_t sz = s.size();
		hash_t h = 0;
		for (size_t i = 0; i < K; i++) {
			h = mult(h) + s[i];
		}
		ha[0] = h;
		for (size_t i = 0; i + K < sz; i++) {
			ha[i + 1] = mult(ha[i]) + s[i + K] - s[i] * kXK;
		}
		TRACE("forward pass - ok");
		h = 0;
		for (size_t i = sz - 1; i + K >= sz; i--) {
			h = mult(h) + (s[i] ^ 3);
		}
		for (size_t i = sz - K;; i--) {
			ha[i] ^= h ^ kXor;
			if (i == 0) {
				break;
			}
			h = mult(h) + (s[i - 1] ^ 3) - (s[i + K - 1] ^ 3) * kXK;
		}
		TRACE("hashing k-mers of " << s << " - done");
	}
};

template< typename T >
struct EqSym {
public:
	unsigned int operator() (const T &s1, const T &s2) const {
		return (s1 == s2) || (s1 == !s2);
	};
};

}

#endif /* HASH_HPP_ */
