#ifndef HASH_HPP_
#define HASH_HPP_

#include "logging.hpp"
#include "parameters.hpp"

namespace hashing {

LOGGER("a.hashing");

extern hash_t H;

template< typename T >
struct Hash {
public:
	hash_t operator() (const T &seq) const {
		hash_t h = 0;
		for (size_t i = 0; i < seq.size(); i++) {
			h = HASH_X(h) + seq[i];
		}
		return h ^ HASH_XOR;
	}
};

template< typename T >
struct HashSym {
	Hash<T> h;
public:
	hash_t operator() (const T &s) const {
		return h(s) ^ h(!s) ^ HASH_XOR;
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
			h = HASH_X(h) + s[i];
		}
		ha[0] = h;
		for (size_t i = 0; i + K < sz; i++) {
			ha[i + 1] = HASH_X(ha[i]) + s[i + K] - s[i] * H;
		}
		TRACE("forward pass - ok");
		h = 0;
		for (size_t i = sz - 1; i + K >= sz; i--) {
			h = HASH_X(h) + (s[i] ^ 3);
		}
		for (size_t i = sz - K;; i--) {
			ha[i] ^= h ^ HASH_XOR;
			if (i == 0) {
				break;
			}
			h = HASH_X(h) + (s[i - 1] ^ 3) - (s[i + K - 1] ^ 3) * H;
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
