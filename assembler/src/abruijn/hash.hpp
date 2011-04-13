#ifndef HASH_HPP_
#define HASH_HPP_

#include "parameters.hpp"

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
public:
	hash_t operator() (const T &seq) const {
		Hash<T> h;
		return h(seq) ^ h(!seq) ^ HASH_XOR;
	};
};

template< typename T >
struct EqSym {
public:
	unsigned int operator() (const T &s1, const T &s2) const {
		return (s1 == s2) || (s1 == !s2);
	};
};

#endif /* HASH_HPP_ */
