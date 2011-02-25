#include "hash.hpp"
#include "parameters.hpp"

// seq[0]*31^(n-1) + seq[1]*31^(n-2) + ... + seq[n-1]*31^0
template< typename T >
unsigned int Hash<T>::operator() (const T &seq) const {
    unsigned int h = HASH_SEED;
	for (int i = 0; i < seq.len(); i++) {
		h = ((h << 5) - h) + seq[i];
	}
	return h;
}

template< typename T >
unsigned int HashSym<T>::operator() (const T &seq) const {
	Hash<T> h;
	return h(seq);// TODO ^ h(!seq);
}

template< typename T >
struct HashSymWeighted {
  unsigned int operator() (const T &seq) {
    // Will use frequency of seq
    return 0;
  }
};
