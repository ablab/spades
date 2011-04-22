#include "hash.hpp"

namespace hashing {

hash_t power(size_t k) {
	hash_t h = 1;
	for (size_t i = 0; i < k; i++) {
		h = mult(h);
	}
	return h;
}

}
