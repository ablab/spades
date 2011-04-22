#include "hash.hpp"

namespace hashing {

hash_t initH() {
	hash_t H = 1;
	for (int i = 0; i < K; i++) {
		H = HASH_X(H);
	}
	return H;
}

hash_t H = initH();

}
