/* test for cuckoo.hpp */
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include "cuckoo.hpp"

void CuckooTest() {
	struct Hasher {
		size_t operator()(int value, int hash_num) {
			long k = 1;
			long l = 3 + hash_num;
			long t = 1 << 16;
			k = (value * l) % t;
			return (size_t) k;
		}
	};
	srand(239);
	/*cuckoo<int, int, Hasher, std::equal_to<int>, 4, 10, 100> Cuckoo;
	for (int i = 0; i < 10000; ++i) {
		Cuckoo.insert(std::make_pair(rand(), 17));
	}
	std::cerr << "Number of elements is " << Cuckoo.size() << std::endl
			<< "Actual size is " << Cuckoo.length() << std::endl;*/
}

cute::suite CuckooSuite() {
	cute::suite s;
	s.push_back(CUTE(TestIFastaStreamNoFile));
	s.push_back(CUTE(TestIFastaStreamSingleRead));
	s.push_back(CUTE(TestIFastaStreamFull));
	return s;
}
