#include <iostream>
#include <vector>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include "trie.hpp"

void TestTrie() {
  ASSERT_EQUAL(1, 1);
}

cute::suite TrieSuite(){
	cute::suite s;
	s.push_back(CUTE(TestTrie));
	return s;
}
