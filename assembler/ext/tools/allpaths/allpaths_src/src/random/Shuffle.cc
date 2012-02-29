// Copyright (c) 2003 Broad Institute/Massachusetts Institute of Technology

#include "math/Functions.h"
#include "system/System.h"
#include "Vec.h"
#include "random/Random.h"
#include "random/Shuffle.h"
#include "random/MersenneTwister.h"

#include <algorithm>
#include <functional>
#include <stdlib.h>

/*
 * Shuffle 32bit
 */

class ShuffleRandom: public unary_function<int, int> {
 public:
  ShuffleRandom(longlong seed = 0) { srandomx(seed); }
  int operator()(int range = RAND_MAX) {
    if (range == RAND_MAX) return randomx();//save a modulo operation.
    return randomx() % range;
  }
};

void Shuffle( int N, vec<int> &shuffled, int seed )
{
  ForceAssertGt( N , 0 );

  shuffled.resize(N);
  for (int i=0; i<N; i++) shuffled[i] = i;

  ShuffleRandom generator(seed);

  random_shuffle(shuffled.begin(), shuffled.end(), generator);
}

/*
 * Shuffle 64bit
 */

// this uses the MersenneTwister
class ShuffleRandom64: public unary_function<uint64_t, uint64_t> {
 public:
  ShuffleRandom64(uint64_t seed = 0) { init_genrand64(seed); }
  uint64_t operator()(uint64_t range = -1) {
    if (range == static_cast<uint64_t>(-1)) return genrand64_int64();//save a modulo operation.
    return genrand64_int64() % range;
  }
};

void Shuffle64( uint64_t N, vec<uint64_t> &shuffled, uint64_t seed )
{
  ForceAssertGt( N , 0u );

  shuffled.resize(N);
  for (uint64_t i=0; i<N; i++) shuffled[i] = i;

  ShuffleRandom64 generator(seed);

  random_shuffle(shuffled.begin(), shuffled.end(), generator);
}
