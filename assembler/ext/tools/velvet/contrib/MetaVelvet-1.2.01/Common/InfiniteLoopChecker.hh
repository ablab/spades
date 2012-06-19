#ifndef _INFINITE_LOOP_CHECKER_HH_
#define _INFINITE_LOOP_CHECKER_HH_
#include <stdlib.h>
#include <iostream>
#include <string>
#define INFINITE_LOOP_CHECKER_MAX_NUM_LOOPS 200

using namespace std;

class InfiniteLoopChecker {
  long numLoops;
  long maxNumLoops;
public:
  InfiniteLoopChecker( long num ) : numLoops(0), maxNumLoops(num) {}
  long getNumLoops() const { return numLoops; }
  bool check();
};

#endif // _INFINITE_LOOP_CHECKER_HH_
