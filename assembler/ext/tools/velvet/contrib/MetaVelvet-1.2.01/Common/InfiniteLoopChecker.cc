#include "InfiniteLoopChecker.hh"

bool InfiniteLoopChecker::check(){
  if( ++numLoops >= maxNumLoops ){
    cerr << "[InfiniteLoopChecker] " << "Warning: Infinite loop was found: (num. loops = " << numLoops << ")" << endl;
    return false;
  }
  return true;
}
