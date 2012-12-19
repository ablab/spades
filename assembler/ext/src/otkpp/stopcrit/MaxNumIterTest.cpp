
#include "MaxNumIterTest.h"
#include "NativeSolver.h"

MaxNumIterTest::MaxNumIterTest(int n) : n_(n) { }

double MaxNumIterTest::getTestValue(const NativeSolver &s) const
{
  return s.getNumIter();
}

bool MaxNumIterTest::test(const NativeSolver &s) const
{
  return (s.getNumIter() >= n_);
}
