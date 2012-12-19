
#include "GradNormTest.h"
#include <otkpp/solvers/native/NativeSolver.h>

#include <exception>
#include <stdexcept>

GradNormTest::GradNormTest(double eps) : eps_(eps) { }

double GradNormTest::getTestValue(const NativeSolver &s) const
{
  return norm_2(s.getGradient());
}

bool GradNormTest::test(const NativeSolver &s) const
{
  if(!s.usesGradient())
    throw std::runtime_error("no gradient information available");
  
  return getTestValue(s) < eps_;
}
