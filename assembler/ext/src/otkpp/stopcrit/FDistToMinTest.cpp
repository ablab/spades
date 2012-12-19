
#include "FDistToMinTest.h"

FDistToMinTest::FDistToMinTest(double fMin, double eps, bool relative) : 
    eps_(eps), fMin_(fMin), relative_(relative) { }

double FDistToMinTest::getEps() const
{
  return eps_;
}

double FDistToMinTest::getFMin() const
{
  return fMin_;
}

double FDistToMinTest::getTestValue(const NativeSolver &s) const
{
  if(relative_ == false)
    return s.getFVal() - fMin_;
  else
  {
    double v = (s.getFVal() - fMin_) / fabs(fMin_);
    return v > 0.0 ? v : NAN;
  }
}

bool FDistToMinTest::isRelative() const
{
  return relative_;
}

bool FDistToMinTest::test(const NativeSolver &s) const
{
  return (getTestValue(s) < eps_);
}
