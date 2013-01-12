
#include "XDistToMinTest.h"

using namespace boost::numeric::ublas;

XDistToMinTest::XDistToMinTest(const vector< double > &xMin, double eps, bool relative) : 
    eps_(eps), relative_(relative), xMin_(xMin) { }

double XDistToMinTest::getTestValue(const NativeSolver &s) const
{
  vector< double > dx = s.getX() - xMin_;
  
  if(relative_ == false)
    return norm_2(dx);
  else
    return norm_2(dx) / norm_2(xMin_);
}

const vector< double > &XDistToMinTest::getXMin() const
{
  return xMin_;
}

bool XDistToMinTest::test(const NativeSolver &s) const
{
  return (getTestValue(s) < eps_);
}
