
#include "LineMinimizer.h"

#include <cfloat>
#include <typeinfo>

using namespace boost::numeric::ublas;

double LineMinimizer::fletcherInitStep(double dphi0, double deltaF)
{
  if(dphi0 < -10.0 * DBL_EPSILON && deltaF > 10.0 * DBL_EPSILON)
    return -2.0 * deltaF / dphi0;
  else
    return 1.0;
}

void LineMinimizer::setup(const Function &f, const LineMinimizer::Setup &s)
{
  if(typeid(s) != typeid(const LineMinimizer::DefaultSetup &) && 
     !s.isCompatibleWith(*this))
    throw std::invalid_argument("the line minimizer and its setup are incompatible");
  
  f_ = &f;
  doSetup_(s);
}
