
#include "HessianEvaluator.h"

using namespace boost::numeric::ublas;

matrix< double > &HessianEvaluator::operator()(const vector< double > &x,
                                               matrix< double > &H) const
{
  if(evalCounting_ && !usesFiniteDifference())
    evalCounter_++;
  
  return eval_(x, H);
}
