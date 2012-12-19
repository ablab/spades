
#include "GradientEvaluator.h"

using namespace boost::numeric::ublas;

double *GradientEvaluator::operator()(const double *x,
                                      double *g) const
{
  if(evalCounting_ && !usesFiniteDifference())
    evalCounter_++;
  
  eval_(x, g);
  return g;
}

vector< double > &GradientEvaluator::operator()(const vector< double > &x,
                                                vector< double > &g) const
{
  if(evalCounting_ && !usesFiniteDifference())
    evalCounter_++;
  
  eval_(&x[0], &g[0]);
  return g;
}
