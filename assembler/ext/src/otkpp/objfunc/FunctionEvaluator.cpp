
#include "FunctionEvaluator.h"

using namespace boost::numeric::ublas;

double FunctionEvaluator::operator()(const double *x) const
{
  if(evalCounting_)
    evalCounter_++;
  
  return eval_(x);
}

double FunctionEvaluator::operator()(const boost::numeric::ublas::vector< double > &x) const
{
  if(evalCounting_)
    evalCounter_++;
  
  return eval_(&x[0]);
}

bool FunctionEvaluator::hasSymbolicExpression() const
{
  return false;
}
