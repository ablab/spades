
#include "ModelFunction.h"

using namespace boost::numeric::ublas;

ModelFunction::ModelFunction(const matrix< double > &A, vector< double > &b) : 
  A_(A), b_(b) { }

double ModelFunction::operator()(const vector< double > &x) const
{
  // NOTE: the constant is not included
  return inner_prod(x, prod(A_, x))/2.0 + inner_prod(b_, x);
}

const vector< double > &ModelFunction::g(const vector< double > &x, vector< double > &g)
{
  g = prod(A_, x) + b_;
  return g;
}
