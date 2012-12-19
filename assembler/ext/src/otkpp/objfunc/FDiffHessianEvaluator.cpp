
#include "FDiffHessianEvaluator.h"
#include "OTK.h"

using namespace boost::numeric::ublas;

FDiffHessianEvaluator::FDiffHessianEvaluator(FunctionEvaluator *fe) : fEval_(fe)
{
  f_hi_.resize(fe->getN());
}

matrix< double > &FDiffHessianEvaluator::eval_(const vector< double > &x,
                                               matrix< double > &H) const
{
  const int N = getN();
  
  double fx_ = (*fEval_)(x); // TODO: use a previously computed value
  
  int i, j;

  x_hi_   = x;
  x_hj_   = x;
  x_hihj_ = x;
  
  for(i = 0; i < N; i++)
  {
    x_hi_[i] = x[i] + OTK::SQRT_EPS;
    f_hi_[i] = (*fEval_)(x_hi_);
    x_hi_[i] = x[i];
  }
  
  for(i = 0; i < N; i++)
  {
    x_hi_[i] = x[i] + OTK::SQRT_EPS;
    
    for(j = 0; j <= i; j++)
    {
      x_hj_[j] = x[j] + OTK::SQRT_EPS;
      
      if(i != j)
      {
        x_hihj_[i] = x_hi_[i];
        x_hihj_[j] = x_hj_[j];
      }
      else
        x_hihj_[i] = x[i] + 2.0 * OTK::SQRT_EPS;

      H(i, j) = ((*fEval_)(x_hihj_) - f_hi_[i] - (f_hi_[j] - fx_)) / 
        (OTK::SQRT_EPS * OTK::SQRT_EPS);
      
      x_hj_[j]   = x[j];
      x_hihj_[i] = x[i];
      x_hihj_[j] = x[j];
    }
    
    x_hi_[i] = x[i];
  }
  
  for(i = 0; i < N; i++)
    for(j = i + 1; j < N; j++)
      H(i, j) = H(j, i);
  
  return H;
}
