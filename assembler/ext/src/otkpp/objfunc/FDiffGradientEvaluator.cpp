
#include "FDiffGradientEvaluator.h"
#include "OTK.h"

#include <cstring>

using namespace boost::numeric::ublas;

FDiffGradientEvaluator::FDiffGradientEvaluator(Type type, FunctionEvaluator *fe) : 
  fEval_(fe), type_(type), x_hi_(new double[fe->getN()]) { }

FDiffGradientEvaluator::FDiffGradientEvaluator(const FDiffGradientEvaluator &e)
{
  fEval_ = e.fEval_;
  type_ = e.type_;
  x_hi_ = new double[fEval_->getN()];
}

FDiffGradientEvaluator::~FDiffGradientEvaluator()
{
  delete[] x_hi_;
}

int FDiffGradientEvaluator::getN() const
{
  return fEval_->getN();
}

void FDiffGradientEvaluator::setFEvaluator(FunctionEvaluator *fEval)
{
  fEval_ = fEval;
}

bool FDiffGradientEvaluator::usesFiniteDifference() const
{
  return true;
}

double *FDiffGradientEvaluator::eval_(const double *x,
                                      double *g) const
{
  if(type_ == FDiffGradientEvaluator::CENTRAL_2)
    return eval_central_2(x, g);
  else
    return eval_central_4(x, g);
}

double *FDiffGradientEvaluator::eval_central_2(const double *x,
                                               double *g) const
{
  const int N = getN();
  
  double fxm, fxp;
  double xi;
  
  memcpy(x_hi_, x, N * sizeof(double));
  
  for(int i = 0; i < N; i++)
  {
    xi = x[i];
    x_hi_[i] = xi + OTK::ROOT3_EPS;
    fxp = (*fEval_)(x_hi_);
    x_hi_[i] = xi - OTK::ROOT3_EPS;
    fxm = (*fEval_)(x_hi_);
    x_hi_[i] = xi;
    
    g[i] = (fxp - fxm) / (2.0 * OTK::ROOT3_EPS);
  }
  
  return g;
}

double *FDiffGradientEvaluator::eval_central_4(const double *x,
                                               double *g) const
{
  const int N = getN();
  
  double fxm, fxp, fxm2, fxp2;
  double xi;
  double r3;
  
  memcpy(x_hi_, x, N * sizeof(double));
  
  for(int i = 0; i < N; i++)
  {
    xi = x[i];
    x_hi_[i] = xi + 2.0*OTK::ROOT3_EPS;
    fxp = (*fEval_)(x_hi_);
    x_hi_[i] = xi - 2.0*OTK::ROOT3_EPS;
    fxm = (*fEval_)(x_hi_);
    x_hi_[i] = xi + OTK::ROOT3_EPS;
    fxp2 = (*fEval_)(x_hi_);
    x_hi_[i] = xi - OTK::ROOT3_EPS;
    fxm2 = (*fEval_)(x_hi_);
    x_hi_[i] = xi;
    
    r3 = 0.5 * (fxp - fxm);
    g[i] = ((4.0 / 3.0) * (fxp2 - fxm2) - (1.0 / 6.0) * r3) / OTK::ROOT3_EPS;
  }
  
  return g;
}
