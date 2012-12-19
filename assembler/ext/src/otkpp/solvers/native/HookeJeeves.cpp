
#include "HookeJeeves.h"

#include <cmath>

using namespace boost::numeric::ublas;

std::string HookeJeeves::getName() const
{
  return "HookeJeeves";
}

bool HookeJeeves::hasBuiltInStoppingCriterion() const
{
  return true;
}

bool HookeJeeves::usesGradient() const
{
  return false;
}

bool HookeJeeves::usesHessian() const
{
  return false;
}

NativeSolver::IterationStatus HookeJeeves::iterate_()
{
  vector< double > d = zero_vector< double >(setup_->n);
  double ft, fy;
  vector< double > yt;
  
  for(int j = 0; j < setup_->n; j++)
  {
    d[j] = 1.0;
    
    fy = setup_->f(y_);
    yt = y_ + delta_*d;
    ft = setup_->f(yt);
    
    if(!std::isnan(ft) && !std::isinf(ft) && ft < fy)
      y_ = yt;
    else
    {
      yt = y_ - delta_*d;
      ft = setup_->f(yt);
      if(ft < fy)
        y_ = yt;
    }
    
    d[j] = 0.0;
  }
  
  if(!std::isnan(fy) && !std::isinf(fy) && fy < state_.fx)
  {
    xPlus_ = y_;
    y_ = xPlus_ + alpha_ * (xPlus_ - state_.x);
    state_.x = xPlus_;
    state_.fx = setup_->f(state_.x);
  }
  else
  {
    if(delta_ < eps_)
      return NativeSolver::ITERATION_SUCCESS;
    else
    {
      delta_ *= rho_;
      y_ = state_.x;
    }
  }
  
  return NativeSolver::ITERATION_CONTINUE;
}

void HookeJeeves::doSetup_(const Function &objFunc,
                           const vector< double > &x0,
                           const Solver::Setup &solverSetup,
                           const Constraints &C)
{
  NativeSolver::doSetup_(objFunc, x0, solverSetup, C);
  state_.x = x0;
  state_.fx = objFunc(x0);

  delta_ = 1.0;
  eps_ = 1e-10;
  alpha_ = 1.0;
  rho_ = 0.75;
  y_ = x0;
}
