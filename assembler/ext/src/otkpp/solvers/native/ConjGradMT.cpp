
#include "ConjGradMT.h"

using namespace boost::numeric::ublas;

ConjGradMT::ConjGradMT(Type type)
{
  type_ = type;
}

std::string ConjGradMT::getName() const
{
  if(type_ == FLETCHER_REEVES)
    return "ConjGrad-FR/MT";
  else
    return "ConjGrad-PR/MT";
}

bool ConjGradMT::usesGradient() const
{
  return true;
}

bool ConjGradMT::usesHessian() const
{
  return false;
}

NativeSolver::IterationStatus ConjGradMT::iterate_()
{
  double alpha, alpha0;
  double gamma;
  double gd;
  double gg;
  
  gg = inner_prod(state_.g, state_.g);
  gd = inner_prod(state_.g, state_.d);
  if(gd >= 0.0)
  {
    state_.d = -state_.g;
    gd = -gg;
  }
  
  alpha0 = LineMinimizer::fletcherInitStep(gd, deltaF_);
  lineMinimizer_.minimize(state_.x, state_.d, alpha0, state_.fx, state_.g,
                          alpha, xPlus_, fPlus_, gPlus_);
  
  if(type_ == FLETCHER_REEVES)
  {
    i_++;
    if(i_ == setup_->n + 1)
    {
      gamma = 0.0;
      i_ = 0;
    }
    else
      gamma = inner_prod(gPlus_, gPlus_ ) / gg;
  }
  else
  {
    gamma = std::max(inner_prod(gPlus_, gPlus_ - state_.g) / gg, 0.0);
    lineMinimizer_.setGamma(gamma);
  }
  state_.d = -gPlus_ + gamma*state_.d;
  
  deltaF_ = state_.fx - fPlus_;
  
  state_.x = xPlus_;
  state_.fx = fPlus_;
  state_.g = gPlus_;
  
  return NativeSolver::ITERATION_CONTINUE;
}

void ConjGradMT::doSetup_(const Function &objFunc,
                          const vector< double > &x0,
                          const Solver::Setup &solverSetup,
                          const Constraints &C)
{
  NativeGradientSolver::doSetup_(objFunc, x0, solverSetup, C);
  lineMinimizer_.setup(setup_->f, MoreThuente::Setup());
  state_.d = -state_.g;
  deltaF_ = 0.0;
  i_ = 0;
}
