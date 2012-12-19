
#include "MoreThuente.h"
#include "PARTAN.h"

using namespace boost::numeric::ublas;

PARTAN::PARTAN()
{
  lineMinimizer_ = new MoreThuente();
}

PARTAN::~PARTAN()
{
  delete lineMinimizer_; 
}

std::string PARTAN::getName() const
{
  return "PARTAN";
}

bool PARTAN::usesGradient() const
{
  return true;
}

bool PARTAN::usesHessian() const
{
  return false;
}

NativeSolver::IterationStatus PARTAN::iterate_()
{
  double alpha;
  
  if(iterState_ == 1 || iterState_ == 2)
  {
    d_ = -state_.g;
    lineMinimizer_->minimize(state_.x, d_, 1.0, state_.fx, state_.g,
                             alpha, xPlus_, fPlus_, gPlus_);
    
    if(iterState_ == 1)
    {
      xMinus_ = state_.x;
      fMinus_ = state_.fx;
      gMinus_ = state_.g;
    }
  }
  else
  {
    d_ = state_.x - xMinus_;
    lineMinimizer_->minimize(xMinus_, d_, 1.0, fMinus_, gMinus_,
                             alpha, xPlus_, fPlus_, gPlus_);
    j_++;
    if(j_ == setup_->n)
    {
      j_ = 1;
      iterState_ = 0;
    }
    else
    {
      xMinus_ = xPlus_;
      fMinus_ = fPlus_;
      gMinus_ = gPlus_;
    }
  }
  
  iterState_++;
  if(iterState_ == 4)
    iterState_ = 1;

  state_.x = xPlus_;
  state_.fx = fPlus_;
  state_.g = gPlus_;
  
  return NativeSolver::ITERATION_CONTINUE;
}

void PARTAN::doSetup_(const Function &objFunc,
                      const vector< double > &x0,
                      const Solver::Setup &solverSetup,
                      const Constraints &C)
{
  NativeGradientSolver::doSetup_(objFunc, x0, solverSetup, C);
  lineMinimizer_->setup(setup_->f, MoreThuente::Setup());
  
  j_ = 1;
  d_ = -state_.g;
  iterState_ = 1;
}
