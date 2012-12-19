
#include "BFGSUpdater.h"
#include "Fletcher.h"
#include "LinminBFGS.h"
#include "MoreThuente.h"

#include <typeinfo>

using namespace boost::numeric::ublas;

LinminBFGS::Setup::Setup(const LineMinimizer::Setup &lmSetup,
                         const matrix< double > &H0)
{
  this->lmSetup = boost::shared_ptr< LineMinimizer::Setup >(lmSetup.clone());
  this->H0 = H0;
}

bool LinminBFGS::Setup::isCompatibleWith(const Solver &s) const
{
  return (typeid(s) == typeid(const LinminBFGS &));
}

LinminBFGS::LinminBFGS(LinMinType lmType, int iterHistLen) : 
    iterHistLen_(iterHistLen), lmType_(lmType)
{
  if(lmType == LinminBFGS::FLETCHER)
    lineMinimizer_ = new Fletcher();
  else
    lineMinimizer_ = new MoreThuente();
  matrixUpdater_ = new BFGSUpdater(BFGSUpdater::INVERSE);
}

LinminBFGS::~LinminBFGS()
{
  delete lineMinimizer_;
  delete matrixUpdater_;
}

std::string LinminBFGS::getName() const
{
  if(lmType_ == LinminBFGS::FLETCHER)
  {
    if(iterHistLen_ == 0)
      return "BFGS/F";
    else
      return "L-BFGS/F";
  }
  else if(lmType_ == LinminBFGS::MORE_THUENTE)
  {
    if(iterHistLen_ == 0)
      return "BFGS/MT";
    else
      return "L-BFGS/MT";
  }

  return "";
}

bool LinminBFGS::usesGradient() const
{
  return true;
}

bool LinminBFGS::usesHessian() const
{
  return false;
}

NativeSolver::IterationStatus LinminBFGS::iterate_()
{
  double D;
  
  if(iterHistLen_ > 0)
  {
    if(state_.nIter >= 1)
      dirUpdater_.update(state_.g, d_);
    else
      d_ = -state_.g;
  }
  else
    d_ = -prod(state_.H, state_.g);
  
  lineMinimizer_->minimize(state_.x, d_, 1.0, state_.fx, state_.g,
                           state_.alpha, xPlus_, fPlus_, gPlus_);
  
  p_ = xPlus_ - state_.x;
  q_ = gPlus_ - state_.g;
  
  D = inner_prod(p_, q_);
  if(iterHistLen_ > 0)
  {
    if(D <= 0.0)
      d_ = -state_.g;
    
    dirUpdater_.storePair(p_, q_);
    dirUpdater_.removeOldestPair();
  }
  else
  {
    if(D <= 0.0)
      state_.H = identity_matrix< double >(setup_->n);
    else
      matrixUpdater_->update(p_, q_, state_.H);
  }
  
  state_.x = xPlus_;
  state_.fx = fPlus_;
  state_.g = gPlus_;
  
  return NativeSolver::ITERATION_CONTINUE;
}

void LinminBFGS::doSetup_(const Function &objFunc,
                          const vector< double > &x0,
                          const Solver::Setup &solverSetup,
                          const Constraints &C)
{
  const int n = objFunc.getN();
  
  NativeGradientSolver::doSetup_(objFunc, x0, solverSetup, C);
  
  if(typeid(solverSetup) == typeid(const Solver::DefaultSetup &))
  {
    state_.H = identity_matrix< double >(n);
    
    if(lmType_ == LinminBFGS::FLETCHER)
      lineMinimizer_->setup(setup_->f, Fletcher::Setup());
    else if(lmType_ == LinminBFGS::MORE_THUENTE)
      lineMinimizer_->setup(setup_->f, MoreThuente::Setup());
  }
  else
  {
    const LinminBFGS::Setup &setup = 
      dynamic_cast< const LinminBFGS::Setup & >(solverSetup);
    
    if(setup.H0.size1() == n && setup.H0.size2() == n)
      state_.H = setup.H0;
    else if(setup.H0.size1() == 0 && setup.H0.size2() == 0)
      state_.H = identity_matrix< double >(n);
    else
      throw std::invalid_argument("dimension mismatch");
    
    if(lmType_ == LinminBFGS::FLETCHER)
      lineMinimizer_->setup(setup_->f, *setup.lmSetup);
    else if(lmType_ == LinminBFGS::MORE_THUENTE)
      lineMinimizer_->setup(setup_->f, *setup.lmSetup);
  }
  
  if(iterHistLen_ > 0)
    dirUpdater_ = InvLBFGSUpdater(iterHistLen_, objFunc.getN());
}
