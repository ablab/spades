
#include "BFGSUpdater.h"
#include "DoglegBFGS.h"
#include "DoglegSolver.h"
#include "SteihaugSolver.h"

#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric::ublas;

DoglegBFGS::DoglegBFGS()
{
  //lmatrixUpdater_ = new LBFGSUpdater(2, 5);
}

DoglegBFGS::~DoglegBFGS()
{
  
}

std::string DoglegBFGS::getName() const
{
  return "BFGS/DLG";
}

bool DoglegBFGS::usesGradient() const
{
  return true;
}

bool DoglegBFGS::usesHessian() const
{
  return false;
}

NativeSolver::IterationStatus DoglegBFGS::iterate_()
{
  bool nonzeroStep;
  
  /*if(nIter_ > 5)
  {*/
  /*r = lmatrixUpdater_->computeProduct(p_, q_, g_, Hg_);
    if(r == false)
      Hg_ = g_;*/
  /*}
  else
    Hg_ = g_;*/
  
  /*std::cout<<"Hg: "<<Hg_<<std::endl;
  std::cout<<"Hg(correct): "<<prod(H_,g_)<<std::endl;*/
  
  trSolver_.computeStep(state_.x, state_.fx, state_.g, H_, p_, nonzeroStep, xPlus_, fPlus_);
  
  if(!nonzeroStep)
  {
    H_ = identity_matrix< double >(setup_->n);
    return NativeSolver::ITERATION_SUCCESS;
  }
  
  setup_->f.g(xPlus_, gPlus_);
  q_ = gPlus_ - state_.g;
  
  matrixUpdater_.update(p_, q_, H_);
  //lmatrixUpdater_->updateVectors(p_, q_);
  
  state_.x = xPlus_;
  state_.fx = fPlus_;
  state_.g = gPlus_;
  
  return NativeSolver::ITERATION_CONTINUE;
}

void DoglegBFGS::doSetup_(const Function &objFunc,
                          const vector< double > &x0,
                          const Solver::Setup &solverSetup,
                          const Constraints &C)
{
  const int n = objFunc.getN();
  
  NativeGradientSolver::doSetup_(objFunc, x0, solverSetup, C);
  
  matrixUpdater_ = BFGSUpdater(BFGSUpdater::DIRECT);
  trSolver_.setup(setup_->f);
  
  H_ = identity_matrix< double >(n);
  state_.g.resize(n);
  gPlus_.resize(n);
}
