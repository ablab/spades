
#include "cholesky.hpp"
#include "MNewton.h"
#include "MoreThuente.h"

using namespace boost::numeric::ublas;

MNewton::MNewton()
{
  lineMinimizer_ = new MoreThuente();
}

MNewton::~MNewton()
{
  delete lineMinimizer_;
}

std::string MNewton::getName() const
{
  return "MNewton";
}

bool MNewton::usesGradient() const
{
  return true;
}

bool MNewton::usesHessian() const
{
  return true;
}

NativeSolver::IterationStatus MNewton::iterate_()
{
  double alpha;
  
  setup_->f.H(state_.x, H_);
  d_ = -state_.g;
  int r = cholesky_decompose(H_, H_chol_);
  while(r != 0)
  {
    H_ = H_ + 1e-3 * identity_matrix< double >(setup_->n);
    r = cholesky_decompose(H_, H_chol_);
  }
  inplace_solve(H_chol_, d_, lower_tag());
  inplace_solve(trans(H_chol_), d_, upper_tag());
  
  lineMinimizer_->minimize(state_.x, d_, 1.0, state_.fx, state_.g,
                           alpha, xPlus_, fPlus_, gPlus_);
  
  state_.x = xPlus_;
  state_.fx = fPlus_;
  state_.g = gPlus_;
  
  return NativeSolver::ITERATION_CONTINUE;
}

void MNewton::doSetup_(const Function &objFunc,
                       const vector< double > &x0,
                       const Solver::Setup &solverSetup,
                       const Constraints &C)
{
  const int N = objFunc.getN();
  NativeGradientSolver::doSetup_(objFunc, x0, solverSetup, C);
  lineMinimizer_->setup(setup_->f, MoreThuente::Setup());
  H_.resize(N, N);
  H_chol_.resize(N, N);
}
