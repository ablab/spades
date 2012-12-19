
#include "GradientSolverBase.h"

using namespace boost::numeric::ublas;

GradientSolverBase::GradientSolverBase(Function::DerivEvalType gEvalType) : 
  gEvalType_(gEvalType) { }

void GradientSolverBase::doSetup_(const Function &objFunc,
                                  const vector< double > &x0,
                                  const Solver::Setup &solverSetup,
                                  const Constraints &C)
{
#ifdef WITH_LIBMATHEVAL
  if(gEvalType_ != Function::DERIV_SYMBOLIC)
#endif
    setup_->f = objFunc.createCopy(gEvalType_);
  setup_->f.resetEvalCounters();
  setup_->f.enableEvalCounting();
  NativeSolver::doSetup_(setup_->f, x0, solverSetup, C);
}
