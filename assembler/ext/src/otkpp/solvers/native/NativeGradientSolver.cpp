
#include "NativeGradientSolver.h"

using namespace boost::numeric::ublas;

#ifdef WITH_LIBMATHEVAL
NativeGradientSolver::NativeGradientSolver(bool useFDiffGradient) : 
  GradientSolverBase(useFDiffGradient == false ? Function::DERIV_SYMBOLIC : Function::DERIV_FDIFF_CENTRAL_2) { }
#else
NativeGradientSolver::NativeGradientSolver(bool useFDiffGradient) : 
  GradientSolverBase(Function::DERIV_FDIFF_CENTRAL_2) { }
#endif

const vector< double > NativeGradientSolver::getGradient() const
{
  return getState().g;
}

NativeGradientSolver::State &NativeGradientSolver::getState_()
{
  return const_cast< NativeGradientSolver::State & >(getState());
}

bool NativeGradientSolver::hasBuiltInStoppingCriterion() const
{
  return false;
}

void NativeGradientSolver::doSetup_(const Function &objFunc,
                                    const vector< double > &x0,
                                    const Solver::Setup &solverSetup,
                                    const Constraints &C)
{
  GradientSolverBase::doSetup_(objFunc, x0, solverSetup, C);
  getState_().g.resize(objFunc.getN());
  objFunc.g(x0, getState_().g);
}

