
#ifndef GRADIENTSOLVERBASE_H

#include <otkpp/objfunc/Function.h>
#include <otkpp/solvers/native/NativeSolver.h>

/// Defines the interface for gradient-based solvers.
class GradientSolverBase : public NativeSolver
{
 public:
  /*/// Returns the current gradient vector \f$\nabla f(\mathbf{x}_{k})\f$.
  virtual const vector< double > getGradient() const = 0;*/
 protected:
  Function::DerivEvalType gEvalType_;
  
  GradientSolverBase(Function::DerivEvalType gEvalType);
  virtual void doSetup_(const Function &objFunc,
                        const boost::numeric::ublas::vector< double > &x0,
                        const NativeSolver::Setup &solverSetup,
                        const Constraints &C);
};

#define GRADIENTSOLVERBASE_H

#endif
