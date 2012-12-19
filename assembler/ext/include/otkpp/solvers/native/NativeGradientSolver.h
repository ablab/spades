
#ifndef NATIVEGRADIENTSOLVER_H

#include <otkpp/solvers/GradientSolverBase.h>

class LineMinimizer;

/// Defines a gradient-based solver.
/**
 * This class defines additional methods and attributes 
 * used by all native gradient-based solvers.
 */
class NativeGradientSolver : public GradientSolverBase
{
 public:
  struct State : public NativeSolver::State
  {
    double alpha;
    boost::numeric::ublas::vector< double > g;
  };
  
  /// Returns the current gradient boost::numeric::ublas::vector \f$\nabla f(\mathbf{x}_{k})\f$.
  virtual const boost::numeric::ublas::vector< double > getGradient() const;
  
  virtual const State &getState() const = 0;
  bool hasBuiltInStoppingCriterion() const;
 protected:
  double fPlus_;
  boost::numeric::ublas::vector< double > gPlus_;
  boost::numeric::ublas::vector< double > xPlus_;
  
  NativeGradientSolver(bool useFDiffGradient = false);
  
  virtual void doSetup_(const Function &objFunc,
                        const boost::numeric::ublas::vector< double > &x0,
                        const Solver::Setup &solverSetup,
                        const Constraints &C);
  State &getState_();
};

#define NATIVEGRADIENTSOLVER_H

#endif

