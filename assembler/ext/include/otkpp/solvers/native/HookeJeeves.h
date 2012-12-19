
#ifndef HOOKEJEEVES_H

#include <otkpp/solvers/native/NativeSolver.h>

/// Implements the derivative-free Hooke-Jeeves algorithm.
class HookeJeeves : public NativeSolver
{
 public:
  struct State : public Cloneable< State, NativeSolver::State > { };
  
  std::string getName() const;
  const State &getState() const { return state_; }
  bool hasBuiltInStoppingCriterion() const;
  bool usesGradient() const;
  bool usesHessian() const;
 private:
  double alpha_;
  double delta_;
  double eps_;
  double rho_;
  State state_;
  boost::numeric::ublas::vector< double > xPlus_;
  boost::numeric::ublas::vector< double > y_;
  
  void doSetup_(const Function &objFunc,
                const boost::numeric::ublas::vector< double > &x0,
                const Solver::Setup &solverSetup,
                const Constraints &C);
  IterationStatus iterate_();
};

#define HOOKEJEEVES_H

#endif
