
#ifndef DOGLEGBFGS_H

#include <otkpp/linalg/BFGSUpdater.h>
#include <otkpp/linalg/LBFGSUpdater.h>
#include <otkpp/solvers/native/trsreg/DoglegSolver.h>
#include <otkpp/solvers/native/NativeGradientSolver.h>

/// Implements the BFGS algorithm with the "dogleg" strategy.
class DoglegBFGS : public NativeGradientSolver
{
 public:
  struct State : public Cloneable< State, NativeGradientSolver::State > { };
  
  DoglegBFGS();
  ~DoglegBFGS();
  
  std::string getName() const;
  const State &getState() const { return state_; }
  bool usesGradient() const;
  bool usesHessian() const;
 private:
  boost::numeric::ublas::matrix< double > H_;
  BFGSUpdater matrixUpdater_;
  //LBFGSUpdater *lmatrixUpdater_;
  boost::numeric::ublas::vector< double > p_;
  boost::numeric::ublas::vector< double > q_;
  State state_;
  DoglegSolver trSolver_;
  
  boost::numeric::ublas::vector< double > Hg_;

  IterationStatus iterate_();
  void doSetup_(const Function &objFunc,
                const boost::numeric::ublas::vector< double > &x0,
                const Solver::Setup &solverSetup,
                const Constraints &C);
};

#define DOGLEGBFGS_H

#endif
