
#ifndef MNEWTON_H

#include <otkpp/solvers/native/NativeGradientSolver.h>

class Solver;

/// Implements the Newton method with line searches.
class MNewton : public NativeGradientSolver
{
 public:
   struct State : public Cloneable< State, NativeGradientSolver::State > { };
  
  MNewton();
  ~MNewton();
  
  std::string getName() const;
  const State &getState() const { return state_; }
  bool usesGradient() const;
  bool usesHessian() const;
 private:
  boost::numeric::ublas::vector< double > d_;
  boost::numeric::ublas::matrix< double > H_;
  boost::numeric::ublas::matrix< double > H_chol_;
  LineMinimizer *lineMinimizer_;
  State state_;
  
  void doSetup_(const Function &objFunc,
                const boost::numeric::ublas::vector< double > &x0,
                const Solver::Setup &solverSetup,
                const Constraints &C);
  IterationStatus iterate_();
};

#define MNEWTON_H

#endif
