
#ifndef PARTAN_H

#include <otkpp/solvers/native/NativeGradientSolver.h>

class LineMinimizer;

/// Implements the PARTAN (parallel tangents) algorithm.
class PARTAN : public NativeGradientSolver
{
 public:
  struct State : public Cloneable< State, NativeGradientSolver::State > { };
  
  PARTAN();
  ~PARTAN();
  
  std::string getName() const;
  const State &getState() const { return state_; }
  bool usesGradient() const;
  bool usesHessian() const;
 private:
  boost::numeric::ublas::vector< double > d_;
  LineMinimizer *lineMinimizer_;
  int j_;
  int iterState_;
  boost::numeric::ublas::vector< double > xMinus_;
  double fMinus_;
  boost::numeric::ublas::vector< double > gMinus_;
  State state_;
  boost::numeric::ublas::vector< double > y_;
  
  void doSetup_(const Function &objFunc,
                const boost::numeric::ublas::vector< double > &x0,
                const Solver::Setup &solverSetup,
                const Constraints &C);
  IterationStatus iterate_();
};

#define PARTAN_H

#endif
