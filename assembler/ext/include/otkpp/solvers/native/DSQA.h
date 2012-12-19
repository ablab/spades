
#ifndef DSQA_H

#include <otkpp/interpolation/QuadInterp.h>
#include <otkpp/solvers/native/NativeSolver.h>

/// An experimental direct search algorithm with quadratic interpolation (similar to Powell's UOBYQA).
class DSQA : public NativeSolver
{
 public:
  struct State : public Cloneable< State, NativeSolver::State >
  {
    double delta;
    QuadInterp model;
  };
  
  std::string getName() const;
  const State &getState() const { return state_; }
  bool hasBuiltInStoppingCriterion() const;
  bool usesGradient() const;
  bool usesHessian() const;
 private:
  int m_;
  boost::numeric::ublas::vector< double > p_;
  State state_;
  boost::numeric::ublas::vector< double > xPlus_;
  
  double computeReduction_(const boost::numeric::ublas::vector< double > &x,
                           const boost::numeric::ublas::vector< double > &xPlus,
                           double f,
                           double fPlus,
                           const boost::numeric::ublas::vector< double > &p) const;
  double computeTau_(const boost::numeric::ublas::vector< double > &d,
                     const boost::numeric::ublas::vector< double > &p);
  boost::numeric::ublas::vector< double > &computeTrsRegStep_(const boost::numeric::ublas::vector< double > &g,
                                       const boost::numeric::ublas::matrix< double > &H,
                                       double delta,
                                       boost::numeric::ublas::vector< double > &p);
  NativeSolver::IterationStatus iterate_();
  void doSetup_(const Function &objFunc,
                const boost::numeric::ublas::vector< double > &x0,
                const Solver::Setup &solverSetup,
                const Constraints &C);
};

#define DSQA_H

#endif
