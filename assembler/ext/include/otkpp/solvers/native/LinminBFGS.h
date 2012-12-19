
#ifndef LINMINBFGS_H

#include <otkpp/linalg/InvLBFGSUpdater.h>
#include <otkpp/linalg/QuasiNewtonUpdater.h>
#include <otkpp/solvers/native/NativeGradientSolver.h>
#include <otkpp/solvers/native/linmin/LineMinimizer.h>

class LineMinimizer;

/// Implements the BFGS algorithm with line searches.
class LinminBFGS : public NativeGradientSolver
{
 public:
  /// Defines the parameters of a LinminBFGS solver.
  struct Setup : public Cloneable< Setup, Solver::Setup >
  {
    boost::shared_ptr< LineMinimizer::Setup > lmSetup; ///< line minimizer parameters
    boost::numeric::ublas::matrix< double > H0;                               ///< initial Hessian or inverse Hessian approximation
  
    Setup(const LineMinimizer::Setup &lmSetup = LineMinimizer::DefaultSetup(),
          const boost::numeric::ublas::matrix< double > &H0 = boost::numeric::ublas::zero_matrix< double >(0, 0));
    
    bool isCompatibleWith(const Solver &s) const;
  };

  struct State : public Cloneable< State, NativeGradientSolver::State >
  {
    boost::numeric::ublas::matrix< double > H;
  };
  
  /// Line search algorithm type.
  enum LinMinType
  {
    FLETCHER,     ///< Fletcher's line search.
    MORE_THUENTE  ///< The More and Thuente line search
  };
  
  /// Constructs a new (L-)BFGS solver.
  /**
   * Constructs a new BFGS solver with the specified line 
   * search algorithm and iteration history length.
   * @param lmType line search algorithm
   * @param iterHistLen the length of iteration history
   *                    (0 means full inverse Hessian approximation)
   */
  LinminBFGS(LinMinType lmType = MORE_THUENTE, int iterHistLen = 0);
  
  ~LinminBFGS();
  
  std::string getName() const;
  const State &getState() const { return state_; }
  bool usesGradient() const;
  bool usesHessian() const;
 private:
  boost::numeric::ublas::vector< double > d_;
  InvLBFGSUpdater dirUpdater_;
  int iterHistLen_;
  LineMinimizer *lineMinimizer_;
  LinMinType lmType_;
  QuasiNewtonUpdater *matrixUpdater_;
  boost::numeric::ublas::vector< double > p_;
  boost::numeric::ublas::vector< double > q_;
  State state_;
  
  void doSetup_(const Function &objFunc,
                const boost::numeric::ublas::vector< double > &x0,
                const Solver::Setup &solverSetup,
                const Constraints &C);
  IterationStatus iterate_();
};

#define LINMINBFGS_H

#endif
