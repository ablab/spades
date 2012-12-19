
#ifndef NATIVESOLVER_H

#include <otkpp/solvers/Solver.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/shared_ptr.hpp>
#include <list>
#include <ostream>

/// Defines the base class for the "native" solvers implemented in OTK++.
/**
 * The solvers derived from this class are the "native" solvers that are written 
 * in C++ by the author. These algorithms find a find a (local) minimum of a 
 * function \f$f:\mathbb{R}^n\rightarrow\mathbb{R}\f$ with simple constraints.
 */
class NativeSolver : public Solver
{
 public:
  /// Defines the status of the iteration.
  enum IterationStatus
  {
    ITERATION_CONTINUE,       ///< iteration is expected to progress towards a solution
    ITERATION_SUCCESS,        ///< iteration finished, success
    ITERATION_NO_PROGRESS,    ///< iteration is not making any progress
    ITERATION_OUT_OF_CONTROL  ///< iteration diverges or has reached infinity
  };
  
  /// Defines the state of the solver.
  struct State
  {
    double fx;            ///< the function value \f$f(\mathbf{x}_{k})\f$ at the current iterate
    unsigned int nIter;   ///< the number of iterations used so far
    boost::numeric::ublas::vector< double > x;   ///< the current iterate \f$\mathbf{x}_{k}\f$
    boost::numeric::ublas::matrix< double > X;   ///< the current iterates \f$\mathbf{x}_{k}^{i}\f$, \f$i=1,\dots,m\f$ (if more than one point is used at each step)
    
    virtual ~State() { }
    virtual State *clone() const = 0;
  };
  
  /// Defines the results produced by the solver.
  /**
   * In addition to the information produced by a generic solver, the solvers 
   * derived from this class also store the state of the solver after each 
   * iteration step.
   */
  struct Results : public Solver::Results
  {
    std::vector< boost::shared_ptr< State > > states;
  };
  
  NativeSolver() { }
  virtual ~NativeSolver() { }
  
  /// Returns the current iterate \f$\mathbf{x}_{k}\f$.
  virtual const boost::numeric::ublas::vector< double > getX() const;
  
  /// Returns an array of the current iterates.
  /**
   * If this solver generates multiple points per iteration, 
   * this method returns a matrix of them with each column 
   * containing one point. If this solver generates only 
   * one point per iteration step, this method returns 
   * a nx1 matrix identical to the vector returned by getX.
   */
  virtual const boost::numeric::ublas::matrix< double > getXArray() const;
  
  /// Returns the current function value \f$f(\mathbf{x}_{k})\f$.
  virtual double getFVal() const;
  
  /// Returns the current gradient \f$\nabla f(\mathbf{x}_{k})\f$.
  virtual const boost::numeric::ublas::vector< double > getGradient() const;
  
  /// Returns the current Hessian \f$H_{f}(\mathbf{x}_{k})\f$.
  virtual const boost::numeric::ublas::matrix< double > getHessian() const;
  
  /// Returns the number of iterations since the last setup of this solver.
  unsigned int getNumIter() const;
  
  /// Returns the number of function evaluations since the last setup of this solver.
  unsigned int getNumFuncEval() const;
  
  /// Returns the number of gradient evaluations since the last setup of this solver.
  unsigned int getNumGradEval() const;
  
  /// Returns the number of Hessian evaluations since the last setup of this solver.
  unsigned int getNumHessEval() const;
  
  /// Returns the state of this solver.
  virtual const State &getState() const = 0;
  
  /// Returns the objective function associated with this solver.
  const Function &getObjectiveFunction() const;
  
  /// Does this solver have a built-in stopping criterion, or are custom stopping criteria allowed.
  virtual bool hasBuiltInStoppingCriterion() const = 0;
  
  /// Takes one iteration step and returns the status of the algorithm.
  IterationStatus iterate();
  
  boost::shared_ptr< Solver::Results > solve(Function &objFunc,
                                             const boost::numeric::ublas::vector< double > &x0,
                                             const StoppingCriterion &stopCrit,
                                             const Setup &solverSetup = DefaultSetup(),
                                             const Constraints &C = NoConstraints(),
                                             bool timeTest = false);
 protected:
  virtual void doSetup_(const Function &objFunc,
                        const boost::numeric::ublas::vector< double > &x0,
                        const Setup &solverSetup = DefaultSetup(),
                        const Constraints &C = NoConstraints());
  State &getState_();
  virtual IterationStatus iterate_() = 0;
};

#define NATIVESOLVER_H

#endif
