
#ifndef SOLVER_H

#include <otkpp/constraints/Constraints.h>
#include <otkpp/objfunc/Function.h>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/shared_ptr.hpp>

class Function;
class Solver;
class StoppingCriterion;

/// The base class for all solvers.
class Solver
{
 public:
  /// Defines the parameters of a solver.
  struct Setup
  {
    boost::shared_ptr< Constraints > C;
    Function f;
    int m;
    int n;
    
    virtual Setup *clone() const = 0;
    
    /// Is this solver setup compatible with the given solver.
    virtual bool isCompatibleWith(const Solver &s) const = 0;
    
    virtual ~Setup() {};
  };
  
  /// Defines a default solver setup specifying that default parameters are used.
  struct DefaultSetup : public Cloneable< DefaultSetup, Setup >
  {
    bool isCompatibleWith(const Solver &s) const { return true; }
  };
  
  /// Defines the results of a solver.
  struct Results
  {
    bool converged;                   ///< was the chosen stopping criterion satisfied
    double fMin;                      ///< final function value
    unsigned int numFuncEval;         ///< the number of used function evaluations
    unsigned int numGradEval;         ///< the number of used gradient evaluations
    unsigned int numIter;             ///< the number of used iterations
    boost::shared_ptr< Setup > setup; ///< the used input parameters
    double termVal;                   ///< the final termination test value
    double time;                      ///< used CPU time
    boost::numeric::ublas::vector< double > xMin;            ///< the final iterate
    
    virtual ~Results() { };
  };
  
  virtual ~Solver() { }
  
  /// Returns the number of points this solver produces each iteration (default = 1).
  virtual unsigned int getM() const { return 1; }
  
  /// Returns the number of dimensions of the problem associated with this solver.
  unsigned int getN() const { return setup_->n; }
  
  /// Returns the name of this solver.
  virtual std::string getName() const = 0;
  
  /// Returns the objective function associated with this solver.
  const Function &getObjFunc() const { return setup_->f; }
  
  /// Is this solver implemented natively in C/C++, or does it call some external routine.
  virtual bool isExternalSolver() const { return false; }
  
  /// Initializes this solver.
  /**
   * Initializes this solver with objective function objFunc and starting 
   * point x0. The dimensions of objFunc and x0 must match.
   * @param objFunc objective function
   * @param x0 starting point
   * @param solverSetup a structure containing algorithm-specific parameters
   */
  void setup(const Function &objFunc,
             const boost::numeric::ublas::vector< double > &x0,
             const Setup &solverSetup = DefaultSetup(),
             const Constraints &C = NoConstraints());
  
  /// Solves the given problem by using this solver and returns the results.
  /**
   * @param objFunc objective function
   * @param x0 starting point
   * @param solverSetup solver parameters
   * @param C constraints
   * @param stopCrit stopping criterion
   * @return a SolverResults structure containing the results
   */
  virtual boost::shared_ptr< Solver::Results > solve(
    Function &objFunc,
    const boost::numeric::ublas::vector< double > &x0,
    const StoppingCriterion &stopCrit,
    const Setup &solverSetup = DefaultSetup(),
    const Constraints &C = NoConstraints(),
    bool timeTest = false) = 0;
  
  /// Does this solver support the given constraints.
  virtual bool supportsConstraints(const Constraints &C) { return false; }
  
  /// Does this solver use gradient information.
  virtual bool usesGradient() const = 0;
  
  /// Does this solver use Hessian information.
  virtual bool usesHessian() const = 0;
 protected:
  boost::shared_ptr< Solver::Setup > setup_;
  
  virtual void doSetup_(const Function &objFunc,
                        const boost::numeric::ublas::vector< double > &x0,
                        const Setup &solverSetup,
                        const Constraints &C) = 0;
};

#define SOLVER_H

#endif
