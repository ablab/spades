
#ifndef FUNCTION_H

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <list>
#include <string>

#include <otkpp/objfunc/FunctionEvaluator.h>
#include <otkpp/objfunc/GradientEvaluator.h>
#include <otkpp/objfunc/HessianEvaluator.h>

class FDiffGradientEvaluator;

/// Defines an objective function \f$f:\mathbb R^{n}\rightarrow\mathbb R\f$.
class Function
{
 public:
  enum FuncEvalType
  {
#ifdef WITH_LIBMATHEVAL
    SYMBOLIC,  ///< symbolically evaluated function
#endif
    COMPILED   ///< a precompiled function written in C/C++
  };
  
  enum DerivEvalType
  {
#ifdef WITH_LIBMATHEVAL
    DERIV_SYMBOLIC,        ///< symbolic derivative
#endif
    DERIV_FDIFF_FORWARD,   ///< two-point forward difference approximation
    DERIV_FDIFF_CENTRAL_2, ///< two-point central difference approximation
    DERIV_FDIFF_CENTRAL_4  ///< four-point central difference approximation
  };
  
  Function() { }

#ifdef WITH_LIBMATHEVAL
  /// Constructs a new function with the given expression.
  Function(const std::string &expr,
           DerivEvalType gEvalType = DERIV_SYMBOLIC);
#endif
  
  /// Constructs a new function from an existing evaluator.
  Function(const FunctionEvaluator &fEval,
           DerivEvalType gEvalType);
  
  /*/// Constructs a projection of an existing function.
  Function(const Function &f,
           const std::list< int > &idx,
           const vector< double > &x,
           DerivEvalType gEvalType = DERIV_SYMBOLIC);*/
  
  virtual ~Function() { }
  
  /// Returns a copy of this function.
  /**
   * Returns a copy of this function. Different gradient 
   * and Hessian evaluators can be specified to the new 
   * instance.
   */
  Function createCopy(DerivEvalType gEvalType) const;
  
  // TODO: this is a temporary hack
  const FunctionEvaluator &getEvaluator() const;
  
  /// Returns the number of variables.
  int getN() const;

#ifdef WITH_LIBMATHEVAL
  /// Returns the expression/name of this function (if applicable).
  const std::string getSymbolicExpression() const;
#endif
  
  /// Returns the names of the variables of this function (if a symbolic expression is available)
  std::vector< std::string > getVariableNames() const;
  
  /// Returns true if this function has a symbolic expression, false if not.
  bool hasSymbolicExpression() const;
  
  /// Evaluates this function at a given point.
  /**
   * Evaluates this function at a given point.
   * @param x an array of size n, where n is the dimension of this function
   */
  double operator()(const double *x) const;
  
  /// Evaluates this function at a given point.
  /**
   * Evaluates this function at a given point.
   * @param x an n-dimensional vector, where n is the dimension of this function
   */
  double operator()(const boost::numeric::ublas::vector< double > &x) const;
  
  /// Evaluates the gradient of this function.
  /**
   * Evaluates the gradient of this function at a given point.
   * @param x the point to evaluate this function at.
   * @param g a gradient vector, where the result is assigned to.
   * @return the computed gradient vector.
   */
  double *g(const double *x, double *g) const;
  
  /// Evaluates the gradient of this function.
  /**
   * Evaluates the gradient of this function at a given point. 
   * Each parameter must have the same dimensions as this function. 
   * @param x the point to evaluate this function at. 
   * @param g a gradient vector, where the result is assigned to.
   * @return the computed gradient vector.
   */
  boost::numeric::ublas::vector< double > &g(const boost::numeric::ublas::vector< double > &x, boost::numeric::ublas::vector< double > &g) const;
  
  /// Evaluates the Hessian of this function.
  /**
   * Evaluates the Hessian of this function at a given point x.
   * Each parameter must have the same dimensions as this function.
   * @param x the point to evaluate this function at.
   * @param H an nxn -matrix, where the result is assigned to.
   * @return the computed Hessian matrix.
   */
  boost::numeric::ublas::matrix< double > &H(const boost::numeric::ublas::vector< double > &x, boost::numeric::ublas::matrix< double > &H) const;
  
  /// Returns the number of function evaluations of this function. 
  int getFuncEvalCounter() const;
  
  /// Returns the number of gradient evaluations of this function. 
  int getGradEvalCounter() const;
  
  /// Returns the number of Hessian evaluations of this function. 
  int getHessEvalCounter() const;
  
  /// Enables evaluation counting for this function. 
  void enableEvalCounting() const;
  
  /// Disables evaluation counting for this function. 
  void disableEvalCounting() const;
  
  /// Resets all evaluation counters. 
  void resetEvalCounters();
 protected:
  struct Evaluators
  {
    FunctionEvaluator *f;
    GradientEvaluator *g;
    HessianEvaluator  *H;
    
    Evaluators();
    Evaluators(const Evaluators &e);
    Evaluators &operator=(const Evaluators &e);
    ~Evaluators();
  };
  
  Evaluators eval_;
  
  /*FunctionEvaluator *evaluator_;
  GradientEvaluator *gEvaluator_;
  HessianEvaluator  *HEvaluator_;*/
  
#ifdef WITH_LIBMATHEVAL
  void constructSymbolicFunction_(const std::string &expr,
                                  DerivEvalType gEvalType = DERIV_SYMBOLIC);
#endif
};

#define FUNCTION_H

#endif
