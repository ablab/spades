
#ifndef FUNCTIONEVALUATOR_H

#include <otkpp/objfunc/EvalCountable.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <list>

/// Defines the interface for function value evaluators.
class FunctionEvaluator : public EvalCountable
{
 public:
  /// Default constructor.
  FunctionEvaluator() { }
  
  /// Destructor.
  virtual ~FunctionEvaluator() { }
  
  /// Evaluates the function value at the given point.
  /// Evaluates this function at a given point.
  /**
   * Evaluates the function value at a given point.
   * @param x an array of size n, where n is the dimension of this function
   */
  double operator()(const double *x) const;
  
  /// Evaluates the function value at the given point.
  /**
   * Evaluates the function value at the given point.
   * @param x an array of size n, where n is the dimension of this function
   */
  double operator()(const boost::numeric::ublas::vector< double > &x) const;
  
  /// Creates a deep copy.
  virtual FunctionEvaluator *clone() const = 0;
  
  /// Returns true if this function evaluator has a symbolic expression, false if not.
  virtual bool hasSymbolicExpression() const;
  
  /// Returns the number of variables.
  virtual int getN() const = 0;
 protected:
  virtual double eval_(const double *x) const = 0;
};

#define FUNCTIONEVALUATOR_H

#endif
