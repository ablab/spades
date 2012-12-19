
#ifndef GRADIENTEVALUATOR_H

#include <otkpp/objfunc/EvalCountable.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

/// Defines the interface for gradient evaluators.
class GradientEvaluator : public EvalCountable
{
 public:
  /// Default constructor.
  GradientEvaluator() { }
  
  /// Destructor.
  virtual ~GradientEvaluator() { }
  
  /// Creates a deep copy.
  virtual GradientEvaluator *clone() const = 0;
  
  /// Evaluates the gradient.
  /**
   * Evaluates the gradient at a given point.
   * @param x the point to evaluate the gradient at.
   * @param g a gradient vector, where the result is assigned to.
   * @return the computed gradient vector.
   */
  double *operator()(const double *x, double *g) const;
  
  /// Evaluates the gradient.
  /**
   * Evaluates the gradient at a given point.
   * @param x the point to evaluate the gradient at.
   * @param g a gradient vector, where the result is assigned to.
   * @return the computed gradient vector.
   */
  boost::numeric::ublas::vector< double > &operator()(const boost::numeric::ublas::vector< double > &x,
                               boost::numeric::ublas::vector< double > &g) const;
  
  /// Returns the number of variables.
  virtual int getN() const = 0;
  
  /// Returns true, if finite difference approximation is used, false if not.
  virtual bool usesFiniteDifference() const = 0;
 protected:
  virtual double *eval_(const double *x,
                        double *g) const = 0;
};

#define GRADIENTEVALUATOR_H

#endif
