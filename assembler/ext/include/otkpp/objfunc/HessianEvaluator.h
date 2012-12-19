
#ifndef HESSIANEVALUATOR_H

#include <otkpp/objfunc/EvalCountable.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

/// Defines the interface for Hessian evaluators.
class HessianEvaluator : public EvalCountable
{
 public:
  /// Default constructor.
  HessianEvaluator() { }
  
  /// Destructor.
  virtual ~HessianEvaluator() { }
  
  /// Creates a deep copy.
  virtual HessianEvaluator *clone() const = 0;
  
  /// Evaluates the Hessian.
  /**
   * Evaluates the Hessian at a given point x.
   * @param x the point to evaluate this function at.
   * @param H an nxn -matrix, where the result is assigned to.
   * @return the computed Hessian matrix.
   */
  boost::numeric::ublas::matrix< double > &operator()(const boost::numeric::ublas::vector< double > &x,
                               boost::numeric::ublas::matrix< double > &H) const;
  
  /// Returns the number of variables.
  virtual int getN() const = 0;
  
  /// Returns true, if finite difference approximation is used, false if not.
  virtual bool usesFiniteDifference() const = 0;
 protected:
  virtual boost::numeric::ublas::matrix< double > &eval_(const boost::numeric::ublas::vector< double > &x,
                                  boost::numeric::ublas::matrix< double > &H) const = 0;
};

#define HESSIANEVALUATOR_H

#endif
