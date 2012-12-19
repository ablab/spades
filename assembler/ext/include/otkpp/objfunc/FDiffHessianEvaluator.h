
#ifndef FDIFFHESSIANEVALUATOR_H

#include <otkpp/lib/Cloneable.h>
#include <otkpp/objfunc/FunctionEvaluator.h>
#include <otkpp/objfunc/HessianEvaluator.h>

class FunctionEvaluator;

/// Implements a finite-difference Hessian evaluator.
class FDiffHessianEvaluator : public Cloneable< FDiffHessianEvaluator, HessianEvaluator >
{
 public:
  /// Constructs a new finite-difference Hessian evaluator by using an existing function evaluator.
  FDiffHessianEvaluator(FunctionEvaluator *fe);
  
  /// Returns the number of variables.
  int getN() const { return fEval_->getN(); }
  
  bool usesFiniteDifference() const { return true; }
 private:
  FunctionEvaluator *fEval_;
  mutable boost::numeric::ublas::vector< double > f_hi_;
  mutable boost::numeric::ublas::vector< double > x_hi_;
  mutable boost::numeric::ublas::vector< double > x_hj_;
  mutable boost::numeric::ublas::vector< double > x_hihj_;
  
  boost::numeric::ublas::matrix< double > &eval_(const boost::numeric::ublas::vector< double > &x,
                          boost::numeric::ublas::matrix< double > &H) const;
};

#define FDIFFHESSIANEVALUATOR_H

#endif
