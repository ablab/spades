
#ifndef SYMBOLICGRADIENTEVALUATOR_H

#include <otkpp/lib/Cloneable.h>
#include <otkpp/objfunc/GradientEvaluator.h>
#include <otkpp/objfunc/SymbolicFunctionEvaluator.h>

/// Implements a symbolic gradient evaluator.
/**
 * This class implements a symbolic gradient evaluator using 
 * GNU libmatheval. Allocated evaluators are reference-counted.
 */
class SymbolicGradientEvaluator : public Cloneable< SymbolicGradientEvaluator, GradientEvaluator >
{
 public:
  /// Constructs a symbolic gradient evaluator by using an existing symbolic function evaluator.
  /**
   * This constructor takes an existing symbolic function evaluator 
   * and creates n partial derivative evaluators, where n is the number 
   * of variables in the function evaluator. The reference count 
   * of the function evaluator is incremented.
   */
  SymbolicGradientEvaluator(SymbolicFunctionEvaluator *fe);
  
  /// Constructs a new evaluator as a copy of an existing symbolic gradient evaluator.
  /**
   * This constructor generates a copy of en existing symbolic 
   * gradient evaluator. The reference count of the source 
   * gradient evaluator is incremented.
   */
  SymbolicGradientEvaluator(const SymbolicGradientEvaluator &eval);
  
  /// Destroys this evaluator and deallocates memory associated with it.
  ~SymbolicGradientEvaluator();
  
  /// Assigns an existing evaluator to this evaluator.
  SymbolicGradientEvaluator &operator=(const SymbolicGradientEvaluator &eval);
  
  /// Returns the libmatheval evaluator for the ith partial derivative.
  const void *getEvalPtr(int i);
  
  /// Returns the number of variables in the function expression.
  int getN() const { return n_; }
  
  bool usesFiniteDifference() const { return false; }
 private:
  SymbolicFunctionEvaluator *fEval_;
  void **gEval_;
  int n_;
  int *refCount_;
  
  double *eval_(const double *x, double *g) const;
};

#define SYMBOLICGRADIENTEVALUATOR_H

#endif
