
#ifndef SYMBOLICHESSIANEVALUATOR_H

#include <otkpp/lib/Cloneable.h>
#include <otkpp/objfunc/HessianEvaluator.h>
#include <otkpp/objfunc/SymbolicFunctionEvaluator.h>
#include <otkpp/objfunc/SymbolicGradientEvaluator.h>

/// Implements a symbolic Hessian evaluator.
/**
 * This class implements a symbolic Hessian evaluator using 
 * GNU libmatheval. Allocated evaluators are reference-counted.
 */
class SymbolicHessianEvaluator : public Cloneable< SymbolicHessianEvaluator, HessianEvaluator >
{
 public:
  /// Constructs a symbolic Hessian evaluator by using an existing symbolic gradient evaluator.
  /**
   * This constructor takes an existing symbolic gradient evaluator 
   * and creates nxn second order partial derivative evaluators, 
   * where n is the number of variables in the function evaluator. 
   * The reference count of the gradient evaluator is incremented.
   */
  SymbolicHessianEvaluator(SymbolicFunctionEvaluator *fe, SymbolicGradientEvaluator *ge);
  
  /// Constructs a new evaluator as a copy of an existing symbolic Hessian evaluator.
  /**
   * This constructor generates a copy of en existing symbolic 
   * Hessian evaluator. The reference count of the source Hessian 
   * evaluator is incremented.
   */
  SymbolicHessianEvaluator(const SymbolicHessianEvaluator &eval);
  
  /// Destroys this evaluator and deallocates memory associated with it.
  ~SymbolicHessianEvaluator();
  
  /// Assigns an existing evaluator to this evaluator.
  SymbolicHessianEvaluator &operator=(const SymbolicHessianEvaluator &eval);
  
  /// Returns the number of variables.
  /**
   * Returns the number of variables in the function evaluator 
   * this gradient evaluator is associated with.
   */
  int getN() const { return n_; }
  
  bool usesFiniteDifference() const { return false; }
 private:
  SymbolicFunctionEvaluator *fEval_;
  SymbolicGradientEvaluator *gEval_;
  void ***HEval_;
  int n_;
  int *refCount_;
  
  matrix< double > &eval_(const vector< double > &x,
                          matrix< double > &H) const;
};

#define SYMBOLICHESSIANEVALUATOR_H

#endif
