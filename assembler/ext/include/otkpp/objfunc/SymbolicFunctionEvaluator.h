
#ifndef SYMBOLICFUNCTIONEVALUATOR_H

#include <otkpp/lib/Cloneable.h>
#include <otkpp/objfunc/FunctionEvaluator.h>

#include <boost/numeric/ublas/vector.hpp>
#include <list>
#include <string>

/// Implements a symbolic function evaluator.
/**
 * This class implements a symbolic function evaluator using 
 * GNU libmatheval. Allocated evaluators are reference-counted.
 */
class SymbolicFunctionEvaluator : public Cloneable< SymbolicFunctionEvaluator, FunctionEvaluator >
{
 public:
  /// Constructs a new symbolic function evaluator from the given function expression.
  SymbolicFunctionEvaluator(const std::string &expr);
  
  /// Constructs a new evaluator as a copy of an existing symbolic function evaluator.
  SymbolicFunctionEvaluator(const SymbolicFunctionEvaluator &eval);
  
  /// Destroys this evaluator and deallocates memory associated with it.
  virtual ~SymbolicFunctionEvaluator();
  
  /// Assigns an existing evaluator to this evaluator.
  SymbolicFunctionEvaluator &operator=(const SymbolicFunctionEvaluator &eval);
  
  /// Returns the symbolic function expression of this evaluator.
  std::string getExpression() const;
  
  /// Returns the libmatheval object allocated by this evaluator.
  const void *getEvaluator() const { return fEval_; }
  
  /// Returns the variable names parsed by libmatheval.
  const char **getVarNames() const { return (const char **)varNames_; }
  
  /// Returns the number of variables parsed by libmatheval.
  int getN() const;
  
  bool hasSymbolicExpression() const;
 private:
  std::string expression_;
  void *fEval_;
  int n_;
  int *refCount_;
  char **varNames_;
  
  double eval_(const double *x) const;
};

#define SYMBOLICFUNCTIONEVALUATOR_H

#endif
