
#ifndef COMPILEDFUNCTIONEVALUATOR_H

#include "FunctionEvaluator.h"

// DEPRECATED?
class CompiledFunctionEvaluator : public FunctionEvaluator
{
 public:
  CompiledFunctionEvaluator(double (*f)(const boost::numeric::ublas::vector< double > &x));
  
  double f(const boost::numeric::ublas::vector< double > &x) const;
  
  virtual int getN() const = 0;
 protected:
  double (*f_)(const boost::numeric::ublas::vector< double > &x);
};

#define COMPILEDFUNCTIONEVALUATOR_H

#endif
