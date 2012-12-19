
#include "SymbolicGradientEvaluator.h"

#include <matheval.h>

using namespace boost::numeric::ublas;

SymbolicGradientEvaluator::SymbolicGradientEvaluator(SymbolicFunctionEvaluator *fe) : 
  refCount_(new int())
{
  fEval_ = fe;
  n_ = fe->getN();
  gEval_ = new void*[fe->getN()];
  
  for(int i = 0; i < fe->getN(); i++)
    gEval_[i] = evaluator_derivative(const_cast< void * >(fe->getEvaluator()),
                                     const_cast< char * >(fe->getVarNames()[i]));
  
  *refCount_ = 1;
  // TODO: this should increment the reference count of the function evaluator!
}

SymbolicGradientEvaluator::SymbolicGradientEvaluator(const SymbolicGradientEvaluator &eval)
{
  fEval_ = eval.fEval_;
  gEval_ = eval.gEval_;
  n_ = eval.n_;
  refCount_ = eval.refCount_;
  (*refCount_)++;
}

SymbolicGradientEvaluator::~SymbolicGradientEvaluator()
{
  if(--*refCount_ <= 0)
  {
    for(int i = 0; i < n_; i++)
      evaluator_destroy(gEval_[i]);
    delete[] gEval_;
    delete refCount_;
  }
}

SymbolicGradientEvaluator &SymbolicGradientEvaluator::operator=(const SymbolicGradientEvaluator &eval)
{
  fEval_ = eval.fEval_;
  gEval_ = eval.gEval_;
  n_ = eval.n_;
  refCount_ = eval.refCount_;
  (*refCount_)++;
  
  return *this;
}

const void *SymbolicGradientEvaluator::getEvalPtr(int i)
{
  return gEval_[i];
}

double *SymbolicGradientEvaluator::eval_(const double *x, double *g) const
{
  for(int i = 0; i < getN(); i++)
    g[i] = evaluator_evaluate(gEval_[i], getN(),
      const_cast< char ** >(fEval_->getVarNames()), (double *)x);
  
  return g;
}
