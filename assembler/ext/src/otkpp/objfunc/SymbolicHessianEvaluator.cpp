
#include "SymbolicHessianEvaluator.h"

#include <matheval.h>

using namespace boost::numeric::ublas;

SymbolicHessianEvaluator::SymbolicHessianEvaluator(SymbolicFunctionEvaluator *fe, SymbolicGradientEvaluator *ge) : 
  refCount_(new int())
{
  int i, j;
  n_ = fe->getN();
  const int n = fe->getN();
  
  fEval_ = fe;
  gEval_ = ge;
  HEval_ = new void**[n];
  
  for(i = 0; i < n; i++)
  {
    HEval_[i] = new void*[n];
    for(j = 0; j <= i; j++)
      HEval_[i][j] = evaluator_derivative(const_cast< void * >(gEval_->getEvalPtr(i)),
        const_cast< char * >(fEval_->getVarNames()[j]));
  }
  
  *refCount_ = 1;
  // TODO: this should increment the reference count of the gradient evaluator!
}

SymbolicHessianEvaluator::SymbolicHessianEvaluator(const SymbolicHessianEvaluator &eval)
{
  fEval_ = eval.fEval_;
  gEval_ = eval.gEval_;
  HEval_ = eval.HEval_;
  n_ = eval.n_;
  refCount_ = eval.refCount_;
  (*refCount_)++;
}

SymbolicHessianEvaluator::~SymbolicHessianEvaluator()
{
  int i, j;
  
  if(--*refCount_ <= 0)
  {
    for(i = 0; i < n_; i++)
      for(j = 0; j <= i; j++)
        if(j <= i)
          evaluator_destroy(HEval_[i][j]);
    
    for(i = 0; i < n_; i++)
      delete[] HEval_[i];
    delete[] HEval_;
    delete refCount_;
  }
}

SymbolicHessianEvaluator &SymbolicHessianEvaluator::operator=(const SymbolicHessianEvaluator &eval)
{
  fEval_ = eval.fEval_;
  gEval_ = eval.gEval_;
  HEval_ = eval.HEval_;
  n_ = eval.n_;
  refCount_ = eval.refCount_;
  (*refCount_)++;
  
  return *this;
}

matrix< double > &SymbolicHessianEvaluator::eval_(const vector< double > &x,
                                                  matrix< double > &H) const
{
  const int n = fEval_->getN();
  int i, j;
  
  for(i = 0; i < n; i++)
  {
    for(j = 0; j <= i; j++)
    {
      H(i, j) = evaluator_evaluate(HEval_[i][j], n,
        const_cast< char ** >(fEval_->getVarNames()),
        const_cast< double * >(&x[0]));
    }
  }
  
  for(i = 0; i < n; i++)
    for(j = i + 1; j < n; j++)
      H(i, j) = H(j, i);
  
  return H;
}
