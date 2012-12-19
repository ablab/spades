
#include <matheval.h>
#include <stdexcept>

#include "SymbolicFunctionEvaluator.h"

using namespace boost::numeric::ublas;

SymbolicFunctionEvaluator::SymbolicFunctionEvaluator(const std::string &expr) : 
  refCount_(new int())
{
  expression_ = expr;
  fEval_ = evaluator_create((char *)expr.c_str());
  if(fEval_ == NULL)
    throw std::invalid_argument("Invalid function expression");
  evaluator_get_variables(fEval_, &varNames_, &n_);
  *refCount_ = 1;
}

SymbolicFunctionEvaluator::SymbolicFunctionEvaluator(const SymbolicFunctionEvaluator &eval)
{
  fEval_ = eval.fEval_;
  evaluator_get_variables(fEval_, &varNames_, &n_);
  refCount_ = eval.refCount_;
  (*refCount_)++;
}

SymbolicFunctionEvaluator::~SymbolicFunctionEvaluator()
{
  if(--*refCount_ <= 0)
  {
    if(fEval_ != NULL)
      evaluator_destroy(fEval_);
    delete refCount_;
  }
}

SymbolicFunctionEvaluator &SymbolicFunctionEvaluator::operator=(const SymbolicFunctionEvaluator &eval)
{
  fEval_ = eval.fEval_;
  evaluator_get_variables(fEval_, &varNames_, &n_);
  refCount_ = eval.refCount_;
  (*refCount_)++;
  
  return *this;
}

double SymbolicFunctionEvaluator::eval_(const double *x) const
{
  return evaluator_evaluate(fEval_, n_, varNames_, (double *)x);
}

std::string SymbolicFunctionEvaluator::getExpression() const
{
  return expression_;
}

int SymbolicFunctionEvaluator::getN() const
{
  return n_;
}

bool SymbolicFunctionEvaluator::hasSymbolicExpression() const
{
  return true;
}
