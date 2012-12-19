
#include "CompiledFunctionEvaluator.h"
#include "FDiffGradientEvaluator.h"
#include "FDiffHessianEvaluator.h"
#include "Function.h"
//#include "ProjectedFunctionEvaluator.h"
#ifdef WITH_LIBMATHEVAL
#include "SymbolicFunctionEvaluator.h"
#include "SymbolicGradientEvaluator.h"
#include "SymbolicHessianEvaluator.h"
#endif

#include <sstream>

using namespace boost::numeric::ublas;

#ifdef WITH_LIBMATHEVAL
Function::Function(const std::string &expr,
                   DerivEvalType gEvalType)
{
  constructSymbolicFunction_(expr, gEvalType);
}
#endif

Function::Function(const FunctionEvaluator &fEval,
                   DerivEvalType gEvalType)
{
  eval_.f = fEval.clone();
  
#ifdef WITH_LIBMATHEVAL
  if(gEvalType == Function::DERIV_SYMBOLIC)
  {
    if(dynamic_cast< SymbolicFunctionEvaluator * >(eval_.f) != NULL)
    {
      eval_.g = new SymbolicGradientEvaluator(
        dynamic_cast< SymbolicFunctionEvaluator * >(eval_.f));
      eval_.H = new SymbolicHessianEvaluator(
        dynamic_cast< SymbolicFunctionEvaluator * >(eval_.f),
        dynamic_cast< SymbolicGradientEvaluator * >(eval_.g));
    }
    else
      throw std::runtime_error("symbolic gradient evaluator requires a symbolic function evaluator");
  }
  else
#endif
  if(gEvalType == Function::DERIV_FDIFF_CENTRAL_2)
  {
    eval_.g = new FDiffGradientEvaluator(FDiffGradientEvaluator::CENTRAL_2, eval_.f);
    eval_.H = new FDiffHessianEvaluator(eval_.f);
  }
  else
    throw std::runtime_error("unsupported gradient evaluator type");
}

double Function::operator()(const double *x) const
{
  return (*eval_.f)(x);
}

double Function::operator()(const vector< double > &x) const
{
  return (*eval_.f)(x);
}

#ifdef WITH_LIBMATHEVAL
void Function::constructSymbolicFunction_(const std::string &expr,
                                          DerivEvalType gEvalType)
{
  eval_.f  = new SymbolicFunctionEvaluator(expr);
  if(gEvalType == Function::DERIV_SYMBOLIC)
  {
    eval_.g = new SymbolicGradientEvaluator(
      dynamic_cast< SymbolicFunctionEvaluator * >(eval_.f));
    eval_.H = new SymbolicHessianEvaluator(
      dynamic_cast< SymbolicFunctionEvaluator * >(eval_.f),
      dynamic_cast< SymbolicGradientEvaluator * >(eval_.g));
  }
  else if(gEvalType == Function::DERIV_FDIFF_CENTRAL_2)
  {
    eval_.g = new FDiffGradientEvaluator(FDiffGradientEvaluator::CENTRAL_2, eval_.f);
    eval_.H = new FDiffHessianEvaluator(eval_.f);
  }
  else
    throw std::runtime_error("unsupported gradient evaluator type");
}
#endif

double *Function::g(const double *x, double *g) const
{
  (*eval_.g)(x, g);
  return g;
}

vector< double > &Function::g(const vector< double > &x, 
                              vector< double > &g) const
{
  (*eval_.g)(x, g);
  return g;
}

matrix< double > &Function::H(const vector< double > &x, 
                              matrix< double > &H) const
{
  (*eval_.H)(x, H);
  return H;
}

Function Function::createCopy(DerivEvalType gEvalType) const
{
  Function f;
  f.eval_.f = eval_.f->clone();
#ifdef WITH_LIBMATHEVAL
  if(gEvalType == Function::DERIV_SYMBOLIC)
    f.eval_.g = new SymbolicGradientEvaluator(
      dynamic_cast< SymbolicFunctionEvaluator * >(eval_.f));
  else
#endif
  if(gEvalType == Function::DERIV_FDIFF_CENTRAL_2)
    f.eval_.g = new FDiffGradientEvaluator(
      FDiffGradientEvaluator::CENTRAL_2, eval_.f);
  else
    throw std::runtime_error("unsupported gradient evaluator type");
  f.eval_.H = new FDiffHessianEvaluator(eval_.f);
  
  return f;
}

const FunctionEvaluator &Function::getEvaluator() const
{
  return *eval_.f;
}

int Function::getN() const
{
  return eval_.f->getN();
}

#ifdef WITH_LIBMATHEVAL
const std::string Function::getSymbolicExpression() const
{
  if(hasSymbolicExpression())
    return dynamic_cast< SymbolicFunctionEvaluator * >(eval_.f)->getExpression();
}
#endif

std::vector< std::string > Function::getVariableNames() const
{
  std::vector< std::string > result;
  
#ifdef WITH_LIBMATHEVAL
  if(hasSymbolicExpression())
  {
    result.clear();
    const char **vn = dynamic_cast< 
      SymbolicFunctionEvaluator * >(eval_.f)->getVarNames();
    for(int i = 0; i < getN(); i++)
      result.push_back(vn[i]);
  }
  else
  {
#endif
    result.clear();
    for(int i = 0; i < getN(); i++)
    {
      std::ostringstream vn;
      vn<<"x"<<i+1;
      result.push_back(vn.str());
    }
#ifdef WITH_LIBMATHEVAL
  }
#endif
  
  return result;
}

bool Function::hasSymbolicExpression() const
{
#ifdef WITH_LIBMATHEVAL
  if(dynamic_cast< SymbolicFunctionEvaluator * >(eval_.f) != NULL)
    return true;
  else
#endif
    return false;
}

int Function::getFuncEvalCounter() const
{
  return eval_.f->getEvalCounter();
}

int Function::getGradEvalCounter() const
{
  if(eval_.g != NULL)
    return eval_.g->getEvalCounter();
  else
    return 0;
}

int Function::getHessEvalCounter() const
{
  if(eval_.H != NULL)
    return eval_.H->getEvalCounter();
  else
    return 0;
}

void Function::enableEvalCounting() const
{
  eval_.f->enableEvalCounting();
  if(eval_.g != NULL)
    eval_.g->enableEvalCounting();
  //HEvaluator_->enableEvalCounting();
}

void Function::disableEvalCounting() const
{
  eval_.f->disableEvalCounting();
  if(eval_.g != NULL)
    eval_.g->disableEvalCounting();
  //HEvaluator_->disableEvalCounting();
}

void Function::resetEvalCounters()
{
  eval_.f->resetEvalCounter();
  if(eval_.g != NULL)
    eval_.g->resetEvalCounter();
  //HEvaluator_->resetEvalCounter();
}

Function::Evaluators::Evaluators()
{
  f = NULL;
  g = NULL;
  H = NULL;
}

Function::Evaluators::Evaluators(const Evaluators &e)
{
  if(e.f != NULL)
    f = e.f->clone();
  else
    f = NULL;
  if(e.g != NULL)
  {
    g = e.g->clone();
    if(g->usesFiniteDifference())
      dynamic_cast< FDiffGradientEvaluator * >(g)->setFEvaluator(f);
  }
  else
    g = NULL;
  if(e.H != NULL)
    H = e.H->clone();
  else
    H = NULL;
  // TODO: do this to Hessian also
  /*if(gEvaluator_->usesFiniteDifference())
    dynamic_cast< FDiffGradientEvaluator * >(gEvaluator_)->setFEvaluator(evaluator_);*/
}

Function::Evaluators &Function::Evaluators::operator=(const Evaluators &e)
{
  Function::Evaluators e_(e);
  
  FunctionEvaluator *oldF = f;
  f  = e_.f;
  GradientEvaluator *oldG = g;
  g = e_.g;
  HessianEvaluator *oldH = H;
  H = e_.H;
  
  e_.f = oldF;
  e_.g = oldG;
  e_.H = oldH;
  
  return *this;
}

Function::Evaluators::~Evaluators()
{
  delete H;
  delete g;
  delete f;
}
