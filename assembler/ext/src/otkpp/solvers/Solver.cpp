
#include "Solver.h"

#include <typeinfo>

using namespace boost::numeric::ublas;

void Solver::setup(const Function &objFunc,
                   const vector< double > &x0,
                   const Solver::Setup &solverSetup,
                   const Constraints &C)
{
  if(typeid(solverSetup) != typeid(const Solver::DefaultSetup &) && 
     !solverSetup.isCompatibleWith(*this))
    throw std::invalid_argument("the solver and its setup are incompatible");
  
  if(objFunc.getN() != x0.size())
    throw std::runtime_error("dimension mismatch: f!=x0");
  
  setup_ = boost::shared_ptr< Solver::Setup >(solverSetup.clone());
  setup_->n = objFunc.getN();
  
  setup_->f = objFunc;
  setup_->f.resetEvalCounters();
  setup_->f.enableEvalCounting();
  
  doSetup_(objFunc, x0, solverSetup, C);
}
