
#include "NativeSolver.h"
#include "StoppingCriterion.h"

#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <typeinfo>

#include <time.h>

using namespace boost::numeric::ublas;

static const unsigned int MAX_NUM_ITER = 50000;
static const unsigned long long MIN_TOTAL_TIME = 1e9;

static unsigned long long getTime()
{
  /*timespec tp;
  clock_gettime(CLOCK_REALTIME, &tp);
  return 1e9 * tp.tv_sec + tp.tv_nsec;*/
  return 0;
}

/*NativeSolver::State *NativeSolver::State::clone() const
{
  return new NativeSolver::State(*this);
}*/

/*NativeSolver::NativeSolver(const NativeSolver &s) : baseState_(s.baseState_->clone()) { }

NativeSolver &NativeSolver::operator=(const NativeSolver &s)
{
  baseState_ = std::auto_ptr< NativeSolver::State >(s.baseState_->clone());
}*/

const vector< double > NativeSolver::getX() const
{
  return getState().x;
}

const matrix< double > NativeSolver::getXArray() const
{
  matrix< double > X(setup_->n, 1);
  matrix_column< matrix< double > >(X, 0) = getState().x;
  
  return X;
}

double NativeSolver::getFVal() const
{
  return getState().fx;
}

const vector< double > NativeSolver::getGradient() const
{
  throw std::runtime_error("the algorithm does not use gradient information");
}

const matrix< double > NativeSolver::getHessian() const
{
  throw std::runtime_error("the algorithm does not use Hessian information");
}

unsigned int NativeSolver::getNumIter() const
{
  return getState().nIter;
}

unsigned int NativeSolver::getNumFuncEval() const
{
  return setup_->f.getFuncEvalCounter();
}

unsigned int NativeSolver::getNumGradEval() const
{
  return setup_->f.getGradEvalCounter();
}

unsigned int NativeSolver::getNumHessEval() const
{
  return setup_->f.getHessEvalCounter();
}

NativeSolver::State &NativeSolver::getState_()
{
  return const_cast< NativeSolver::State & >(getState());
}

NativeSolver::IterationStatus NativeSolver::iterate()
{
  NativeSolver::IterationStatus status = iterate_();
  getState_().nIter++;
  return status;
}

boost::shared_ptr< Solver::Results > NativeSolver::solve(Function &objFunc,
                                                         const vector< double > &x0,
                                                         const StoppingCriterion &stopCrit,
                                                         const Solver::Setup &solverSetup,
                                                         const Constraints &C,
                                                         bool timeTest)
{
  bool converged = false;
  unsigned int k;
  unsigned int numRuns = 0;
  NativeSolver::Results *results;
  unsigned long long startTime = 0;
  IterationStatus status = ITERATION_CONTINUE;
  unsigned long long totalTime = 0;
  
  results = new NativeSolver::Results();
  
  do
  {
    k = 0;
    converged = false;
    status = NativeSolver::ITERATION_CONTINUE;
    
    setup(objFunc, x0, solverSetup, C);
    
    if(timeTest == true)
      startTime = getTime();
    else
      objFunc.resetEvalCounters();
    
    do
    {
      if(timeTest == false)
      {
        results->states.push_back(
          boost::shared_ptr< NativeSolver::State >(getState().clone()));
        results->states.back()->X = getXArray();
      }
      if(status == ITERATION_CONTINUE && converged)
        break;
      
      status = iterate();
      
      if(status != ITERATION_CONTINUE && stopCrit.test(*this) == false)
      {
        converged = false;
        break;
      }
      k++;
      
      converged = stopCrit.test(*this);
    }
    while(status == ITERATION_CONTINUE && 
          converged == false && k < MAX_NUM_ITER);
    
    if(timeTest == false)
    {
      results->states.push_back(
        boost::shared_ptr< NativeSolver::State >(getState().clone()));
      results->states.back()->X = getXArray();
    }
    else
    {
      totalTime += getTime() - startTime;
      numRuns++;
    }
  }
  while(timeTest == true && totalTime < MIN_TOTAL_TIME);
  
  results->converged   = converged;
  results->xMin        = getX();
  results->fMin        = getFVal();
  results->numIter     = std::min(k + 1, MAX_NUM_ITER);
  results->numFuncEval = getNumFuncEval();
  results->numGradEval = getNumGradEval();
  //results.numHessEval = getNumHessEval();
  results->setup = 
    boost::shared_ptr< Solver::Setup >(solverSetup.clone());
  results->termVal     = stopCrit.getTestValue(*this);
  
  results->setup->C = boost::shared_ptr< Constraints >(C.clone());
  results->setup->f = objFunc;
  results->setup->m = getM();
  results->setup->n = objFunc.getN();
  
  if(timeTest == true)
  {
    if(converged)
      results->time = ((double)totalTime) / numRuns / 1000.0;
    else
      results->time = 0;
  }
  else
    results->time = 0;
  
  return boost::shared_ptr< NativeSolver::Results >(results);
}

/*void NativeSolver::allocateState_()
{
  state_ = new NativeSolver::State();
}*/

void NativeSolver::doSetup_(const Function &objFunc,
                            const vector< double > &x0,
                            const Solver::Setup &solverSetup,
                            const Constraints &C)
{
  getState_().x = x0;
  getState_().fx = objFunc(x0);
  getState_().nIter = 0;
}
