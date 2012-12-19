#include "OTK.h"

/*std::list< std::string > OTK::getAlgorithmList()
{
  std::list< std::string > l;
  
  l.push_back("gsl_conjugate_fr");
  l.push_back("gsl_conjugate_pr");
  l.push_back("gsl_nmsimplex");
  l.push_back("gsl_steepest_descent");
  l.push_back("gsl_vector_bfgs");
  l.push_back("gsl_vector_bfgs2");
  
  return l;
}

Solver *OTK::getSolverInstance(const std::string &algoName)
{
  if(algoName == "gsl_conjugate_fr")
    return new GSLFDFSolver(gsl_multimin_fdfminimizer_conjugate_fr);
  else if(algoName == "gsl_conjugate_pr")
    return new GSLFDFSolver(gsl_multimin_fdfminimizer_conjugate_pr);
  else if(algoName == "gsl_steepest_descent")
    return new GSLFDFSolver(gsl_multimin_fdfminimizer_steepest_descent);
  else if(algoName == "gsl_vector_bfgs")
    return new GSLFDFSolver(gsl_multimin_fdfminimizer_vector_bfgs);
  else if(algoName == "gsl_vector_bfgs2")
    return new GSLFDFSolver(gsl_multimin_fdfminimizer_vector_bfgs2);
  else
    throw std::invalid_argument("Unknown minimization algorithm");
}

std::list< std::string > OTK::getStopCritList()
{
  std::list< std::string > l;
  
  l.push_back("f-f*      < eps");
  l.push_back("|nabla f| < eps");
  l.push_back("|x-x*|    < eps");
  
  return l;
}*/
