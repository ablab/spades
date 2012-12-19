
#include "DoglegSolver.h"

#include <boost/numeric/ublas/lu.hpp>

using namespace boost::numeric::ublas;

DoglegSolver::DoglegSolver()
{
  deltaMax_ = 1e6;
  eps_      = 1e-14;
  eta_      = 0.1;
}

vector< double > &DoglegSolver::computeStep(const vector< double > &x,
                                            double fx,
                                            const vector< double > &g,
                                            const matrix< double > &H,
                                            vector< double > &p,
                                            bool &nonzeroStep,
                                            vector< double > &xPlus,
                                            double &fxPlus,
                                            const vector< double > *Hg)
{
  const int n = x.size();
  
  int j;
  vector< double > pB, pU;
  permutation_matrix< double > P(n);
  matrix< double > S(n, n);
  //vector< double > Hg(n);
  
  //Hg = prod(H, g);
  
  j = 0;
  do
  {
    // compute the Newton point
    S.assign(H);
    lu_factorize(S, P);
    pB = -g;
    lu_substitute(S, P, pB);
    
    // compute the minimizer of the quadratic 
    // model along the negative gradient
    pU = -inner_prod(g, g) / inner_prod(g, prod(H, g)) * g;
    
    if(norm_2(pB) <= delta_)
      p = pB;
    else
    {
      if(norm_2(pU) > delta_)
        p = delta_ * pU / norm_2(pU);
      else
      {
        double tau = computeTau_(pB, pU);
        p = pU + (tau-1.0) * (pB - pU);
      }
    }
    
    xPlus = x + p;
    fxPlus = (*f_)(xPlus);
    
    nonzeroStep = updateRadius_(fx, fxPlus, x, g, H, p);
    
    j++;
  }
  while(!nonzeroStep && j < 100);
  
  if(!nonzeroStep)
  {
    xPlus = x;
    fxPlus = fx;
  }
  
  return p;
}

double DoglegSolver::computeTau_(const vector< double > &pB,
                                 const vector< double > &pU)
{
  vector< double > pB_pU = pB - pU;
  vector< double > v = 2.0*pU-pB;
  
  double a = inner_prod(pB_pU, pB_pU);
  double b = 2.0 * (inner_prod(v, pB_pU));
  double c = inner_prod(v, v) - delta_*delta_;
  
  return (-b + sqrt(b*b-4.0*a*c)) / (2.0*a);
}

void DoglegSolver::doSetup_()
{
  delta_ = 1.0;
}
