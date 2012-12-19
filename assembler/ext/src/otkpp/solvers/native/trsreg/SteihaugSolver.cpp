
#include "SteihaugSolver.h"

using namespace boost::numeric::ublas;

SteihaugSolver::SteihaugSolver()
{
  deltaMax_ = 1e8;
  eps_      = 1e-16;
  eta_      = 0.15;
}

vector< double > &SteihaugSolver::computeStep(const vector< double > &x,
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
  
  vector< double > r;
  vector< double > rPlus;
  vector< double > d;
  double dHd;
  double alpha, beta;
  double tau;
  //double rho;
  vector< double > pPlus;
  int j, k;
  double r0_norm = norm_2(g);
  vector< double > Hd;
  double rr;
  double rPlus2;
  
  j = 0;
  do
  {
    if(r0_norm < eps_)
    {
      return p;
      pPlus = p;
      goto end;
    }
    
    r = g;
    d = -r;
    p = zero_vector< double >(n);
    
    k = 0;
    while(k < 2*n)
    {
      Hd = prod(H, d);
      dHd = inner_prod(d, Hd);
      if(dHd <= 0.0)
      {
        tau = computeTau_(d, p);
        pPlus = p + tau*d;
        break;
      }
      rr = inner_prod(r, r);
      alpha = rr / dHd;
      pPlus = p + alpha * d;
      if(norm_2(pPlus) >= delta_)
      {
        tau = computeTau_(d, p);
        pPlus = p + tau*d;
        break;
      }
      rPlus = r + alpha * Hd;
      rPlus2 = inner_prod(rPlus, rPlus);
      if(rPlus2 < (eps_*r0_norm)*(eps_*r0_norm))
        break;
      beta = rPlus2 / rr;
      d = -rPlus + beta * d;
      r = rPlus;
      p = pPlus;
      k++;
    }
    
    end:
    p = pPlus;
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

double SteihaugSolver::computeTau_(const vector< double > &d,
                                   const vector< double > &p)
{
  double pd = inner_prod(p, d);
  double d2 = inner_prod(d, d);
  double p2 = inner_prod(p, p);
  
  return (-pd + sqrt(pd*pd - d2*(p2 - delta_*delta_))) / d2;
}

void SteihaugSolver::doSetup_()
{
  delta_ = 1.0;
}
