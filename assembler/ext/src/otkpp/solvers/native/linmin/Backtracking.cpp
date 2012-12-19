
#include "Backtracking.h"
#include "OTK.h"

using namespace boost::numeric::ublas;

Backtracking::Backtracking()
{
  alpha0_ = 1.0;
  mu_ = 0.01;
  rho_ = 2.0;
  sigma_ = 0.5;
}

int Backtracking::minimize(const vector< double > &x,
                           const vector< double > &d,
                           double alpha0,
                           double fx,
                           const vector< double > &gx,
                           double &alpha,
                           vector< double > &x_plus,
                           double &f_plus,
                           vector< double > &g_plus)
{
  alpha = alpha0_;
  double gd = inner_prod(gx, d);
  
  while((f_plus = (*f_)(x + alpha*d)) >= fx + mu_*alpha*gd && alpha > OTK::SQRT_EPS)
    alpha *= sigma_;
  f_plus = (*f_)(x + alpha*d);
  
  if(alpha == alpha0_)
    alpha0_ *= rho_;
  else
    alpha0_ = alpha;
  
  x_plus = x + alpha * d;
  
  return OTK::SUCCESS;
}
