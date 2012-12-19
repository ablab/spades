
#include "TrustRegionSolver.h"

using namespace boost::numeric::ublas;

double TrustRegionSolver::computeReduction_(double fx,
                                            double fxPlus,
                                            const vector< double > &x,
                                            const vector< double > &g,
                                            const matrix< double > &H,
                                            const vector< double > &p)
{
  double ared = fx - fxPlus;
  double pred = -inner_prod(g, p) - 0.5 * inner_prod(p, prod(H, p));
  
  return ared / pred;
}

double TrustRegionSolver::getDelta() const
{
  return delta_;
}

void TrustRegionSolver::setDelta(double delta)
{
  delta_ = delta;
}
  
void TrustRegionSolver::setup(const Function &f)
{
  f_ = &f;
  doSetup_();
}

bool TrustRegionSolver::updateRadius_(double fx,
                                      double fxPlus,
                                      const vector< double > &x,
                                      const vector< double > &g,
                                      const matrix< double > &H,
                                      const vector< double > &p)
{
  double rho = computeReduction_(fx, fxPlus, x, g, H, p);
  if(rho < 0.25)
    delta_ = 0.25 * norm_2(p);
  else
  {
    if(rho > 0.75 && fabs(norm_2(p) - delta_) < 1e-14)
      delta_ = std::min(2.0 * delta_, deltaMax_);
  }
  if(rho <= eta_)
    return false;
  else
    return true;
}
