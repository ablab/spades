
#include "Fletcher.h"
#include "OTK.h"
#include "PolyFit.h"

#include <typeinfo>

using namespace boost::numeric::ublas;

static const int I = 20;
static const int J = 10;

Fletcher::Setup::Setup(double eta, double mu, double chi, double tau) : 
  eta(eta), mu(mu), chi(chi), tau(tau) { }

bool Fletcher::Setup::isCompatibleWith(const LineMinimizer &s) const
{
  return (typeid(s) == typeid(const Fletcher &));
}

int Fletcher::minimize(const vector< double > &x,
                       const vector< double > &d,
                       double alpha0,
                       double fx,
                       const vector< double > &gx,
                       double &alpha,
                       vector< double > &x_plus,
                       double &f_plus,
                       vector< double > &g_plus)
{
  alphal_ = 0.0;
  alphau_ = 1e10;
  alphat_ = alpha0;
  
  int i, j;
  int status = OTK::SUCCESS;
  
  double alphat_plus;
  double phi0;
  double dphi0;
  
  gt_ = vector< double >(f_->getN());
  
  phi0 = fx;
  dphi0 = inner_prod(gx, d);
  
  phil_ = fx;
  dphil_ = dphi0;
  
  j = 0;
  do
  {
    xt_ = x + alphat_ * d;
    phit_ = (*f_)(xt_);
    
    i = 0;
    while(phit_ > phi0 + setup_.mu * alphat_ * dphi0)
    {
      status = quad_minimizer1(alphal_, phil_, dphil_,
                               alphat_, phit_, &alphat_plus);
      if(status != OTK::SUCCESS)
        alphat_plus = alphal_ + 0.5 * (alphat_ - alphal_);
      alphat_plus = interp_safequard_(alphat_plus);
      
      alphau_ = alphat_;
      alphat_ = alphat_plus;
      
      xt_ = x + alphat_ * d;
      phit_ = (*f_)(xt_);
      
      i++;
      if(i == I)
        break;
    }
    
    f_->g(xt_, gt_);
    dphit_ = inner_prod(gt_, d);
    
    if(j == J || fabs(alphal_ - alphau_) < 10.0 * OTK::EPS)
      break;
    
    if(dphit_ < setup_.eta * dphi0)
    {
      status = quad_minimizer2(alphal_, dphil_, alphat_, dphit_, &alphat_plus);
      if(status != OTK::SUCCESS)
        alphat_plus = 0.5 * (alphau_ + alphat_);
      else   
        alphat_plus = extrap_safequard_(alphat_plus);
      
      alphal_ = alphat_;
      phil_   = phit_;
      dphil_  = dphit_;
      alphat_ = alphat_plus;
    }
    else
      break;
    
    j++;
  }
  while(1);
  
  alpha = alphat_;
  x_plus = xt_;
  f_plus = phit_;
  g_plus = gt_;
  
  return OTK::SUCCESS;
}

double Fletcher::extrap_safequard_(double alpha)
{
  double delta = alphat_ - alphal_;
  double alpha_min = alphat_ + setup_.tau * delta;
  double alpha_max = alphat_ + setup_.chi * delta;
  
  alpha = std::max(alpha, alpha_min);
  alpha = std::min(alpha, alpha_max);
  
  alpha = std::min(alpha, alphat_ + 0.5 * (alphau_ - alphat_));
  
  return alpha;
}

double Fletcher::interp_safequard_(double alpha)
{
  double delta_alpha = alphau_ - alphal_;
  
  alpha = std::max(alpha, alphal_ + setup_.tau * delta_alpha);
  alpha = std::min(alpha, alphau_ - setup_.tau * delta_alpha);
  
  return alpha;
}

void Fletcher::doSetup_(const LineMinimizer::Setup &s)
{
  if(typeid(s) != typeid(LineMinimizer::DefaultSetup))
    setup_ = dynamic_cast< const Fletcher::Setup & >(s);
  else
    setup_ = Fletcher::Setup();
}
