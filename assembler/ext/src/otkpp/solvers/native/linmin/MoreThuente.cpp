
#include "MoreThuente.h"
#include "OTK.h"
#include "PolyFit.h"

#include <typeinfo>

using namespace boost::numeric::ublas;

MoreThuente::Setup::Setup(double eta, double mu, double gamma, double chi) : 
  eta(eta), mu(mu), gamma(gamma), chi(chi) { }

bool MoreThuente::Setup::isCompatibleWith(const LineMinimizer &s) const
{
  return (typeid(s) == typeid(const MoreThuente &));
}

MoreThuente::MoreThuente() : gamma_(0.0) { }

int MoreThuente::minimize(const vector< double > &x,
                          const vector< double > &d,
                          double alpha0,
                          double fx,
                          const vector< double > &gx,
                          double &alpha,
                          vector< double > &x_plus,
                          double &f_plus,
                          vector< double > &g_plus)
{
  double alphat_plus;
  double phi0, dphi0;
  int i = 0;
  int terminated = 0;
  double gNorm2;
  int suff_desc;
  int stage = 1;
  
  gt_ = vector< double >(f_->getN());
  gu_ = vector< double >(f_->getN());
  
  alphal_ = 0.0;
  alphau_ = 1e3;
  alphat_ = alpha0;
  
  alphau_evaluated_ = 0;
  
  phi0 = fx;
  dphi0 = inner_prod(gx, d);
  
  phil_ = fx;
  dphil_ = dphi0;
  
  while(1)
  {
    xt_ = x + alphat_ * d;
    
    phit_ = (*f_)(xt_);
    f_->g(xt_, gt_);
    dphit_ = inner_prod(gt_, d);
    
    if(stage == 1)
    {
      psil_ = psi(phil_, alphal_, phi0, dphi0, setup_.mu);
      psit_ = psi(phit_, alphat_, phi0, dphi0, setup_.mu);
      psiu_ = psi(phiu_, alphau_, phi0, dphi0, setup_.mu);
      
      dpsil_ = dphil_ - setup_.mu * dphi0;
      dpsit_ = dphit_ - setup_.mu * dphi0;
      dpsiu_ = dphiu_ - setup_.mu * dphi0;
    }
    
    if(i == 50 || fabs(alphal_ - alphau_) < 10.0 * OTK::EPS)
    {
      terminated = 1;
      break;
    }
    
    /* test the sufficient descent condition */
    if(dphit_ <= 0.0 || gamma_ == 0.0)
      suff_desc = 1;
    else
    {
      gNorm2 = inner_prod(gt_, gt_);
      if(-gNorm2 + gamma_ * dphit_ < -setup_.chi * gNorm2)
        suff_desc = 1;
      else
        suff_desc = 0;
    }
    
    /* test the strong Wolfe conditions */
    if(phit_ <= phi0 + setup_.mu * alphat_ * dphi0 && 
       fabs(dphit_) <= setup_.eta * fabs(dphi0) && suff_desc)
      break;
    
    /* Step 1.:trial value selection */
    if(stage == 1)
      alphat_plus = trialstep(alphat_, alphal_, alphau_,
                              psil_, dpsil_, psit_, dpsit_,
                              psiu_, dpsiu_,
                              *f_, x, d);
    else
      alphat_plus = trialstep(alphat_, alphal_, alphau_,
                              phil_, dphil_, phit_, dphit_,
                              phiu_, dphiu_,
                              *f_, x, d);
    
    assert(!std::isnan(alphat_plus));
    
    if(psit_ <= 0.0 && dphit_ > 0.0)
      stage = 2;
    
    /* Step 2.:bracketing */
    if(stage == 1)
      bracket(alphal_, alphat_, alphau_,
              psil_, psit_, dpsit_,
              phil_, phit_, phiu_,
              dphil_, dphit_, dphiu_);
    else
      bracket(alphal_, alphat_, alphau_,
              phil_, phit_, dphit_,
              phil_, phit_, phiu_,
              dphil_, dphit_, dphiu_);
    
    alphat_ = alphat_plus;
    
    i++;
  }
  
  alpha = alphat_;
  x_plus = xt_;
  f_plus = phit_;
  g_plus = gt_;
  
  if(terminated)
    return OTK::DOMAIN_ERROR;
  else
    return OTK::SUCCESS;
}

void MoreThuente::setGamma(double gamma)
{
  gamma_ = gamma;
}

double MoreThuente::trialstep(double alphat, double alphal, double alphau,
                              double fl, double dfl,
                              double ft, double dft,
                              double fu, double dfu,
                              const Function &f,
                              const vector< double > &x,
                              const vector< double > &d)
{
  double alphac, alphaq, alphas, alphae;
  double alphat_plus;
  int status1, status2;
  
  if(ft > fl)
  {
    status1 = cubic_minimizer(alphal, fl, dfl, alphat, ft, dft, 0, &alphac);
    status2 = quad_minimizer1(alphal, fl, dfl, alphat, ft, &alphaq);
    
    if(status1 == OTK::SUCCESS && status2 == OTK::SUCCESS)
    {
      if(fabs(alphac - alphal) < fabs(alphaq - alphal))
        alphat_plus = alphac;
      else
        alphat_plus = 0.5 * (alphaq + alphac);
      
      assert(!std::isnan(alphat_plus));
    }
    else
      goto bisection1;
  }
  else if(dft * dfl < 0.0)
  {
    status1 = cubic_minimizer(alphal, fl, dfl, alphat, ft, dft, 0, &alphac);
    status2 = quad_minimizer2(alphal, dfl, alphat, dft, &alphas);
    
    if(status1 == OTK::SUCCESS && status2 == OTK::SUCCESS)
    {
      if(fabs(alphac - alphat) >= fabs(alphas - alphat))
        alphat_plus = alphac;
      else
        alphat_plus = alphas;
    }
    else
      goto bisection1;
  }
  else if(fabs(dft) <= fabs(dfl))
  {
    status1 = cubic_minimizer(alphal, fl, dfl, alphat, ft, dft, 1, &alphac);
    status2 = quad_minimizer2(alphal, dfl, alphat, dft, &alphas);
    
    if(status1 == OTK::SUCCESS && status2 == OTK::SUCCESS)
    {
      if(fabs(alphac - alphat) < fabs(alphas - alphat))
        alphat_plus = alphac;
      else
        alphat_plus = alphas;
      
      if(alphat > alphal)
        alphat_plus = std::min(alphat + 0.66 * (alphau - alphat), alphat_plus);
      else
        alphat_plus = std::max(alphat + 0.66 * (alphau - alphat), alphat_plus);
    }
    else
      goto bisection2;
  }
  else
  {
    if(!alphau_evaluated_)
    {
      xt_ = x + alphau * d;
      phiu_ = f(xt_);
      f.g(xt_, gu_);
      dphiu_ = inner_prod(gu_, d);
      alphau_evaluated_ = 1;
    }
    
    if(cubic_minimizer(alphat, ft, dft, alphau, fu, dfu, 1, &alphae) == OTK::SUCCESS)
      alphat_plus = alphae;
    else
      goto bisection2;
  }
  
  return alphat_plus;
  
  bisection1:
  //assert(!isnan(alphal) && !isnan(alphat) && !isnan(alphau));
  return alphal + 0.25 * (alphat - alphal);
  bisection2:
  //assert(!isnan(alphal) && !isnan(alphat) && !isnan(alphau));
  return alphat + 0.75 * (alphau - alphat);
}

void MoreThuente::bracket(double &alphal, double alphat, double &alphau,
                          double fl, double ft, double dft,
                          double &phil, double phit, double &phiu,
                          double &dphil, double dphit, double &dphiu)
{
  if(ft > fl)                            /* Case 1. */
  {
    alphau = alphat;
    phiu = phit;
    dphiu = dphit;
  }
  else
  {
    if(dft * (alphat - alphal) < 0.0)    /* Case 2. */
    {
      alphal = alphat;
      phil = phit;
      dphil = dphit;
    }
    else
    {
      alphau = alphal;                   /* Case 3. */
      phiu = phil;
      dphiu = dphil;
      alphal = alphat;
      phil = phit;
      dphil = dphit;
    }
  }
}

double MoreThuente::psi(double phi, double alpha, double phi0, double dphi0, double mu)
{
  return phi - phi0 - mu * dphi0 * alpha;
}

void MoreThuente::doSetup_(const LineMinimizer::Setup &s)
{
  if(typeid(s) != typeid(LineMinimizer::DefaultSetup))
  {
    setup_ = dynamic_cast< const MoreThuente::Setup & >(s);
    gamma_ = setup_.gamma;
  }
  else
  {
    setup_ = MoreThuente::Setup();
    gamma_ = 0.0;
  }
}
