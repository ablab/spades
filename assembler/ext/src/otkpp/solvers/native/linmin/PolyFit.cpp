
#include "OTK.h"
#include "PolyFit.h"

#include <algorithm>
#include <cassert>
#include <cmath>

int cubic_minimizer (double a, double fa, double dfa,
                     double b, double fb, double dfb,
                     int extrapolate, double *alpha_min)
{
  double z, w;
  double denom;
  double D;
  
  if(std::isnan(a) || std::isinf(a) || std::isnan(fa) || std::isinf(fa) || std::isnan(dfa) || std::isinf(dfa) || 
     std::isnan(b) || std::isinf(b) || std::isnan(fb) || std::isinf(fb) || std::isnan(dfb) || std::isinf(dfb))
    return OTK::DOMAIN_ERROR;
  
  z = 3.0 * (fa - fb) / (b - a) + dfa + dfb;
  D = z * z - dfa * dfb;
  
  if(extrapolate == 0)
  {
    if(D < 0.0)
      return OTK::DOMAIN_ERROR;
    
    if(a < b)
      w = sqrt(D);
    else
      w = -sqrt(D);
  }
  else
  {
    if(a < b)
      w = sqrt(std::max(D, 0.0));
    else
      w = -sqrt(std::max(D, 0.0));
  }
  
  denom = dfb - dfa + 2.0 * w;
  
  if(fabs(denom) < OTK::EPS)
    return OTK::ZERO_DIVISION;
  
  *alpha_min = b - (b - a) * (dfb + w - z) / denom;
  /*if(std::isnan(*alpha_min) || std::isinf(*alpha_min))
     return GSL_FAILURE;*/
  assert(!std::isnan(*alpha_min) && !std::isinf(*alpha_min));
  
  return OTK::SUCCESS;
}

int quad_minimizer1 (double a, double fa, double dfa,
                     double b, double fb,
                     double *alpha_min)
{
  if(std::isnan(a) || std::isinf(a) || std::isnan(fa) || std::isinf(fa) || std::isnan(dfa) || std::isinf(dfa) || 
     std::isnan(b) || std::isinf(b) || std::isnan(fb) || std::isinf(fb))
    return OTK::DOMAIN_ERROR;
  
  double denom = 2.0 * (fa - fb + (b - a) * dfa);
  
  if(fabs(denom) < OTK::EPS)
    return OTK::ZERO_DIVISION;
  
  *alpha_min = a + (b - a) * (b - a) * dfa / denom;
  
  return OTK::SUCCESS;
}

int quad_minimizer2 (double a, double dfa,
                     double b, double dfb,
                     double *alpha_min)
{
  if(std::isnan(a) || std::isinf(a) || std::isnan(dfa) || std::isinf(dfa) || 
     std::isnan(b) || std::isinf(b) || std::isnan(dfb) || std::isinf(dfb))
    return OTK::DOMAIN_ERROR;
  
  double denom = dfa - dfb;
  
  if(fabs(denom) < OTK::EPS)
    return OTK::ZERO_DIVISION;
  
  *alpha_min = b + (b - a) * dfb / denom;
  
  return OTK::SUCCESS;
}
