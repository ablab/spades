
#include "OTK.h"
#include "PolyFit.h"

#include <algorithm>
#include <cassert>
#include <math.h>

int cubic_minimizer (double a, double fa, double dfa,
                     double b, double fb, double dfb,
                     int extrapolate, double *alpha_min)
{
  double z, w;
  double denom;
  double D;
  
  if(isnan(a) || isinf(a) || isnan(fa) || isinf(fa) || isnan(dfa) || isinf(dfa) || 
     isnan(b) || isinf(b) || isnan(fb) || isinf(fb) || isnan(dfb) || isinf(dfb))
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
  /*if(isnan(*alpha_min) || isinf(*alpha_min))
     return GSL_FAILURE;*/
  assert(!isnan(*alpha_min) && !isinf(*alpha_min));
  
  return OTK::SUCCESS;
}

int quad_minimizer1 (double a, double fa, double dfa,
                     double b, double fb,
                     double *alpha_min)
{
  if(isnan(a) || isinf(a) || isnan(fa) || isinf(fa) || isnan(dfa) || isinf(dfa) || 
     isnan(b) || isinf(b) || isnan(fb) || isinf(fb))
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
  if(isnan(a) || isinf(a) || isnan(dfa) || isinf(dfa) || 
     isnan(b) || isinf(b) || isnan(dfb) || isinf(dfb))
    return OTK::DOMAIN_ERROR;
  
  double denom = dfa - dfb;
  
  if(fabs(denom) < OTK::EPS)
    return OTK::ZERO_DIVISION;
  
  *alpha_min = b + (b - a) * dfb / denom;
  
  return OTK::SUCCESS;
}
