// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 


#include <math.h>

#include "system/Assert.h"
#include "random/NormalDistribution.h"


/*
 * NormalDensity
 */
float NormalDensity( float alpha, float mu, float sigma )
{
  ForceAssert( sigma > 0 );
  float const_part = 1.0 / ( sqrt( 2.0 * M_PI ) * sigma );
  float exp_part = - ( ( alpha-mu ) * ( alpha-mu ) ) / ( 2.0 * sigma * sigma );
  return const_part * exp( exp_part );
}


/*
 * StandardNormalDistribution
 */
float StandardNormalDistribution( float alpha )
{
  float t, z, ans;
  
  // Return 0 or 1 before -5 or after 5.
  if (alpha > 5)
    return 1;
  if (alpha < -5)
    return 0;
  
  z = fabs(alpha) / M_SQRT2;
  t = 1.0 / ( 1.0 + 0.5*z );
  ans = 0.5 * t *
    exp( -z*z - 1.26551223 + t*(1.00002368 +
				t*(0.37409196 + 
				   t*(0.09678418 + 
				      t*(-0.18628806 +
					 t*(0.27886807 +
					    t*(-1.13520398 +
					       t*(1.48851587 +
						  t*(-0.82215223 +
						     t*0.17087277)))))))));
  
  // Return.
  return alpha >= 0.0 ? 1.0 - ans : ans;
}


/*
 * NormalDistribution
 */
float NormalDistribution( float alpha, float mu, float sigma )
{
  ForceAssert ( sigma > 0 );
  return StandardNormalDistribution( ( alpha - mu ) / sigma );
}


/*
 * NormalDeviate
 */
bool NormalDeviate( double U, double V, double &X )
{
  ForceAssert ( U > 0 );
  double res = ( sqrt( 8.0 / M_E ) * ( V - 0.5 ) ) / U;
  if ( res * res <= -4.0 * log( U ) ) {
    X = res;
    return true;
  }
  return false;
}
