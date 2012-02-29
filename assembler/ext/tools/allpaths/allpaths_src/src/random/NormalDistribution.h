// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 



/*
 * Cambridge, March 27, 2001.
 *
 * NormalDistribution.cc
 *
 * Calculate normal distribution, normal deviate, etc.
 */
#ifndef NORMALDISTRIBUTION
#define NORMALDISTRIBUTION



/*
 * NormalDensity
 *
 * Normal random variable probability density function.
 */
float NormalDensity( float alpha, float mu, float sigma );



/*
 * StandardNormalDistribution
 *
 * Just the standard normal distribution.
 * Returns: the standard normal distribution by Chebyshev fitting.
 */
float StandardNormalDistribution( float alpha );



/*
 * NormalDistribution
 *
 * Normal distribution (with mean mu and standard deviation sigma.)
 * Returns: integral from -infinity to alpha of Normal(mu, sigma).
 */
float NormalDistribution( float alpha, float mu, float sigma );



/*
 * NormalDeviate
 *
 * Generate normally distributed random numbers (normal deviates), as per
 * algorithm R, p.125 in par.3.4.1 of "The Art of Computer Programming,
 * vol.2 (Seminumerical Algorithms)", by D. Knuth.
 *
 * Legenda:
 *  U: strictly positive real number;
 *  V: real number;
 *  X: result.
 *
 * Return:
 *  If U and V belong to a specific region in the (U, V) plane, then
 *  it returns true and fills X with a value depending on U, and V.
 *  Otherwise (i.e. if the point does not belong to the region)
 *  it returns false. The idea is that we can deduce a normally
 *  distributed random sequence from two given uniformly distributed 
 *  and independent random sequences U_n and V_n.
 * 
 * Remark 1: if X_m is normally distributed with mean 0 and standard 
 *  deviation 1, then (mu + sigma*X_m) is normally distributed with
 *  mean mu and standard deviation sigma.
 *
 * Remark 2: the region in the (U, V) plane is roughly contained in the
 *  rectangle (0, 1] x [-sqrt(2/e), +sqrt(2/e)]. More exactly, the region
 *  is given by:
 *  U in (0, 1];
 *  V in [-2u * sqrt(ln(1/U)), +2u * sqrt(ln(1/U))].
 */
bool NormalDeviate( double U, double V, double &X );

/// Given intervals of lengths len1 and len2, whose left ends
/// are separated (left2-left1) by distance sep_mean +- sep_sd, 
/// what is the probability that they overlap?
inline float IntervalOverlapProbability( int len1, int len2, 
					 int sep_mean, float sep_sd ) {
  return ( NormalDistribution( len1, sep_mean, sep_sd ) +
	   NormalDistribution( len2, -sep_mean, sep_sd ) - 1 );
}

#endif
