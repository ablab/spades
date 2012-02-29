///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "STLExtensions.h"
#include "math/Functions.h"




/*  CombineNormalDistribution              Filipe Ribeiro 2009-06-22
 *
 *   When you have different measurements (normaly distributed) of the same quantity, 
 *   you can compute the resultant normal distribution.   
 * 
 *   The resultant normal distribution (mu, sigma) satisfies:
 * 
 *      1 / sigma^2  =  sum(     1 / sigma[i]^2 )
 * 
 *     mu / sigma^2  =  sum( mu[i] / sigma[i]^2 )
 *
 */

normal_distribution CombineNormalDistributions(const vec<normal_distribution> & ds,
                                               float * score)
{
  size_t n = ds.size();

  // do everything with 'double' and convert to 'float' only at the end 
  double one_sig2 = 0; //  1   / sigma^2
  double mu_sig2  = 0; //  mu  / sigma^2 
  
  for (size_t i = 0; i != n; i++) {
    const double x_i   = ds[i].mu_; 
    const double sig_i = ds[i].sigma_;

    const double one_sig2_i = 1.0 / (sig_i * sig_i);

    one_sig2 +=        one_sig2_i;
    mu_sig2  +=  x_i * one_sig2_i;
  }
  
  const double sig2 = 1.0 / one_sig2;


  // if a pointer for the score is provided calculate the score
  // score is a ratio between the observed spread and the expected spread
  // so it should be close to 1.
  // If score is much larger than 1 (say 3), then you probably don't have a good gaussian. 
  if (score) {
    const double mu2_nsig2 = mu_sig2 * mu_sig2 * sig2 / static_cast<double>(n);
      
    double s2_sig2 = 0;
    for (size_t i = 0; i != n; i++) {
      const double x_i   = ds[i].mu_; 
      const double sig_i = ds[i].sigma_;
      
      // add x2 and subtract mu2/n to keep round off error low.
      s2_sig2 += (x_i * x_i) / (sig_i * sig_i) - mu2_nsig2;
    }
    
    // s2_sig2 should always be >= 0, but you never know...
    *score = (s2_sig2 <= 0.0) ? 0.0 : sqrt(s2_sig2 / static_cast<double>(n));
  }

  return normal_distribution(mu_sig2 * sig2, sqrt(sig2));
}






/* ClusterNormalDistribution              Filipe Ribeiro 2010-04-16
 *  
 *  Extract clusters based on overlaps of n_sigma x sigma.
 */
vec< vec<normal_distribution> >
ClusterNormalDistributions(const vec<normal_distribution> & ds,
                           const float n_sigma)
{
  const size_t n = ds.size();
  bool sorted = true;
  for (size_t i = 1; i != n && sorted; i++) 
    sorted = (ds[i-1].mu_ <= ds[i].mu_);

  vec<normal_distribution> dss;
  if (!sorted) {
    dss = ds;
    sort(dss.begin(), dss.end());
  }
  
  const vec<normal_distribution> & myds = (sorted) ? ds : dss;

  vec< vec<normal_distribution> > clusters(1, vec<normal_distribution>(1, myds[0]));
  size_t i_cluster = 0;

  for (size_t i = 1; i != n; i++) {
    
    const bool overlap = (myds[i-1].mu_ + n_sigma * myds[i-1].sigma_ > 
                          myds[i].mu_   - n_sigma * myds[i].sigma_); 

    if (overlap) {
      clusters[i_cluster].push_back(myds[i]);
    }
    else {
      clusters.push_back(vec<normal_distribution>(1, myds[i]));
      i_cluster++;
    }
  }
  
  return clusters;
}






/*
 * SafeMeanStdev
 *
 * Calculate Mean and Standard deviation of a sef ot points, in a roundoff
 * error free way (see Knuth, The Art of Computer Programming, vol. 2.)
 *
 * Legenda:
 *  points: input data.
 *
 * Return:
 *  mu and sigma for the given data. If points has only one entry, then
 *  it returns (mu=points[0]; sigma=0). If any error occurs, it returns
 *  (mu=-1; sigma=-1).  NaN and inf values in the input data are silently
 *  skipped.
 */
normal_distribution SafeMeanStdev(const vec<float> & points)
{
  unsigned int num_points = points.size();
  
  normal_distribution result;
  result.mu_ = -1;
  result.sigma_ = -1;

  if ( 0 == num_points )
    return result;
  
  //find the first valid i and initialize the recurrence
  unsigned int i = 0;
  while (i<num_points && !isfinite(points[i]))
    ++i;
  if (i==num_points)
    return result;
  double mu = points[i];
  double sigma = 0;
  unsigned int validPoints = 1;

  //for remaining valid points, use Knuth recurrence.
  for ( unsigned int j = i+1; j < num_points; ++j ) {
    const double point = points[j];
    if (isfinite(point)) {
      ++validPoints;
      const double last_mu = mu;
      mu += (point - last_mu) / validPoints;
      sigma += ( point - last_mu ) * ( point - mu );
    }
  }
  
  if ( sigma < 0 )
    return result;
  
  double n_points = validPoints;
  if ( sigma > 0 ) {
    result.mu_ = mu;
    result.sigma_ = sqrt( sigma / ( n_points - 1 ) );

    return result;
  }

  // Last case, if sigma=0;
  result.mu_ = mu;
  result.sigma_ = 0;

  return result;
}

normal_distribution SafeMeanStdev(int count, double sum, double sumsq) {
  normal_distribution ret(-1,-1);
  if (!count) return ret;
  ret.mu_ = sum/count;
  if (1 == count) {
    ret.sigma_ = 0;
  } else {
    ret.sigma_ = sqrt((sumsq - (sum*sum/count))/(count - 1));
  }
  return ret;
}

//For the general case, we solve the following equation:
//Minimize the sum of the areas under the curves given the cutoff x
//We multiply the area of d1 by the ratio, and end up with this equation:
//log(ratio) - ((x-d1.mu_)/d1.sigma_)2 = - ((x-d2.mu_)/d2.sigma_)2
//where the 2 means squared.
//if ratio is 1 then we can get rid of the squares and have the simple case.
//otherwise we end up with the quadratic equation in the general case.  
double OptimalCutoff(const normal_distribution & d1,
		     const normal_distribution & d2,
		     const double ratio) {
  Assert(d1.sigma_ >= 0);
  Assert(d2.sigma_ >= 0);
  Assert(d1.mu_ < d2.mu_);

  //return weighted mean of means if either sigma is 0.
  if (0 == d1.sigma_ || 0 == d2.sigma_) { 
    return (ratio * d1.mu_ + d2.mu_) / (1+ratio);
  }

  if (1 == ratio) {//simplest case
    return ((d1.sigma_ * d2.mu_ + d1.mu_ * d2.sigma_) 
	    / (d1.sigma_ + d2.sigma_));
  } 
  else { //quadratic parameters
    double s = (d1.sigma_*d1.sigma_) / (d2.sigma_*d2.sigma_);
    double a = 1-s;
    double b = -2*(d1.mu_ - s*d2.mu_);
    double c = d1.mu_*d1.mu_- s*d2.mu_*d2.mu_ 
      - (d1.sigma_*d1.sigma_)*log(ratio);
    QuadraticFunction quad(a,b,c);
    pair<double,double> sol = quad.solutions();
    //if (!isfinite(sol.first) || sol.first == sol.second) return sol.first;
    
    //Now the general case: which of the two solutions is between the means?
    //PRINT4(sol.first, sol.second, d1.mu_, d2.mu_);
    if (sol.first >= d1.mu_ && sol.first <= d2.mu_) return sol.first;
    if (sol.second >= d1.mu_ && sol.second <= d2.mu_) return sol.second;
    return numeric_limits<double>::quiet_NaN();
  }
}

pair<double,double> QuadraticFunction::solutions() const {
  if (0 == a) { //not quadratic!
    double sol = numeric_limits<double>::quiet_NaN();
    if (0 == b) {
      if (0 == c) sol=0;
    } else {
      sol = -c/b;
    }
    return make_pair(sol,sol);
  } 
      
  double discrim = b*b -4*a*c;
  if (discrim < 0) return make_pair(numeric_limits<double>::quiet_NaN(),
				    numeric_limits<double>::quiet_NaN());
  else return make_pair((-b - sqrt(discrim))/(2*a),
			(-b + sqrt(discrim))/(2*a));
}

double Poisson(double lambda, unsigned int n) {
  double ret = exp(-lambda);
  for (unsigned int i = 1; i <= n; ++i) {
    ret *= lambda/(double)i;
  }
  return ret;
}

double PoissonCdf(double lambda, unsigned int n) {
  double term = exp(-lambda);
  double cum = term;
  for (unsigned int i = 1; i <= n; ++i) {
    term *= lambda/(double)i;
    cum += term;
  }
  return cum;
}

int InversePoissonCdf(double lambda, double p) {
  double term = exp(-lambda);
  double cum = term;
  int numEvents = 0;
  while ( ! ( cum > p ) ) {
    ++numEvents;
    term *= lambda/(double)numEvents;
    cum += term;
  }
  return numEvents;
}




template<class T> T N50( const vec<T>& v )
{    
  ForceAssert( v.size( ) > 0 );

  longlong sum = 0, half = 0;
  for ( int i = 0; i < (int)v.size( ); i++ ) {
    ForceAssert( v[i] > 0 );
    sum += v[i];
  }

  // If v is reverse-sorted.
  if ( is_sorted( v.rbegin(), v.rend() ) ) {
    for ( int i = v.size()-1; i >= 0; i-- ) {
      half += v[i];
      if ( 2 * half == sum && i > 0 ) 
        return (v[i] + v[i-1])/2;
      if ( 2 * half >= sum )
	return v[i];
    }

    return 0; // never executed
  }
  
  // If v is not sorted, copy it on a local vector and sort it.
  vec<T> sorted_v;
  bool use_given_v = true;
  if ( !is_sorted( v.begin(), v.end() ) ) {
    sorted_v = v;
    use_given_v = false;
    sort( sorted_v.begin(), sorted_v.end() );
  }
  
  const vec<T> &the_vector = ( use_given_v ) ? v : sorted_v;

  for ( int i = 0; i < (int) the_vector.size( ); i++ ) {
    half += the_vector[i];
    if ( 2 * half == sum && i < (int) the_vector.size( ) - 1 ) 
      return (the_vector[i] + the_vector[i+1])/2;
    if ( 2 * half >= sum ) return the_vector[i];
  }
  
  return 0; // never executed
}

template int N50( const vec<int>& v );
template unsigned int N50( const vec<unsigned int>& v );
template longlong N50( const vec<longlong>& v );
template unsigned long N50( const vec<unsigned long>& v );

// N50_size: Returns the N50 of the vector sizes in the input vec<vec<T> >.
template<class T> int N50_size( const vec<vec<T> >& v )
{
  int N = v.isize( );
  vec<int> sizes( N, 0 );
  for ( int i = 0; i < N; i++ )
    sizes[i] = v[i].isize( );
  
  Sort( sizes );
  return N50( sizes );
}

template int N50_size( const vec<vec<int> >& v );

template<class T> vec<T> NStatistics( const vec<T>& v )
{
  ForceAssert( v.size( ) > 0 );

  vec<T> stats( 9, 0 );

  longlong sum = 0;
  for ( unsigned int i = 0; i < v.size( ); ++i ) {
    ForceAssert( v[i] > 0 );
    sum += v[i];
  }

  vec<longlong> targets( stats.size() );
  for ( unsigned int i = 0; i < stats.size(); ++i )
    targets[i] = sum * (i+1)*10 / 100;

  longlong partial_sum = 0;
  unsigned int target_idx = 0;

  if ( is_sorted( v.begin(), v.end() ) )
  {
    for ( typename vec<T>::const_reverse_iterator v_iter = v.rbegin(); 
          v_iter != v.rend() && target_idx < targets.size(); ++v_iter )
    {
      partial_sum += *v_iter;
      
      while ( partial_sum >= targets[target_idx] )
      {
	stats[target_idx] = *v_iter;
	
	if ( ++target_idx >= targets.size() )
	  break;
      }
    }
  }

  else if ( is_sorted( v.rbegin(), v.rend() ) ) 
  {
    for ( typename vec<T>::const_iterator v_iter = v.begin(); 
          v_iter != v.end() && target_idx < targets.size(); ++v_iter )
    {
      partial_sum += *v_iter;
      
      while ( partial_sum >= targets[target_idx] )
      {
	stats[target_idx] = *v_iter;
	
	if ( ++target_idx >= targets.size() )
	  break;
      }
    }
  }
  
  else
  {
    cerr << "Vector not sorted." << endl;
    ForceAssert( 0 );
  }

  return stats;
}

template vec<int> NStatistics( const vec<int>& v );

template vec<longlong> NStatistics( const vec<longlong>& v );

template vec<size_t> NStatistics( const vec<size_t>& v );

Bool CombineMeasurements( const double g1, const double g2, const double d1,
                          const double d2, const double dmult, double& g, double& d )
{    if ( d1 == 0 && d2 == 0 )
     {    if ( g1 != g2 ) return False;
          g = g1;
          d = d1;    }
     else if ( d1 == 0 )
     {    g = g1;
          d = d1;
          if ( Abs( g2 - g ) > dmult * d2 ) return False;    }
     else if ( d2 == 0 )
     {    g = g2;
          d = d2;
          if ( Abs( g1 - g ) > dmult * d1 ) return False;    }
     else
     {    double w1 = 1.0/(d1*d1), w2 = 1.0/(d2*d2);
          g = ( w1*g1 + w2*g2 ) / (w1+w2);
          // This statement (i.e. the old code)
          //   d = sqrt( w1*w1*d1*d1 + w2*w2*d2*d2 ) / (w1+w2);
          // is equivalent to
          //   d = sqrt( w1*w1*(1/w1) + w2*w2*(1/w2) ) / (w1+w2);
          // which is equivalent to
          //   d = sqrt( w1+w2 ) / (w1+w2);
          // which is equivalent to
          d = 1.0 / sqrt( w1+w2 );
          if ( Abs( g1 - g ) > dmult * d1 ) return False;
          if ( Abs( g2 - g ) > dmult * d2 ) return False;    }
     return True;    }

double InverseNormalCDF(double p, double mu, double sigma) 
{
  // Algorithm by Peter J. Acklam from http://home.online.no/~pjacklam/index.html
  // Documentation claims error has absolute value less than 1.15 * 10^-9 in the entire region.
  const double a1 = -39.69683028665376;
  const double a2 = 220.9460984245205;
  const double a3 = -275.9285104469687;
  const double a4 = 138.3577518672690;
  const double a5 =-30.66479806614716;
  const double a6 = 2.506628277459239;

  const double b1 = -54.47609879822406;
  const double b2 = 161.5858368580409;
  const double b3 = -155.6989798598866;
  const double b4 = 66.80131188771972;
  const double b5 = -13.28068155288572;

  const double c1 = -0.007784894002430293;
  const double c2 = -0.3223964580411365;
  const double c3 = -2.400758277161838;
  const double c4 = -2.549732539343734;
  const double c5 = 4.374664141464968;
  const double c6 = 2.938163982698783;

  const double d1 = 0.007784695709041462;
  const double d2 = 0.3224671290700398;
  const double d3 = 2.445134137142996;
  const double d4 = 3.754408661907416;

  const double p_low =  0.02425;
  const double p_high = 1 - p_low;

  double q, x=0, r;

  ForceAssert(0.0<p);
  ForceAssert(p<1.0);
  if (0 < p && p < p_low) { //Rational approximation for lower region.
    q = sqrt(-2*log(p));
    x = (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / ((((d1*q+d2)*q+d3)*q+d4)*q+1);
  } else if (p_low <= p && p <= p_high) { //Rational approximation for central region.
    q = p - 0.5;
    r = q*q;
    x = (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q / (((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1);
  } else if (p_high < p && p < 1) { //Rational approximation for upper region.
    q = sqrt(-2*log(1-p));
    x = -(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / ((((d1*q+d2)*q+d3)*q+d4)*q+1);
  }
  return mu+(x*sigma);
}

double EstProbability(int passes, int trials)
{
  // We use the Wilson point estimate; see reference below.  Note it's
  // not necessarily the center of the confidence interval described
  // below, but close to it.
  return double(passes+2)/double(trials+4);
}


pair<double, double> BinomialConfidenceInterval(int passes, int trials, double error)
{
  // We use the Wilson score interval, which is the interval of p
  // values for which the z-score of the observation is in the desired
  // range.  It performs well even for very small #trials and is
  // tighter than the "exact" interval.  See Agresti, A., and Coull,
  // B. "Approximate is better than 'exact' for interval estimation of
  // binomial proportions." The American Statistician 52: 119-126,
  // 1998.
  const double z = InverseNormalCDF(1-error/2.0);
  const double oneOverN = 1.0/double(trials);
  const double z2OverN = z*z*oneOverN;
  const double p = double(passes) * oneOverN; // really p-hat==observed p
  const double center = (p + 0.5 * z2OverN) / (1+ z2OverN);
  const double halfwidth = (z*sqrt(oneOverN * (p*(1.0-p)+0.25*z2OverN)))/(1+z2OverN);
  return make_pair(center-halfwidth, center+halfwidth);
}

float QualScoreFloat(int correct, int wrong) {
  double errprob = (1.0 + wrong)/(1.0 + wrong + correct);
  Assert(errprob <= 1);
  Assert(errprob > 0);
  return -10*log(errprob)/log(10);
}
