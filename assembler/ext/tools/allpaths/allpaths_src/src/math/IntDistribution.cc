///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
// author: Filipe Ribeiro      03/2011
// 
//



#include <math.h>

#include "math/IntDistribution.h"
#include "math/IntFrequencies.h"


typedef IntDistribution::real_t real_t;



// the reason to build the cummulative vecs is so that 
// in the end everything adds up to 1 rather than 1 + <epsilon> 
// because of round off.
void IntDistribution::_normalize()
{
  const int x0 = _prob.x_min();
  const int x1 = _prob.x_max();
    
  ForceAssertGe(_prob[x0], 0.0);
  _prob_le[x0] = _prob[x0];
  for (int x = x0 + 1; x <= x1; x++) {
    const real_t p = _prob[x];
    ForceAssertGe(p, 0.0);
    _prob_le[x] = _prob_le[x - 1] + p;
  }
  
  const real_t sum = _prob_le[x1];
  if (sum == 0.0) {
    _prob = IntFunction<real_t>(); // this signals an invalid distribution
  }
  else {
    const real_t norm = 1.0 / sum;
    
    if (norm != 1.0) { // admitedly unlikely
      for (int x = x0; x <= x1; x++) {
        _prob[x] *= norm;
        _prob_le[x] *= norm;
      }
    }
    
    _prob_le_sum[x0] = _prob_le[x0];
    for (int x = x0 + 1; x <= x1; x++)
      _prob_le_sum[x] = _prob_le_sum[x - 1] + _prob_le[x];
  }

}



int IntDistribution::quantile(const real_t Fq) const
{
  ForceAssertGe(Fq, 0.0);
  ForceAssertLe(Fq, 1.0);

  const int x0 = x_min();
  if (Fq == 0.0) return x0;
  const int x1 = x_max();
  if (Fq == 1.0) return x1;

  int x = x0;
  while (x <= x1 && _prob_le[x] < Fq) x++;
  ForceAssertLe(x, x1);
  return x;
}


real_t IntDistribution::mean() const
{
  const int x0 = x_min();
  const int x1 = x_max();
  real_t m = 0.0;
  for (int x = x0; x <= x1; x++) 
    m += x * _prob[x]; 
  return m;
}


real_t IntDistribution::moment(const real_t c, const unsigned order) const 
{
  const int x0 = x_min();
  const int x1 = x_max();
  real_t m = 0.0;
  for (int x = x0; x <= x1; x++) 
    m += pow(real_t(x) - c, static_cast<int>(order)) * _prob[x];
  return m;
} 



real_t IntDistribution::prob_in_sum(const int a, const int b, const int n) const
{
  ForceAssertLt(a, b);
  ForceAssertGt(n, 0);

  const int x0 = x_min();
  const int x1 = x_max();

  const int a0 = a - 2;
  const int an = a0 + n;

  const int b0 = b - 1;
  const int bn = b0 + n;

  // phi(x) = phi_r(x) + phi_i(x)

  const real_t phi_a0_r = (a0 < x0) ? 0 : (a0 >= x1) ? _prob_le_sum[x1] : _prob_le_sum[a0];
  const real_t phi_an_r = (an < x0) ? 0 : (an >= x1) ? _prob_le_sum[x1] : _prob_le_sum[an];
  const real_t phi_b0_r = (b0 < x0) ? 0 : (b0 >= x1) ? _prob_le_sum[x1] : _prob_le_sum[b0];
  const real_t phi_bn_r = (bn < x0) ? 0 : (bn >= x1) ? _prob_le_sum[x1] : _prob_le_sum[bn];

  const int phi_a0_i = (a0 <= x1) ? 0 : a0 - x1;
  const int phi_an_i = (an <= x1) ? 0 : an - x1;
  const int phi_b0_i = (b0 <= x1) ? 0 : b0 - x1;
  const int phi_bn_i = (bn <= x1) ? 0 : bn - x1;

  return (real_t(phi_bn_i - phi_b0_i - phi_an_i + phi_a0_i) +
	  (phi_bn_r - phi_b0_r - phi_an_r + phi_a0_r));
}






void IntDistribution::split(IntDistribution * p_dist_lt,
                            IntDistribution * p_dist_gt,
                            const int x_split) const
{
  IntFunction<real_t> f_lt = probs();
  IntFunction<real_t> f_gt = probs();

  real_t sum_lt = 0.0;
  real_t sum_gt = 0.0;
  for (int x = x_min(); x != x_max(); x++) {
    if (x < x_split) {
      f_gt[x] = 0;
      sum_lt += prob(x);
    }
    else if (x > x_split) {
      f_lt[x] = 0;
      sum_gt += prob(x);
    }
    else {  // (x == x_split)
      f_gt[x] *= 0.5;
      f_lt[x] *= 0.5;
      sum_gt += prob(x);
      sum_lt += prob(x);
    }
  }

  if (sum_lt > 0) *p_dist_lt = IntDistribution(f_lt);
  if (sum_gt > 0) *p_dist_gt = IntDistribution(f_gt);
}
 




IntDistribution IntDistribution::reverse() const 
{
  IntFunction<real_t> f(-x_max(), -x_min(), 0);

  for (int x = x_min(); x != x_max(); x++)
    f[-x] = prob(x);

  return f;
}










void IntDistribution::to_text_file(const String & head) const
{
  const int x0 = x_min();
  const int x1 = x_max();

  const String fn = head + ".prob";
  ofstream os;
  os.open(fn.c_str());
  os << "# x_min = " << x0 << endl;
  os << "# x_max = " << x1 << endl;
  os << "# 1:x  2:p(x) 3:F(x)" << endl;

  os << fixed;
  for (int x = x0; x <= x1; x++) {
    os << setw(10) << x << " "
       << setw(16) << setprecision(12) << _prob[x] << " "
       << setw(16) << setprecision(12) << _prob_le[x] << endl;
  }

  os.close();
}














IntDistribution IntDistribution::gaussian(const int mean, 
                                          const int sigma,
                                          const real_t n_sigma)
{
  const int dx = 0.5 + n_sigma * abs(sigma);
  
  IntFunction<double> func(mean - dx, mean + dx);
  
  const real_t inv_sigma = 1.0 / real_t(sigma);

  // no need to normalize since the normalization is done internally anyway
  //
  //const real_t norm = inv_sigma / sqrt(2 * 3.14159265358979);

  for (int x = -dx; x <= dx; x++) {
    const real_t z = x * inv_sigma; 
    func[mean + x] = exp(-0.5 * z * z);
    //func[mean + x] = norm * exp(-0.5 * z * z));
  }
  // here I rely on the IntDistribution constructor that accepts an IntFunction
  return func;
}



IntDistribution IntDistribution::uniform(const int a,
                                         const int b)
{
  ForceAssertLe(a, b);
  
  IntFunction<double> func(a, b);
  const real_t p = 1.0 / real_t(b - a + 1);
  for (int x = a; x <= b; x++)
    func[x] = p;
  
  // here I rely on the IntDistribution constructor that accepts an IntFunction
  return func;
}









// Compute the distribution of the sum of 2 random variables
//   Z = X + Y    =>   p(Z) = p(X) * p(Y)
//
//  I know, this is O(n^2).  When this becomes a bottleneck somebody will
//  reimplement this with FFT to make it O(n log n).  But I don't care now.
//
IntDistribution distribution_of_sum(const IntDistribution & x_dist,
                                    const IntDistribution & y_dist)
{
  const int x_min = x_dist.x_min();
  const int x_max = x_dist.x_max();

  const int y_min = y_dist.x_min();
  const int y_max = y_dist.x_max();

  const int z_min = x_min + y_min;
  const int z_max = x_max + y_max;
  
  IntFunction<IntDistribution::real_t> z_func(z_min, z_max);

  if (x_max - x_min > y_max - y_min) {      // more xs than ys -> cycle through ys
    for (int z = z_min; z <= z_max; z++) {
      real_t pz = 0.0;
      for (int y = y_min; y <= y_max; y++)
        pz += y_dist.prob(y) * x_dist.prob(z - y);
      z_func[z] = pz;
    }
  }
  else {                                    // more ys than xs -> cycle through xs
    for (int z = z_min; z <= z_max; z++) {
      real_t pz = 0.0;
      for (int x = x_min; x <= x_max; x++)
        pz += x_dist.prob(x) * y_dist.prob(z - x);
      z_func[z] = pz;
    }
  }
  // here I rely on the IntDistribution constructor that accepts an IntFunction
  return z_func;  
}



// Compute the distribution of the difference of 2 random variables
//   Z = X - Y    =>   p(Z) = p(X) * p(-Y)
//
//  I know, this is O(n^2).  See above.
//
IntDistribution distribution_of_difference(const IntDistribution & x_dist,
                                           const IntDistribution & y_dist)
{
  const int x_min = x_dist.x_min();
  const int x_max = x_dist.x_max();

  const int y_min = y_dist.x_min();
  const int y_max = y_dist.x_max();

  const int z_min = x_min - y_max;
  const int z_max = x_max - y_min;
  
  IntFunction<IntDistribution::real_t> z_func(z_min, z_max);

  if (x_max - x_min > y_max - y_min) {      // more xs than ys -> cycle through ys
    for (int z = z_min; z <= z_max; z++) {
      real_t pz = 0.0;
      for (int y = y_min; y <= y_max; y++)
        pz += y_dist.prob(y) * x_dist.prob(z + y);   // x = z + y
      z_func[z] = pz;
    }
  }
  else {                                    // more ys than xs -> cycle through xs
    for (int z = z_min; z <= z_max; z++) {
      real_t pz = 0.0;
      for (int x = x_min; x <= x_max; x++)
        pz += x_dist.prob(x) * y_dist.prob(x - z);   // y = x - z
      z_func[z] = pz;
    }
  }
  // here I rely on the IntDistribution constructor that accepts an IntFunction
  return z_func; 
}






