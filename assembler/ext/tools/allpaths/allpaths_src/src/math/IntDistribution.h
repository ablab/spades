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
//


#ifndef _MATH__INT_DISTRIBUTION_H
#define _MATH__INT_DISTRIBUTION_H

#include "MainTools.h"

#include "math/IntFunction.h"

class IntDistribution 
{
public:
  typedef double real_t;

private: 
  IntFunction<real_t> _prob;
  IntFunction<real_t> _prob_le;
  IntFunction<real_t> _prob_le_sum;

  void _normalize();
public:
  
  IntDistribution()
    : _prob(),
      _prob_le(),
      _prob_le_sum()
  {}
  
  IntDistribution(const IntFunction<real_t> & func) 
    : _prob(func), 
      _prob_le(func.x_min(), func.x_max(), 0.0),
      _prob_le_sum(func.x_min(), func.x_max(), 0.0)
  {
    _normalize();
  }


  ~IntDistribution()
  {
    //cout << "~int_dist(): " << _prob.x_min() << ", " << _prob.x_max() 
    //     << " le: " << _prob_le.x_min() << ", " << _prob_le.x_max() 
  }

  operator bool() const { return (_prob.size() > 1); }
  
  static IntDistribution gaussian(const int mean, 
                                  const int sigma,
                                  const real_t n_sigma = 4.0);
  static IntDistribution uniform(const int a,
                                 const int b);
  
  template <class INT_t>
  void from_hits(const vec<INT_t> & hits)
  {
    IntFunction<real_t> f;
    for (size_t i = 0; i != hits.size(); i++)
      f[hits[i]]++;
    _prob = f;
    _normalize();
  }


  int x_min() const { return _prob.x_min(); }
  int x_max() const { return _prob.x_max(); }


  real_t prob(const int x) const 
  { 
    if (x < x_min() || x > x_max()) return 0.0;
    return _prob[x]; 
  }
  
  real_t prob_le(const int x) const 
  {
    if (x < x_min()) return 0.0;
    if (x > x_max()) return 1.0;
    return _prob_le[x]; 
  }

  real_t prob_lt     (const int x)                const { return prob_le(x - 1); }
  real_t prob_ge     (const int x)                const { return 1.0 - prob_lt(x); }
  real_t prob_gt     (const int x)                const { return 1.0 - prob_le(x); }
  real_t prob_in     (const int x0, const int x1) const { return prob_le(x1) - prob_lt(x0); }
  real_t prob_not_in (const int x0, const int x1) const { return 1.0 - prob_in(x0, x1); }
  real_t prob_max    ()                           const { return prob(x_prob_max()); }
  
  // x_prob_le(random01) is a simple way of getting a random sample from the distribution. 
  int    x_prob_le   (const real_t F)             const { return quantile(F); } 
  int    x_prob_max  ()                           const { return _prob.x_f_max(); } 

  real_t mean        ()                                     const;
  real_t moment      (const real_t c, const unsigned order) const;
  real_t variance    ()                                     const { return moment(mean(), 2); }
  int    quantile    (const real_t Fq)                      const;
  int    median      ()                                     const { return quantile(0.5); }
  int    mode        ()                                     const { return x_prob_max(); }

  real_t prob_in_sum(const int x0, const int x1, const int n) const;
  


  IntFunction<real_t> probs() const { return _prob; }

  void split(IntDistribution * p_dist_left,
             IntDistribution * p_dist_right,
             const int x_split = 0) const;

  IntDistribution reverse() const;
  

  void to_text_file(const String & head) const;


public:
  // ---- SELF_SERIALIZABLE method
  size_t writeBinary(BinaryWriter& writer) const { return _prob.writeBinary(writer); }

  // ---- SELF_SERIALIZABLE method
  void readBinary(BinaryReader& reader) { _prob.readBinary(reader); _normalize(); }

  // ---- SELF_SERIALIZABLE method
  static size_t externalSizeof() { return 0; }

};



SELF_SERIALIZABLE(IntDistribution);









// Compute the distribution of the sum of 2 random variables
//   Z = X + Y    =>   p(Z) = p(X) * p(Y)
IntDistribution distribution_of_sum(const IntDistribution & x_dist,
                                    const IntDistribution & y_dist);

// Compute the distribution of the sum of 2 random variables
//   Z = X - Y    =>   p(Z) = p(X) * p(Y)
IntDistribution distribution_of_difference(const IntDistribution & x_dist,
                                           const IntDistribution & y_dist);










class IntLogDistribution : public IntFunction<double>
{
public:
  typedef double real_t;

public:
  IntLogDistribution(const int x0 = 0, const int x1 = 0, const real_t val = 0) 
    : IntFunction<real_t>(x0, x1, val) 
  {
    this->expand_infinity();
  }

  IntLogDistribution(const IntDistribution & dist)
    : IntFunction<real_t>(dist.x_min(), dist.x_max())
  {
    const int x0 = dist.x_min();
    const int x1 = dist.x_max();
    
    real_t p_min = -1; // undefined
    for (int x = x0; x <= x1; x++) {
      const real_t p = dist.prob(x);
      ForceAssertGe(p, 0);
      if (p > 0.0 && (p_min < 0 || p < p_min))
        p_min = p;
    }

    IntLogDistribution & me = *this;
    for (int x = x0; x <= x1; x++) {
      const real_t p = dist.prob(x);
      me[x] = (p > 0.0) ? log(p) : log(p_min);
    }        
    me.expand_infinity();
  }


  IntLogDistribution(const IntFunction<real_t> & func)
    : IntFunction<real_t>(func.x_min(), func.x_max())
  {
    const int x0 = func.x_min();
    const int x1 = func.x_max();
    
    real_t f_min = -1; // undefined
    for (int x = x0; x <= x1; x++) {
      const real_t f = func[x];
      ForceAssertGe(f, 0);
      if (f > 0.0 && (f_min < 0 || f < f_min))
        f_min = f;
    }

    IntLogDistribution & me = *this;
    for (int x = x0; x <= x1; x++) {
      const real_t f = func[x];
      me[x] = (f > 0.0) ? log(f) : log(f_min);
    }        
    me.expand_infinity();
  }


  ~IntLogDistribution()
  {
    // cout << "~int_log_dist(): " << x_min() << ", " << x_max() << endl;
  }


  IntDistribution int_distribution() const 
  {
    const int x0 = x_min();
    const int x1 = x_max();
    IntFunction<real_t> func(x0, x1);
    const IntLogDistribution & me = *this;

    real_t f_max = me.f_max(); 

    for (int x = x0; x <= x1; x++)
      func[x] = exp(me[x] - f_max);   // = 1 at me[x] = f_max

    // here I rely on the IntDistribution constructor that accepts an IntFunction
    return func;
  }


  void to_text_file(const String & fn) const
  {
    const real_t fmax = f_max();
    
    const int x0 = x_min();
    const int x1 = x_max();

    ofstream os;
    os.open(fn.c_str());
    os << "# x_min = " << x0 << endl;
    os << "# x_max = " << x1 << endl;
    os << "# 1:x  2:f_x" << endl;
    os << fixed;
    
    for (int x = x0; x <= x1; x++) {
      os << setw(10) << x << " "
	 << setw(16) << ((*this)[x] - fmax)
	 << endl;
    }
    os.close();
  }




};





#endif
