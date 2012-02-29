///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
#ifndef PAIR_DIST_MODELS_H
#define PAIR_DIST_MODELS_H

#include "CoreTools.h"
#include "MainTools.h"
#include "PairsManager.h"
#include "Superb.h"
#include "VecTemplate.h"
#include "VecUtilities.h"
#include "math/Functions.h"
#include "paths/Alignlet.h"
#include "paths/FixScaffoldsCore.h"
#include "paths/ScaffoldsUtils.h"
#include "random/Random.h"
#include "random/NormalRandom.h"
#include <map>
#include "math/IntDistribution.h"

// =====================================================
// Class ProbFuncInterface
//
// abstract class for pair distribution correction
// =====================================================
class ProbFuncInterface {
 public:
  virtual double p0(double t1) const = 0; // probability density function, p0(x)
  virtual double f0(double t1,double t2) const = 0; // int_t1^t2 { p0(x) }
  virtual double f1(double t1,double t2) const = 0; // int_t1^t2 { p0(x)*x }
  virtual double f2(double t1,double t2) const = 0; // int_t1^t2 { p0(x)*x*x }
  virtual double Avg() const = 0; // average of the distribution
  virtual double SD() const = 0; // standard deviation of the distribution
  virtual double L() const = 0; // p0(t) = 0 for t < L()
  virtual double U() const = 0; // p0(t) = 0 for t > U()
};

// ==============================================================================
// Class ProbFuncIntDist 
//
// pair distribution from IntDistribution
// ===============================================================================

// Note that in scaffold_reads, the distribution from IntDistribution are between
// the invariant end of the read pairs, different from the end-to-end distance 
// reported in pairsManager.
//
class ProbFuncIntDist : public ProbFuncInterface {
 private:
  double min, max; // minimal and maximum range in t-axis
  double avg, sd;
  double step;
  vec<double> p; // probability distribution
  vec<double> d0, d1, d2; // cumulative integrations

 public:
  ProbFuncIntDist(const IntDistribution& dist, const int shift = 0, bool TRIM = true,  
                  const int MIN_SEP = -1000, const int MAX_SEP = 1000000000 );
  ProbFuncIntDist( const vec<double>& array, const int MIN_SEP = 0 );
  ProbFuncIntDist(){};

  void FromArray( const vec<double>& array, const int MIN_SEP = 0 );

  double f0(double t) const { return seekArray(t, d0); }
  double f1(double t) const { return seekArray(t, d1); }
  double f2(double t) const { return seekArray(t, d2); }

  // the implementation of virtual methods
  double p0(double t) const { return seekArray(t, p); }
  double f0(double t1,double t2) const {return f0(t2) - f0(t1);}
  double f1(double t1,double t2) const {return f1(t2) - f1(t1);}
  double f2(double t1,double t2) const {return f2(t2) - f2(t1);}
  double Avg() const { return avg; }
  double SD() const  { return sd;  }
  double L() const  { return min;  }
  double U() const  { return max;  }
 
 private:
  // some private methods
  double seekArray(double t, const vec<double>& dd) const {
    int index = (t - min) / step;
    if ( index < 0 ) index = 0;
    if ( index >= dd.isize() ) index = dd.isize() -1;
    return dd[index];
  }
};

// Pair distribution from an ascii file 
class ProbFuncDist : public ProbFuncInterface {
  map<double, double> dist;
  map<double, double> dist0;
  map<double, double> dist1;
  map<double, double> dist2;
  double avg, sd;
 public:
  ProbFuncDist(String filename) ;
  double f0(double t) const { return seekArray(t, dist0); }
  double f1(double t) const { return seekArray(t, dist1); }
  double f2(double t) const { return seekArray(t, dist2); }

  // the implementation of virtual methods
  double p0(double t) const { return seekArray(t, dist); }
  double f0(double t1,double t2) const {return f0(t2) - f0(t1);}
  double f1(double t1,double t2) const {return f1(t2) - f1(t1);}
  double f2(double t1,double t2) const {return f2(t2) - f2(t1);}
  double Avg() const { return avg; }
  double SD()  const { return sd;  }
  double L() const  { return dist.begin()->first;  }
  double U() const  { return dist.rbegin()->first;  }

 private:
  // some private methods
  double seekArray(double t, const map<double,double>& dd) const
  {
    map<double,double>::const_iterator lb;
    lb = dd.lower_bound(t);
    if (lb == dist.end()) lb--;
    return lb->second;
  }
};

// Basic statistic functions for a normal distribution(avg, std)
// Note that the distribution has (avg_/sqrt(2), std_/sqrt(2))
// if the pairs separation has (avg_, std_)
//
// f0 = int(exp(-x^2), -inf, t) = (pi^(1/2)*(erf(t) + 1))/2
// f1 = int(exp(-x^2)*x, -inf, t)  = -1/(2*exp(t^2))
// f2 = int(exp(-x^2)*x*x, -inf, t) f = pi^(1/2)/2 - t/(2*exp(t^2)) - (pi^(1/2)*erfc(t))/4
// f3 = int(exp(-x^2)*x*x*x, -inf, t)  = -(t^2 + 1)/(2*exp(t^2))
class ProbFuncNormal: public ProbFuncInterface {
  double pi;
  double sqrt2;
  double avg, std; 
  double min, max;

  double f0_normal(double t) const {
    if ( t < -9 ) return 0.0 ;
    if ( t > 9 ) return sqrt(pi) ;
    return sqrt(pi)*(erf(t)+1.0)/2.0 ;
  }
  double f1_normal(double t) const {
    if ( t < -9 ) return 0;
    if ( t > 9 ) return 0;
    return  -0.5*exp(-t*t);
  }
  double f2_normal(double t) const {
    if ( t < -9 ) return -0.25*sqrt(pi);
    if ( t > 9 ) return 0.25*sqrt(pi);
    return 0.25*(sqrt(pi)*erf(t) - 2.0*t*exp(-t*t));
  }

 public:
  ProbFuncNormal(double avg_, double std_) {  
    pi = acos(-1.0); sqrt2 = sqrt(2.0);
    avg = avg_/sqrt2;
    std = std_/sqrt2;
    min = avg - 9 * std; // p0(t) = 0 if t < min or t > max
    max = avg + 9 * std;
  }
  double f0(double t) const {
    t = ( t - avg )/std;
    return f0_normal(t);
  }
  double f1(double t) const {
    t = (t - avg)/std;
    return f1_normal(t)*std + f0_normal(t)*avg;
  }
  double f2(double t) const {
    t = (t - avg)/std;
    return f2_normal(t)*std*std + f1_normal(t)*2*std*avg + f0_normal(t)*avg*avg;
  }
  // the implementation of virtual methods
  double p0(double t) const { 
    t = ( t - avg )/std;
    return(exp(-t*t)); 
  }
  double f0(double t1,double t2) const {return f0(t2) - f0(t1);}
  double f1(double t1,double t2) const {return f1(t2) - f1(t1);}
  double f2(double t1,double t2) const {return f2(t2) - f2(t1);}
  double Avg() const { return avg; }
  double SD()  const { return std;  }
  double L() const  { return min;  }
  double U() const  { return max;  }
};

template <class T> 
double Sum0( const T& t, double x1, double x2, double x3, double x4) {
  double S0 = t.f1(x1, x2) - x1 * t.f0(x1, x2) + x4 * t.f0(x3, x4) - t.f1(x3, x4) + (x2 - x1) * t.f0(x2, x3); 
  return 2 * S0;
}

template <class T>
double Sum0( const T& t, double x1, double x2, double x3, double x4, double xm ) {
  double h = x2 - x1;
  double S0;
  if ( xm < x1 ) {
    S0 = 0;
  } else if ( xm < x2 ){
    S0 = t.f1(x1, xm) - x1 * t.f0(x1, xm);
  } else if ( xm < x3 ) {
    S0 = t.f1(x1, x2) - x1 * t.f0(x1, x2) +  h * t.f0(x2, xm); 
  } else if ( xm < x4 ) {
    S0 = t.f1(x1, x2) - x1 * t.f0(x1, x2) + x4 * t.f0(x3, xm) - t.f1(x3, xm) + h * t.f0(x2, x3); 
  } else {
    S0 = t.f1(x1, x2) - x1 * t.f0(x1, x2) + x4 * t.f0(x3, x4) - t.f1(x3, x4) + h * t.f0(x2, x3); 
  }
  return 2 * S0;
}

template <class T> 
double Sum1(const T& t, double x1, double x2, double x3, double x4) {
  double S1 = t.f2(x1, x2) - x1 * t.f1(x1, x2) + x4 * t.f1(x3, x4) - t.f2(x3, x4) + (x2 - x1) * t.f1(x2, x3);
  return 2 * S1;
}


template <class T> 
double Sum1(const T& t, double x1, double x2, double x3, double x4, double xm) {
  double h = x2 - x1;
  double S1;
  if ( xm < x1 ) {
    S1 = 0;
  } else if ( xm < x2 ){
    S1 = t.f2(x1, xm) - x1 * t.f1(x1, xm);
  } else if ( xm < x3 ) {
    S1 = t.f2(x1, x2) - x1 * t.f1(x1, x2) + h * t.f1(x2, xm); 
  } else if ( xm < x4 ) {
    S1 = t.f2(x1, x2) - x1 * t.f1(x1, x2) + h * t.f1(x2, x3) + x4 * t.f1(x3, xm) - t.f2(x3, xm) ; 
  } else {
    S1 = t.f2(x1, x2) - x1 * t.f1(x1, x2) + h * t.f1(x2, x3) + x4 * t.f1(x3, x4) - t.f2(x3, x4);
  }
  return 2 * S1;
}

// Observed average pair separation over the region
template <class T>
double MeanDistObserve(const T& t, double x1, double x2, double x3, double x4) {
  double S0 = Sum0( t, x1, x2, x3, x4 );
  double S1 = Sum1( t, x1, x2, x3, x4 );
  if ( S0 > 0 )
    S1 = S1 / S0; // t.first moment
  else {
    // in undefined region. But which side ?
    double s1 = abs( x4 - t.L() );
    double s4 = abs( x1 - t.U() );
    if ( s1 < s4 ) S1 = x4;
    else  S1 = x1;
  }
  return S1;
}

// For a given pair separation distribution <T>
// how much is the observed mean pair separations of the gap?
template <class T>
int MeanSepObserve(T t, int gap0, int len1,int len2);

// Total probability observed over the region
template <class T>
double SumProbObserve(T t, int gap0, int len1,int len2);

// Total probability observed over the region up to a cutoff of sumX1X2
template <class T>
double SumProbObserve(T t, int gap0, int len1, int len2, int sumX1X2Max);

template <class T>
void MostProbableGap( const T& distr, const int len1, const int len2, const vec< pair< int, int > >& links, 
    int& gap, int& std, bool verbose=false );


// ================================================================================= //
// =                  implementations of template functions                        = //
// ================================================================================= //


// For a given pair separation distribution <T>
// how much is the observed mean pair separations of the gap?
template <class T>
int MeanSepObserve(T t, int gap0, int len1,int len2){
  const double scaling = 1.0 / sqrt(2.0);
  double l1p = len1 * scaling; // length projected on the axis
  double l2p = len2 * scaling;
  // the corrdinates transformation
  // x = ( x0 + -lib_avg + gap0 ) * scaling
  //double x1 =  ( - lib_avg + gap0 ) * scaling;
  double x1 = gap0 * scaling;
  double h = min(l1p, l2p);
  double x2 = x1 + h;
  double x4 = x1 + l1p + l2p;
  double x3 = x4 - h;
  double actualSepMean = MeanDistObserve(t, x1, x2, x3, x4)/scaling;
  //cout << "x1x2x3x4= "<<  x1 <<" "<<x2<<" "<<x3<<" "<<x4<<endl;
  //cout << "actualSepMean= " << actualSepMean << endl;
  // the observed shift
  return actualSepMean;
}

// Total probability observed over the region
template <class T>
double SumProbObserve(T t, int gap0, int len1,int len2)
{
  // calculate the coordinates
  const double scaling = 1.0 / sqrt(2.0);
  double l1p = len1 * scaling; // length projected on the axis
  double l2p = len2 * scaling;
  double x1 = gap0 * scaling;
  double h = min(l1p, l2p);
  double x2 = x1 + h;
  double x4 = x1 + l1p + l2p;
  double x3 = x4 - h;
  return Sum0( t, x1, x2, x3, x4 );
}

// Total probability observed over the region, with restriction  x1 + x2 < sumX1X2
template <class T>
double SumProbObserve(T t, int gap0, int len1, int len2, int sumX1X2Max)
{
  // calculate the coordinates
  const double scaling = 1.0 / sqrt(2.0);
  double l1p = len1 * scaling; // length projected on the axis
  double l2p = len2 * scaling;
  double x1 = gap0 * scaling;
  double h = min(l1p, l2p);
  double x2 = x1 + h;
  double x4 = x1 + l1p + l2p;
  double x3 = x4 - h;
  double xm = x1 + sumX1X2Max * scaling;
  return Sum0( t, x1, x2, x3, x4, xm );
}

#endif
