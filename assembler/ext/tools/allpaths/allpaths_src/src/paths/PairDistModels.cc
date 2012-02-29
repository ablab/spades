///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
#include "paths/PairDistModels.h"

// Initialize the probability functon from the IntDistribution. Allow optional shift
// to generate an end-to-end pair separation, also allow distribution to be trimmed
// to contain only single peak, and trimmed for given MIN_SEP and MAX_SEP.
ProbFuncIntDist:: ProbFuncIntDist(const IntDistribution& dist, const int shift, bool TRIM,
                               const int MIN_SEP, const int MAX_SEP) 
{
  bool verbose = false;
  // avg and std of the original distribution
  double s_avg = dist.mean();
  double s_sd = sqrt( dist.moment( s_avg, 2 ) );
  avg = s_avg / sqrt(2.0);
  sd = s_sd / sqrt(2.0);
  //cout << "original avg std " <<  avg << " " << sd << endl;
  ///int MIN_SEP_new = dist.x_min();
  int MIN_SEP_new = Max( dist.x_min(), 0 ); // force MIN_SEP_new > 0;
  int MAX_SEP_new = dist.x_max();
  if ( TRIM)
  {
    // Distribution trimming 
    // Select only one peak from the distribution. Also trim 0.1% ditribution from both side 
    double prob_cutoff = 0.001;
    int maxProbArg = 0; 
    double maxProb = dist.prob(maxProbArg); // no innies
    //int maxProbArg = dist.x_min(), maxProb = dist.prob(maxProbArg);
    for ( int i = 0; i < dist.x_max(); i++ ) {
      if ( dist.prob(i) > maxProb ) {
	maxProb = dist.prob(i);
	maxProbArg = i;
      }
    }
    int iU = maxProbArg + s_sd, iL =  maxProbArg - s_sd; // upperbound and lowbound
    for ( ; iL > dist.x_min(); iL-- )
      if ( dist.prob(iL+1) - dist.prob(iL) < 0 || dist.prob_lt(iL) < prob_cutoff ) break;
    for ( ; iU < dist.x_max(); iU++ )
      if ( dist.prob(iU) - dist.prob(iU+1) < 0 || dist.prob_gt(iU) < prob_cutoff ) break;
    if ( verbose ) cout << "  Distribution trimming min_x= " << iL << " max_x= " << iU << endl;
    MIN_SEP_new = iL;
    MAX_SEP_new = iU;
  }
  // min and max in invariant pair distance 
  int minSep = Max ( MIN_SEP - shift, MIN_SEP_new );
  int maxSep = Min ( MAX_SEP - shift, MAX_SEP_new );
  // Every invarant gap size in [minSep,  maxSep) will be mapped to the array of length len,
  // and step size 1/sqrt(2), in [min, max) region
  int len = maxSep - minSep;
  ForceAssertGt( len, 0 );
  step = 1.0 / sqrt(2.0);
  min = (minSep + shift ) * step;
  max = (maxSep + shift ) * step;
  p.resize(len,0);
  d0.resize(len,0);
  d1.resize(len,0);
  d2.resize(len,0);
  // renormalize p[]
  double sum = 0;
  for(int i=0; i<len; i++)
  {
    int sep = minSep + i; // invariant sep
    p[i] = dist.prob(sep) * step;
    sum += p[i];
  }
  //for( int i=0; i < len; i++)
  //  p[i]/= sum;
  
  // calcualte d0, d1, d2
  double total0 = 0, total1= 0, total2 = 0;
  for(int i=0; i<len; i++)
  {
    double x = min + i * step;
    total0 += p[i] * step;
    total1 += p[i] * x * step;
    total2 += p[i] * x * x * step;
    d0[i] = total0;
    d1[i] = total1;
    d2[i] = total2;
  }
  //avg = total1 / total0;
  //sd = sqrt( total2 / total0 - avg * avg );
  //cout << "new avg std " <<  avg*sqrt(2.0) << " " << sd*sqrt(2.0) << endl;
}


ProbFuncIntDist:: ProbFuncIntDist( const vec<double>& array,
                               const int MIN_SEP )
{
  this->FromArray(array, MIN_SEP);
}


// Initialize the probability functon from an array
void ProbFuncIntDist:: FromArray( const vec<double>& array,
                               const int MIN_SEP )
{
  int len = array.size();
  ForceAssertGt( len, 0 );
  step = 1.0 / sqrt(2.0);
  p.resize(len,0);
  d0.resize(len,0);
  d1.resize(len,0);
  d2.resize(len,0);
  // renormalize p[]
  double sum = 0;
  for(int i=0; i<len; i++)
  {
    p[i] = array[i] * step;
    sum += p[i];
  }
  // calcualte d0, d1, d2
  double total0 = 0, total1= 0, total2 = 0;
  for(int i=0; i<len; i++)
  {
    double x = ( MIN_SEP + i )* step;
    total0 += p[i] * step;
    total1 += p[i] * x * step;
    total2 += p[i] * x * x * step;
    d0[i] = total0;
    d1[i] = total1;
    d2[i] = total2;
  }
  // assign the essential parameters 
  avg = total1 / total0;
  sd = sqrt( total2 / total0 - avg * avg );
  cout << "new avg std " <<  avg*sqrt(2.0) << " " << sd*sqrt(2.0) << endl;
  min = MIN_SEP * step;
  max = min + (double) len * step;
}


// Pair distribution from an ascii file 
ProbFuncDist::ProbFuncDist(String filename) 
{
  ForceAssert(IsRegularFile(filename));
  ifstream ifs(filename.c_str());
  string line;
  double total = 0;
  while(getline(ifs,line))
  {
    if (line[0] == '#') continue;
    istringstream iss(line);
    double x, count;
    iss>> x >> count;
    if ( count < 0 ) continue;
    if ( x < 0 ) continue;
    dist[x/sqrt(2.0)] = count;
    total += count;
  }
  // normalize
  for(map<double,double>::iterator it=dist.begin(); it != dist.end(); it++)
    it->second/=total;

  double total0 = 0, total1= 0, total2 = 0;
  for(map<double,double>::iterator it=dist.begin(); it != dist.end(); it++)
  {
    total0 += it->second;
    total1 += it->second * it->first;
    total2 += it->second * it->first * it->first;
    dist0[it->first] = total0;
    dist1[it->first] = total1;
    dist2[it->first] = total2;
    //cout << "lib " << it->first <<" "<< total0 <<" "<<total1 <<" "<<total2<<endl;
  }
  // average and standard deviation
  avg = dist1.rbegin()->second / dist0.rbegin()->second;
  sd = sqrt( dist2.rbegin()->second / dist0.rbegin()->second - avg*avg ); // < x - <x> >^2 = <x^2 > - <x>^2
}

