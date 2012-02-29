///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
#ifndef PAIR_DIST_FITTING_H
#define PAIR_DIST_FITTING_H

#include "CoreTools.h"
#include "MainTools.h"
#include "PairsManager.h"
#include "Superb.h"
#include "VecTemplate.h"
#include "VecUtilities.h"
#include "math/Functions.h"
#include "paths/PairDistModels.h"

// =====================================================
// Class ProbFuncInterface
//
// abstract class for pair distribution correction
// =====================================================
/// Find the gap (g) between two contigs (of length l1, l2 respectively) if we know 
/// a set of links that connect these two contigs, and the parent distribution that
/// describes the link length (L). Each link is given as a pair of integers (x1, x2)
///  that describe the distances from the linkage point at each contig. 
//      
//                                   | <--- g  --> |
//  contig 1 -------------------------             -------------------------------- contig 2
//                        _______________________________________
//                        | <- x1 -> |             | <-- x2 --> |
//
//                        | <---------   L -------------------- |
// Although x1 + x2 + g == L for each link,  <x1 + x2> + g == <L> is not necesarrily ture because 
// the links are from a biased subset of the library.
//
template <class T>
void MostProbableGap( const T& distr, const int len1, const int len2, const vec< pair< int, int > >& links, 
    int& gap, int& std, bool verbose=false )
{
  ForceAssertGt( links.isize(), 0 );
  // find the averge of x1 + x2 ( <x1+x2> = <L> - g in the ideal case )
  int x12avg = 0;
  for ( size_t i = 0; i < links.size(); i++ )
  {
    x12avg += links[i].first + links[i].second;
  }
  x12avg /= links.size();
  // the initial guess of the gap
  int libSep = distr.Avg() * sqrt(2.0); 
  int libSD = distr.SD() * sqrt(2.0);
  int gap0 = libSep - x12avg;
  int std0 = (double) libSD / sqrt(links.size()); 
  // iteratively find the gap size that will give a <x1 + x2>_theory matching the observation
  // Please note that <x1 + x2>_theory decrease monotonically to zero with increasing gap size,
  // and is approximately linear at small gap size.
  int CUTOFF = 5; // find solution within 5 bases
  int gapTest = gap0, gapTestOld = -1;
  int y = -1, yOld = -1;
  bool flagSignChange = false;
  int nrun=0;
  int iter_counter = 0;
  for(;nrun<20;nrun++){
    y = MeanSepObserve(distr, gapTest, len1, len2) - gapTest - x12avg;
    if ( abs( y ) < CUTOFF ) break;
    // calculate the new gapTest
    int gapTestNew = 0;
    //cout << gapTestOld<<" "<< yOld << " to " << gapTest << " " << y << endl;
    if ( nrun == 0 )
      gapTestNew = gapTest + y;
    else // when previous points are available
    {
      
      // sign change
      if ( y * yOld < 0 ) {
	flagSignChange = true;
	break;
      }
      else
      {
	double slope = abs( y - yOld ) < CUTOFF ? -1.0 : ( gapTest - gapTestOld ) / ( y - yOld );
	// step should not be too large
	int step = - y  * slope;
	if ( step > libSD ) { step = libSD; ++iter_counter; }
	if ( step < -libSD ){ step = -libSD; ++iter_counter; }
	gapTestNew = gapTest + step;
      }
    }
    // save the large run
    yOld = y;
    gapTestOld = gapTest;
    // set the new test position
    gapTest = gapTestNew;
  }
  if ( nrun == 20)
    if (verbose)
      cout << "   Warning: not converging! " << "x12avg= " << x12avg << 
	" y= " << y << " gapTest= " << gapTest << endl;
  if ( flagSignChange )
  {
    // now we know that between  gapTest and gapTestOld, the equation changes sign,
    // Find the the solution gapTestNew using binary searching.
    int gapBoundNeg = yOld < 0 ? gapTestOld : gapTest;
    int gapBoundPos = yOld < 0 ? gapTest: gapTestOld;
    int gapTestNew = -1, yNew = -1; // the test points
    int nrun2 = 0;
    for ( ; nrun2 < 20; nrun2++ ){
      gapTestNew = ( gapBoundNeg + gapBoundPos ) / 2;
      yNew = MeanSepObserve(distr, gapTestNew, len1, len2) - gapTestNew - x12avg;
      if ( abs(yNew) < CUTOFF ) break;
      if ( yNew > 0 ) 
	gapBoundPos = gapTestNew;
      else
	gapBoundNeg = gapTestNew;
      //cout << "   second round " << nrun2 << " " << gapTestNew << " " << yNew << endl;
    }
    gap = gapTestNew;
    if (verbose) cout << "   second round " << gap << " " << yNew << endl;
  }
  else
    gap = gapTest;

  gap = Min( gap, libSep + 5 * libSD + 2000 ); // gap should not be too large
  // the standard deviation
  int gapTest1 = gap + std0;
  int x12avg_theory1 = MeanSepObserve(distr, gapTest1, len1,len2) - gapTest1;
  int gapTest2 = gap - std0;
  int x12avg_theory2 = MeanSepObserve(distr, gapTest2, len1,len2) - gapTest2;
  //  (x12avg_theory1 - x12avg_theory2) / ( 2 * std0 ) is the slope in the x12avg_theory vs. gap curve
  if (  abs( x12avg_theory1 - x12avg_theory2 ) == 0 ) 
    std = libSD; // cannot estimate the uncertainty
  else 
    std = Min( libSD, 2 * std0 * std0 / abs( x12avg_theory1 - x12avg_theory2 ) );
}


// Scoring function for Kolmogorov-Smirnov test
template <class T>
double KSTestScore( const T& distr, const int len1, const int len2, const vec< int > & x12s, 
    const int gap )
{
  int xmin = x12s.front();
  int xmax = x12s.back();
  //double amin = SumProbObserve ( distr, gap, len1, len2, xmin);
  double amin = 0;
  double amax = SumProbObserve ( distr, gap, len1, len2, xmax);
  double total_prob = amax - amin;
  if ( amax <= amin ) return -1;
  //cout << "amin amax " << amin << " " << amax << endl;
  double max_diff = -1.0;
  for ( size_t i = 0; i < x12s.size() ; i++ ) {
    int x12 = x12s[ i ];
    double percent_model  = ( SumProbObserve ( distr, gap, len1, len2, x12) - amin ) / total_prob;
    double percent =  double(i) / double( x12s.size() );
    double percentUp =  double(i+1) / double( x12s.size() );
    //double diff =  abs(percent_model - percent);
    double diff = Max( percent_model - percent, percentUp - percent_model );
    max_diff = Max ( diff, max_diff ); 
    //cout << "gap" << gap << " " << x12 << " " << percent <<" " << percent_model << endl;
  }
  return max_diff;
}


// Scoring function for Kolmogorov-Smirnov test
template <class T>
double KSTestScore( const vec<T>& distrs, const int len1, const int len2, const vec< vec< int > >& x12s, 
    const int gap )
{
  unsigned int nLib = distrs.size();
  // combine all the libraries
  vec< pair<int, int> > xy;
  for ( size_t iLib = 0; iLib < nLib; iLib++ ) 
    for ( size_t i = 0; i < x12s[iLib].size() ; i++ ) 
      xy.push( make_pair( x12s[iLib][i], iLib ) );
  Sort(xy);
  // the range of the data
  int xmin = xy.front().first;
  int xmax = xy.back().first;
  double amin = 0, amax = 0;
  for ( size_t iLib = 0; iLib < nLib; iLib++ ) {
    amax += SumProbObserve ( distrs[iLib], gap, len1, len2, xmax) *  x12s[iLib].size();
  }
  double total_prob = amax - amin;
  if ( total_prob <= 0 ) return -1; // sanity-check
  // find the D value for Kolmogorov-Smirnov distribution
  double max_diff = -1.0;
  for( size_t i = 0; i< xy.size(); i++ ) { 
    int x12 = xy[i].first;
    double val = 0;
    for ( size_t iLib = 0; iLib < nLib; iLib++ ) 
      val += SumProbObserve ( distrs[iLib], gap, len1, len2, x12 ) *  x12s[iLib].size();
    double Y = ( val - amin ) / total_prob ;
    double Yi = double(i) / double( xy.size() ); // actual percentile from the data
    double Yi_up  =  double(i+1) / double( xy.size() );
    double diff = Max( Y - Yi, Yi_up - Y);
    max_diff = Max ( diff, max_diff ); 
    //cout << "gap" << gap << " " << x12 << " " << Yi <<" " << Y << endl;
  }
  return max_diff;
}

// Estimate the gap size by Kolmogorov-Smirnov fitting
template <class T>
void GapByFitKS( const T& distr, const int len1, const int len2, const vec< pair< int, int > >& links, 
    int& gap, int& std )
{
  vec<int> x12s;
  for ( size_t i = 0; i < links.size(); i++ )
    x12s.push( links[i].first + links[i].second );
  Sort( x12s );
  int LibStd = distr.SD() * sqrt(2.0);
  int step = LibStd / 3;
  int gap0 = 0;
  int t1 = gap0;
  double y1 = KSTestScore( distr, len1, len2, x12s, t1 );
  int t2 = gap0 + step;
  double y2 = KSTestScore( distr, len1, len2, x12s, t2 );
  if ( y2 > y1 )  { // t1 -> t2 step decrease the target function
    swap( y1, y2 ); 
    swap( t1, t2 );
    step = -step;
  }
  int nrun = 0;
  for(;nrun<20;nrun++) {
    if ( abs(step) < 5 ) break;
    int t3 = t2 + step;
    double y3 = KSTestScore( distr, len1, len2, x12s, t3 ); 
    if ( y3 >= y2 ) {
      step /= 3;
      t2 = t1 + step;
      y2 = KSTestScore( distr, len1, len2, x12s, t2 );
    } else {
      y1 = y2; t1 = t2;
      y2 = y3; t2 = t3;
    }
  }
  gap = t2;
}

// Estimate the gap size by Kolmogorov-Smirnov fitting of the combined distribution
template <class T>
void GapByFitKSCombine( const vec<T>& distrs, const int len1, const int len2, const vec< vec< pair<int, int> > > & links, 
    int& gap, int& std )
{
  unsigned int nLib =  distrs.size();
  vec< vec<int> > x12s( nLib );
  for( size_t iLib = 0; iLib < nLib; iLib++ )
    for( size_t i = 0; i < links[iLib].size(); i++ )
      x12s[ iLib ].push( links[iLib][i].first + links[iLib][i].second );
  for( size_t iLib = 0; iLib < nLib; iLib++ )
    Sort( x12s[iLib] );
  int gap0 = 0; // initial guess of the gap
  // minimal libStd decides the trial step for the searching
  double minLibStd = distrs[0].SD();
  for( size_t iLib = 1; iLib < nLib; iLib++ )
    minLibStd = Min ( minLibStd, distrs[iLib].SD() );
  minLibStd *= sqrt(2.0);
  int step = minLibStd / 3;
  int t1 = gap0;
  double y1 = KSTestScore( distrs, len1, len2, x12s, t1 );
  int t2 = gap0 + step;
  double y2 = KSTestScore( distrs, len1, len2, x12s, t2 );
  if ( y2 > y1 )  { // t1 -> t2 step decrease the target function
    swap( y1, y2 ); 
    swap( t1, t2 );
    step = -step;
  }
  int nrun = 0;
  for(;nrun<20;nrun++) {
    if ( abs(step) < 5 ) break;
    int t3 = t2 + step;
    double y3 = KSTestScore( distrs, len1, len2, x12s, t3 ); 
    //cout << "t3y3= " << t3 << " " << y3 << endl;
    if ( y3 >= y2 ) {
      step /= 3;
      t2 = t1 + step;
      y2 = KSTestScore( distrs, len1, len2, x12s, t2 );
    } else {
      y1 = y2; t1 = t2;
      y2 = y3; t2 = t3;
    }
  }
  //cout << "t2y2= " << t2 << " " << y2 << endl;
  gap = t2;
}


// Model=estimated the total number of counts in the selected regions
template <class T>
double ModelTotalCounts( int gap, const T& distr, const vec< pair<int, int> >& sel1, 
    const vec< pair<int,int> >& sel2, const double cr, const int max_x1x2 )
{
  double total = 0;
  for( size_t i = 0; i < sel1.size(); i++ ){
    for( size_t j = 0; j < sel2.size(); j++ ){
      int x11 = sel1[i].first;
      int x12 = sel1[i].second;
      int x21 = sel2[j].first;
      int x22 = sel2[j].second;
      int len1_cut = x12 - x11;
      int len2_cut = x22 - x21;
      int shift = x11 + x21;
      total += SumProbObserve( distr, gap + shift, len1_cut, len2_cut, max_x1x2 - shift ) * cr;
    }
  }
  return total;
}

// Model=estimated the total number of counts in the selected regions
template <class T>
double ModelTotalCountsMerge( int gap, const vec<T>& distrs, 
    const vec< pair<int, int> >& sel1, const vec< pair<int,int> >& sel2, 
    const vec<double>& cr, const vec<int>& max_x1x2s )
{
  unsigned int nLib =  distrs.size();
  // Sanity check
  ForceAssertEq( nLib, cr.size() );
  double prob = 0;
  for( size_t iLib = 0; iLib < nLib; iLib++ ) {
    for( size_t i = 0; i < sel1.size(); i++ ){
      for( size_t j = 0; j < sel2.size(); j++ ){
	int x11 = sel1[i].first;
	int x12 = sel1[i].second;
	int x21 = sel2[j].first;
	int x22 = sel2[j].second;
	int len1_cut = x12 - x11;
	int len2_cut = x22 - x21;
	int shift = x11 + x21;
	prob += SumProbObserve( distrs[iLib], gap + shift, len1_cut, len2_cut, max_x1x2s[iLib] - shift ) * cr[iLib];
      }
    }
  }
  return prob;
}


// Model=estimated average seperation between links
template <class T>
double ModelAvgSepMerge( int gap, const vec<T>& distrs, 
    const vec< pair<int, int> >& sel1, const vec< pair<int,int> >& sel2, 
    const vec< double >& cr, const vec<int>& MAX_cuts )
{
  unsigned int nLib =  distrs.size();
  ForceAssertEq( nLib, cr.size() );
  double S1 = 0, S0 = 0;
  const double scaling = 1.0 / sqrt(2.0);
  for( size_t iLib = 0; iLib < nLib; iLib++ ) {
    double cutoff =  ( MAX_cuts[iLib] + gap ) * scaling;
    for( size_t i = 0; i < sel1.size(); i++ ) {
      for( size_t j = 0; j < sel2.size(); j++ ) {
	int x11 = sel1[i].first;
	int x12 = sel1[i].second;
	int x21 = sel2[j].first;
	int x22 = sel2[j].second;
	int len1_cut = x12 - x11;
	int len2_cut = x22 - x21;
	double l1p = len1_cut * scaling; // length projected on the axis
	double l2p = len2_cut * scaling;
	double x1 = ( gap + x11 + x21 ) * scaling;
	double h = min(l1p, l2p);
	double x2 = x1 + h;
	double x4 = x1 + l1p + l2p;
	double x3 = x4 - h;
	S0 += Sum0( distrs[iLib], x1, x2, x3, x4, cutoff ) * cr[iLib];
	S1 += Sum1( distrs[iLib], x1, x2, x3, x4, cutoff ) * cr[iLib];
      }
    }
  }
  if ( S0 == 0 ) return 1000000;
  else return S1 / S0 / scaling;
}

// Model=estimated average seperation between links
template <class T>
double ModelAvgSep( int gap, const T& distr, 
    const vec< pair<int, int> >& sel1, const vec< pair<int,int> >& sel2, 
    int MAX_cut )
{
  double S1 = 0, S0 = 0;
  const double scaling = 1.0 / sqrt(2.0);
  double cutoff =  (MAX_cut + gap ) * scaling;
  for( size_t i = 0; i < sel1.size(); i++ ) {
    for( size_t j = 0; j < sel2.size(); j++ ) {
      int x11 = sel1[i].first;
      int x12 = sel1[i].second;
      int x21 = sel2[j].first;
      int x22 = sel2[j].second;
      int len1_cut = x12 - x11;
      int len2_cut = x22 - x21;
      double l1p = len1_cut * scaling; // length projected on the axis
      double l2p = len2_cut * scaling;
      double x1 = ( gap + x11 + x21 ) * scaling;
      double h = min(l1p, l2p);
      double x2 = x1 + h;
      double x4 = x1 + l1p + l2p;
      double x3 = x4 - h;
      S0 += Sum0( distr, x1, x2, x3, x4, cutoff );
      S1 += Sum1( distr, x1, x2, x3, x4, cutoff );
    }
  }
  if ( S0 == 0 ) return -1000000;
  else return S1 / S0 / scaling;
}


// estimate the gap size using the assumption < x1 + x2 > + gap = < s >
template <class T>
void GapByAvgSep( const vec<T>& distrs, const vec< pair<int, int> >& sel1, const vec< pair<int,int> >& sel2, 
    const vec< vec< pair< int, int > > >& links, 
    int& gap, int& std)
{
  unsigned int nLib =  distrs.size();
  ForceAssertEq( nLib, links.size() );
  vec< normal_distribution > gap_est;
  for( size_t iLib = 0; iLib < nLib; iLib++ ) {
    double sum = 0; 
    int count = 0;
    for( size_t k = 0; k < links[iLib].size(); k++ ) {
      int x1 = links[iLib][k].first;
      int x2 = links[iLib][k].second;
      for( size_t i = 0; i < sel1.size(); i++ ){
	for( size_t j = 0; j < sel2.size(); j++ ){
	  int x11 = sel1[i].first;
	  int x12 = sel1[i].second;
	  int x21 = sel2[j].first;
	  int x22 = sel2[j].second;
	  int shift = x11 + x21;
	  if ( x1 >= sel1[i].first && x1 < sel1[i].second && 
	      x2 >= sel2[j].first && x2 < sel2[j].second ){
	    count++;
	    sum += x1 + x2;
	  }
	}
      }
    }
    if ( count == 0 ) continue;
    double gap = distrs[iLib].Avg() * sqrt(2.0) - sum / count;
    double std = distrs[iLib].SD() * sqrt(2.0) / sqrt( double(count) );
    gap_est.push( normal_distribution( gap, std ) );
  }
  normal_distribution combine = CombineNormalDistributions(gap_est);
  gap = combine.mu_;
  std = combine.sigma_;
  return ;
}

template <class T>
void GapByFitCount( const vec<T>& distrs, const vec< pair<int, int> >& sel1, const vec< pair<int,int> >& sel2, 
    const vec< vec< pair< int, int > > >& links, const vec<double>& cr,
    int& gap, int& std)
{ 
  bool verbose = false;
  unsigned int nLib =  distrs.size();
  ForceAssertEq( nLib, cr.size() );
  int nlinks = 0;
  //int max_x1x2 = 0;
  vec<int> max_x1x2s( nLib, 0 );
  for( size_t iLib = 0; iLib < nLib; iLib++ ) 
    for( size_t i = 0; i < links[iLib].size(); i++ ) {
      int x1 = links[iLib][i].first;
      int x2 = links[iLib][i].second;
      for( size_t i = 0; i < sel1.size(); i++ )
	for( size_t j = 0; j < sel2.size(); j++ )
	  if ( x1 >= sel1[i].first && x1 < sel1[i].second && 
	      x2 >= sel2[j].first && x2 < sel2[j].second )
	  {
	    nlinks++;
	    //if ( x1 + x2 > max_x1x2 ) max_x1x2 = x1 + x2;
	    if ( x1 + x2 > max_x1x2s[iLib] ) max_x1x2s[iLib] = x1 + x2;
	  }
    }
  for( size_t iLib = 0; iLib < nLib; iLib++ ) {
    max_x1x2s[iLib]++; // [ 0, max_x1x2 )
    if (verbose) cout << iLib << "max_x1x2= " << max_x1x2s[iLib]<< endl;
  }
  if (verbose) cout << "nlinks= " << nlinks << endl;
  
  int gap0 = gap;
  int step = 1000;
  int t1 = gap0;
  double y1 = abs ( ModelTotalCountsMerge ( t1, distrs,  sel1, sel2, cr, max_x1x2s ) - nlinks );
  int t2 = gap0 + step;
  double y2 = abs ( ModelTotalCountsMerge ( t2, distrs,  sel1, sel2, cr, max_x1x2s ) - nlinks );
  if ( y2 > y1 )  { // t1 -> t2 step decrease the target function
    swap( y1, y2 ); 
    swap( t1, t2 );
    step = -step;
  }
  int nrun = 0;
  for(;nrun<20;nrun++) {
    if ( abs(step) < 5 ) break;
    int t3 = t2 + step;
    double y3 = abs ( ModelTotalCountsMerge ( t3, distrs,  sel1, sel2, cr, max_x1x2s) - nlinks );
    //cout << "t3y3= " << t3 << " " << y3 << endl;
    if ( y3 >= y2 ) {
      step /= 3;
      t2 = t1 + step;
      y2 = abs ( ModelTotalCountsMerge ( t2, distrs,  sel1, sel2, cr, max_x1x2s) - nlinks );
    } else {
      y1 = y2; t1 = t2;
      y2 = y3; t2 = t3;
    }
  }
  gap = t2;
  if ( nrun == 20 ) cout << "grep=!!GatByFitCount" << endl;
  return;
}


template <class T>
void GapByFitSepMerge( const vec<T>& distrs, const vec< pair<int, int> >& sel1, const vec< pair<int,int> >& sel2, 
    const int& mean_x1x2, const int& max_x1x2_input, const vec<double>& cr,
    int& gap, int& std)
{ 
  unsigned int nLib =  distrs.size();
  ForceAssertEq( nLib, cr.size() );
  int gap0 = gap;
  int step = 1000;
  //double max_x1x2 = sel1.back().second + sel2.back().second;
  double max_x1x2 = max_x1x2_input;
  int t1 = gap0;
  double y1 = abs( ModelAvgSepMerge( t1, distrs,  sel1, sel2, cr, max_x1x2) - t1 - mean_x1x2);
  int t2 = gap0 + step;
  double y2 = abs( ModelAvgSepMerge( t2, distrs,  sel1, sel2, cr, max_x1x2) - t2 - mean_x1x2);
  if ( y2 > y1 )  { // t1 -> t2 step decrease the target function
    swap( y1, y2 ); 
    swap( t1, t2 );
    step = -step;
  }
  int nrun = 0;
  for(;nrun<20;nrun++) {
    if ( abs(step) < 5 || y2 < 5 ) break;
    int t3 = t2 + step;
    double y3 = abs( ModelAvgSepMerge( t3, distrs,  sel1, sel2, cr, max_x1x2) - t3 - mean_x1x2);
    //cout << "t3y3= " << t3 << " " << y3 << endl;
    if ( y3 >= y2 ) {
      step /= 3;
      t2 = t1 + step;
      y2 = abs( ModelAvgSepMerge( t2, distrs,  sel1, sel2, cr, max_x1x2 ) - t2 - mean_x1x2);
      if ( y2 >= y1 ) {
	swap( y1, y2 ); 
	swap( t1, t2 );
	step = -step;
      }
    } else {
      y1 = y2; t1 = t2;
      y2 = y3; t2 = t3;
    }
  }
  gap = t2;
  if ( nrun == 20 ) cout << "grep=!!GatByFitSepMerge" << endl;
  return;
}

template <class T>
double ScoreFitSepCombine( int t, const vec<T>& distrs, const vec< pair<int, int> >& sel1,
    const vec< pair<int, int> >& sel2, const vec<int>& max_x1x2, const vec< double >& means, 
    const vec< double >& stds)
{
  double y = 0;
  for( size_t iLib = 0; iLib < distrs.size(); iLib++ ){
    if ( stds[iLib] <= 0 ) continue;
    double yi = ( ModelAvgSep( t, distrs[iLib], sel1, sel2, max_x1x2[iLib]) - t - means[iLib] ) / stds[iLib];
    y += yi * yi;
  }
  return y;
}

// Scoring function for the average separation differance
template <class T>
double ScoreFitSep( int g, const T& distr, const vec< pair<int, int> >& sel1,
    const vec< pair<int, int> >& sel2, const int max_x1x2, 
    const double mean, const double std)
{
  if ( std <= 0 ) return 0;
  double yi =  ( ModelAvgSep( g, distr, sel1, sel2, max_x1x2 ) - g - mean ) / std;
  return yi * yi;
}

// Scoring function for the count differance
template <class T>
double ScoreFitCount( int g, const T& distr, const vec< pair<int, int> >& sel1, 
    const vec< pair<int,int> >& sel2,  const int max_x1x2, double cr,
    const int nlinks, const double std ) {
  if ( std <= 0 ) return 0;
  double yi = ( ModelTotalCounts( g, distr,  sel1, sel2, cr, max_x1x2 ) - nlinks ) / std;
  return yi * yi;
}

template <class T>
double ScoreFitAll( int g, const vec<T>& distrs, const vec< pair<int, int> >& sel1,
    const vec< pair<int, int> >& sel2, const vec<int>& max_x1x2s, const vec<double> crs,
    const vec<int>& nlinks, const vec<double>& nlink_stds, 
    const vec<double>& means, const vec<double>& mean_stds )
{
  unsigned int nLib =  distrs.size();
  double s = 0;
  for( size_t iLib = 0; iLib < nLib; iLib++ ){
    if ( mean_stds[iLib] > 0 )
      s += ScoreFitSep( g, distrs[iLib], sel1, sel2, max_x1x2s[iLib], means[iLib], mean_stds[iLib] );
    if ( nlink_stds[iLib] > 0 )
      s += ScoreFitCount( g, distrs[iLib], sel1, sel2, max_x1x2s[iLib], crs[iLib], nlinks[iLib], nlink_stds[iLib] );
  }
  return s;
}

template <class T>
void GapByFitSepCombine(
    const vec<T>& distrs, 
    const vec< pair<int, int> >& sel1, 
    const vec< pair<int,int> >& sel2, 
    const vec<int>& cutoffs, 
    const vec<double>& means,
    const vec<double>& stds,
    int& gap, 
    int& std)
{ 
  unsigned int nLib = distrs.size();
  ForceAssertEq( nLib, means.size() );
  int gap0 = gap;
  int step = 1000;
  //double max_x1x2 = sel1.back().second + sel2.back().second;
  int t1 = gap0;
  double y1 = ScoreFitSepCombine( t1, distrs,  sel1, sel2, cutoffs, means, stds);
  int t2 = gap0 + step;
  double y2 = ScoreFitSepCombine( t2, distrs,  sel1, sel2, cutoffs, means, stds);
  if ( y2 > y1 )  { // t1 -> t2 step decrease the target function
    swap( y1, y2 ); 
    swap( t1, t2 );
    step = -step;
  }
  int nrun = 0;
  for(;nrun<20;nrun++) {
    if ( abs(step) < 5 || y2 < 5 ) break;
    int t3 = t2 + step;
    double y3 = ScoreFitSepCombine( t3, distrs,  sel1, sel2, cutoffs, means, stds);
    //cout << "t3y3= " << t3 << " " << y3 << endl;
    if ( y3 >= y2 ) {
      step /= 3;
      t2 = t1 + step;
      y2 = ScoreFitSepCombine( t2, distrs,  sel1, sel2, cutoffs, means, stds);
      if ( y2 >= y1 ) {
	swap( y1, y2 ); 
	swap( t1, t2 );
	step = -step;
      }
    } else {
      y1 = y2; t1 = t2;
      y2 = y3; t2 = t3;
    }
  }
  gap = t2;
  if ( nrun == 20 ) cout << "grep=!!GatByFitSepCombine" << endl;
  return;
}

template <class T>
void GapByFitSep( const T distr, const vec< pair<int, int> >& sel1, const vec< pair<int,int> >& sel2, 
    const int& mean_x1x2, const int& max_x1x2_input, 
    int& gap, int& std)
{ 
  int gap0 = gap;
  int step = 1000;
  //double max_x1x2 = sel1.back().second + sel2.back().second;
  double max_x1x2 = max_x1x2_input;
  int t1 = gap0;
  double y1 = abs( ModelAvgSep( t1, distr,  sel1, sel2, max_x1x2) - t1 - mean_x1x2);
  int t2 = gap0 + step;
  double y2 = abs( ModelAvgSep( t2, distr,  sel1, sel2, max_x1x2) - t2 - mean_x1x2);
  if ( y2 > y1 )  { // t1 -> t2 step decrease the target function
    swap( y1, y2 ); 
    swap( t1, t2 );
    step = -step;
  }
  int nrun = 0;
  for(;nrun<20;nrun++) {
    int t3 = t2 + step;
    double y3 = abs( ModelAvgSep( t3, distr,  sel1, sel2, max_x1x2) - t3 - mean_x1x2);
    if ( abs(step) < 5 || t3 < 5 ) break;
    if ( y3 >= y2 ) {
      step /= 3;
      t2 = t1 + step;
      y2 = abs( ModelAvgSep( t2, distr,  sel1, sel2, max_x1x2 ) - t2 - mean_x1x2);
      if ( y2 >= y1 ) {
	swap( y1, y2 ); 
	swap( t1, t2 );
	step = -step;
      }
    } else {
      y1 = y2; t1 = t2;
      y2 = y3; t2 = t3;
    }
  }
  gap = t2;
  if ( nrun == 20 ) cout << "grep= !! GatByFitSep" << endl;
  return;
}


template <class T>
void GapByFitAllCombine(
    const vec<T>& distrs, 
    const vec< pair<int, int> >& sel1, 
    const vec< pair<int,int> >& sel2, 
    const vec<int>& cutoffs, 
    const vec<double>& crs,
    const vec< vec< pair<int,int> > >& links,
    int& gap, 
    int& std)
{ 
  bool verbose = false;
  unsigned int nLib = distrs.size();
  ForceAssertEq( nLib, links.size() );
  ForceAssertEq( nLib, cutoffs.size() );
  // get all x1x2 of the links, in the selected region only.
  vec< int > x1x2s;
  vec< vec<int> > x1x2s_lib(nLib);
  for( size_t iLib = 0; iLib < nLib; iLib++ )
    for( size_t k = 0; k < links[iLib].size(); k++ ) {
      int x1 = links[iLib][k].first;
      int x2 = links[iLib][k].second;
      for( size_t i = 0; i < sel1.size(); i++ )
	for( size_t j = 0; j < sel2.size(); j++ )
	  if ( x1 >= sel1[i].first && x1 < sel1[i].second && 
	      x2 >= sel2[j].first && x2 < sel2[j].second ) {
	    int xsum = x1 + x2;
	    x1x2s_lib[iLib].push( xsum );
	    x1x2s.push( xsum );
	  }
    }

  // number of links
  vec<int> nlinks(nLib);
  vec<double> nlink_stds(nLib);
  for( size_t i = 0; i < nLib; i++) {
    nlinks[i] = x1x2s_lib[i].size();
    nlink_stds[i] = Max( sqrt( nlinks[i] + 1.0), (nlinks[i] + 1.0) * 0.5 ); // observed 50% coverage flutuation
  }
  int nlink_all = x1x2s.size();
  if ( nlink_all == 0 ) return; 

  // mean and std of x1 + x2 by lib
  vec<double> means(nLib, 0), mean_stds(nLib, -1); // no points if mean_stds[i] < 0
  vec< normal_distribution > mean_x1x2;
  for( size_t iLib = 0; iLib < nLib; iLib++) {
    if ( x1x2s_lib[iLib].size() == 0 ) continue;
    means[iLib] = Mean( x1x2s_lib[iLib] );
    double n = x1x2s_lib[iLib].size();
    mean_stds[iLib] = distrs[iLib].SD() * sqrt(2.0) / sqrt(n);
    mean_x1x2.push( normal_distribution( means[iLib], mean_stds[iLib] ) );
  }

  // debug
  for( size_t iLib = 0; iLib < nLib; iLib++) {
    if (verbose) cout << "lib " << iLib << " " << int(distrs[iLib].Avg()*sqrt(2.0)) << " " 
      <<  int(distrs[iLib].SD()*sqrt(2.0)) << " " <<  nlinks[iLib] << " " << means[iLib] << " " << mean_stds[iLib] << endl;
  }

  // the mean, std of the combined distribution
  double std_all = 0, mean_all = 0;
  mean_all = Mean( x1x2s );
  std_all = StdDev( x1x2s, mean_all ) / sqrt( double(nlink_all) );
  if (verbose) cout << "mean_all= " << mean_all << endl;
  if (verbose) cout << "std_all= " << std_all << endl;
  normal_distribution nd = CombineNormalDistributions( mean_x1x2 );
  mean_all = nd.mu_;
  std_all = nd.sigma_;
  if (verbose) cout << "mean_all= " << mean_all << endl;
  if (verbose) cout << "std_all= " << std_all << endl;

  // the initial guess
  int gap0, gap0_std;
  GapByAvgSep( distrs, sel1, sel2, links,  gap0, gap0_std );
  double diff = ModelAvgSepMerge( gap0, distrs, sel1, sel2, crs, cutoffs ) - gap0 - mean_all;
  double diff2 = ModelTotalCountsMerge( gap0, distrs, sel1, sel2, crs, cutoffs ) - nlink_all;
  double range2 = Max(  sqrt(nlink_all+1.0), (nlink_all + 1.0) * 0.5 ); // allowed fluctuation 
  
  // debug
  if (verbose) cout << "gap0= " << gap0 << endl;
  if (verbose) cout << "diff= " << diff << endl;
  String str1 = abs(diff) < std_all ? "avg" : "";
  String str2 = abs(diff2) < range2 ? "cnt" : "";
  if (verbose) cout << "grep0" << str1 << str2 << " " << endl;
  
  // initial guess looks ok
  if ( abs(diff) < std_all && abs(diff2) < range2 ) {
    gap = gap0;
    std = gap0_std;
    return;
  }

  // only mean sep fit
  int gap2 = gap0, gap2_std;
  GapByFitSepCombine( distrs, sel1, sel2, cutoffs, means, mean_stds, gap2, gap2_std );

  //int step = std_all * sqrt(nlink_all);
  int step = 1000;
  int t1 = gap0;
  double y1 = ScoreFitAll( t1, distrs, sel1, sel2, cutoffs, crs, nlinks, nlink_stds, means, mean_stds );
  int t2 = gap0 + step;
  double y2 = ScoreFitAll( t2, distrs, sel1, sel2, cutoffs, crs, nlinks, nlink_stds, means, mean_stds );
  if ( y2 > y1 )  { // t1 -> t2 step decrease the target function
    swap( y1, y2 ); 
    swap( t1, t2 );
    step = -step;
  }
  int nrun = 0;
  for(;nrun<20;nrun++) {
    if ( abs(step) < 5 || y2 < 0.01 ) break;
    int t3 = t2 + step;
    double y3 = ScoreFitAll( t3, distrs, sel1, sel2, cutoffs, crs, nlinks, nlink_stds, means, mean_stds );
    if ( y3 >= y2 ) {
      step /= 3;
      t2 = t1 + step;
      y2 = ScoreFitAll( t2, distrs, sel1, sel2, cutoffs, crs, nlinks, nlink_stds, means, mean_stds );
      if ( y2 >= y1 ) {
	swap( y1, y2 ); 
	swap( t1, t2 );
	step = -step;
      }
    } else {
      y1 = y2; t1 = t2;
      y2 = y3; t2 = t3;
    }
  }
  gap = t2;
  if ( nrun == 20 ) cout << "grep!!GatByFitSep" << endl;

  if (verbose) cout << "grep= " << gap2 << " " 
    << int(diff) << " " << int(std_all) << " " 
    << int(diff2) << " " << " " << int(range2) << " " 
    << y2 << " " << nlink_all << " " << nrun << endl;

  return;
}


#endif
