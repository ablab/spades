///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
#ifndef PAIR_DIST_CORRECTION_H
#define PAIR_DIST_CORRECTION_H

#include "CoreTools.h"
#include "MainTools.h"
#include "PairsManager.h"
#include "Superb.h"
#include "VecTemplate.h"
#include "VecUtilities.h"
#include "math/Functions.h"
#include "paths/Alignlet.h"
#include <map>
#include "math/IntDistribution.h"
#include "paths/PairDistModels.h"
#include "paths/PairDistFitting.h"

// Find consecutive positive coverage regions. Watch for short contigs where the coverage may drop near the end
void inline FindValidRegion( const vec<int>& total_cr, vec< pair<int,int> >& region, unsigned int rep_cutoff )
{
  int MIN_COUNT = 5;
  int left = -1, right = -1, counter = 0;
  for( size_t p = 0; p < total_cr.size(); p++ ) {
    if ( total_cr[p] == 0 && p < total_cr.size() - 1) continue;
    counter ++;
    if ( left < 0 ) left = right = p; // first encounter
    // if find positive coverage or at the end
    if ( p - right > rep_cutoff ) {
      if ( counter >= MIN_COUNT )
	region.push( make_pair( left, right+1 ) );
      left = p;
      counter = 1;
    } else if ( p == total_cr.size() - 1 ) 
      region.push( make_pair( left, p+1 ) ) ;
    right = p;
  }
}

// clone of fitgap3 just for test
template <class T>
void FitGap3Test( const vec<T>& distrs, const int len1, const int len2, const vec< vec< pair<int, int> > > & links, 
    const vec< vec<int> >& x1s, const vec< vec<int> >& x2s, const vec<int> cutoffs,
    int& gap, int& std )
{
  unsigned nLib = distrs.size();
  // Sanity check
  ForceAssertEq( nLib, links.size() );
  ForceAssertEq( nLib, x1s.size() );
  ForceAssertEq( nLib, x2s.size() );
  // Calculate the average coverage for the libraries. Also identify the 
  // repetative region, for which the links has been filtered out ( low coverage ).
  unsigned int cl1 = Min( 20000, len1 ); // donot sample more than 10000 bases
  unsigned int cl2 = Min( 20000, len2 ); 
  cout << "cl1= " << cl1 << endl;
  cout << "cl2= " << cl2 << endl;
  vec<int> total_cr1( cl1, 0 ); // totoal coverage
  vec<int> total_cr2( cl2, 0 ); 
  // Note that the numbering will be inverted for contig1 coverage x1s
  // Also any coverage for p < 96 is zero for using invarant pair separation
  for( size_t iLib = 0; iLib < nLib; iLib++ ) {
    for( size_t p = 96; p < x1s[iLib].size() && p < cl1; p++ ) 
      total_cr1[p] += x1s[iLib][x1s[iLib].size() -1 - p];
    for( size_t p = 96; p < x2s[iLib].size() && p < cl2; p++ ) 
      total_cr2[p] += x2s[iLib][p];
  }
  // find the region with positive coverage
  vec< pair<int, int> > region1;
  vec< pair<int, int> > region2;
  unsigned int rep_cutoff = 90;
  FindValidRegion( total_cr1, region1, rep_cutoff );
  FindValidRegion( total_cr2, region2, rep_cutoff );
  ForceAssertGt( region1.size(), 0u );
  ForceAssertGt( region2.size(), 0u );
  {
    for( size_t i = 0; i< region1.size(); i++)
      cout << "rep1= " << region1[i].first << " - " << region1[i].second << endl;
    for( size_t i = 0; i< region2.size(); i++)
      cout << "rep2= " << region2[i].first << " - " << region2[i].second << endl;
  }

  // find the average coverage over the selected regions
  vec<double> avg_cov1( nLib, 0 );
  vec<double> avg_cov2( nLib, 0 );
  for( size_t iLib = 0; iLib < nLib; iLib++ ) {
    int sum = 0, counter = 0;
    for( size_t i = 0; i < region1.size(); i++ ){
      for ( int p = region1[i].first; p < region1[i].second; p++) {
	sum += x1s[iLib][ x1s[iLib].size() -1  - p ];
	counter++;
      }
    }
    avg_cov1[iLib] = double(sum) / counter;
  }
  for( size_t iLib = 0; iLib < nLib; iLib++ ) {
    int sum = 0, counter = 0;
    for( size_t i = 0; i < region2.size(); i++ ){
      for ( int p = region2[i].first; p < region2[i].second; p++) {
	sum += x2s[iLib][p];
	counter++;
      }
    }
    avg_cov2[iLib] = double(sum) / counter;
  }
  // average coverage excluding invalid regions
  vec<double> cr( nLib, 0 );
  for( size_t iLib = 0; iLib < nLib; iLib++ ) {
    cr[iLib] = ( avg_cov1[iLib] + avg_cov2[iLib] ) / 4.0;  // also consider 
    cout << "cr"<<iLib << " " <<  avg_cov1[iLib] << "," << avg_cov2[iLib] << endl;
  }

  // vec< int > x1x2s;
  // vec< vec<int> > x1x2s_lib(nLib);
  // for( size_t iLib = 0; iLib < nLib; iLib++ )
  //   for( size_t k = 0; k < links[iLib].size(); k++ ) {
  //     int x1 = links[iLib][k].first;
  //     int x2 = links[iLib][k].second;
  //     for( size_t i = 0; i < region1.size(); i++ )
  //       for( size_t j = 0; j < region2.size(); j++ )
  //         if ( x1 >= region1[i].first && x1 < region1[i].second && 
  //             x2 >= region2[j].first && x2 < region2[j].second ) {
  //           int xsum = x1 + x2;
  //           x1x2s.push( xsum );
  //           x1x2s_lib[iLib].push( xsum );
  //         }
  //   }
  // ForceAssertGt( x1x2s.size(), 0u );
  // int mean_x1x2 = Mean(x1x2s);
  // int max_x1x2 = Max(x1x2s) + 1;
  // //int max_x1x2 = region1.back().second + region2.back().second;
  // int sd_x1x2 = StdDev(x1x2s, mean_x1x2);
  // int nlinks = x1x2s.size();
  //cout << "mean_x1x2= " << mean_x1x2 << endl;
  //cout << "max_x1x2= " << max_x1x2 << endl;
  //cout << "nlinks= " << nlinks << endl;

  // initial guess of the gap by calculating the average
  int gap0 = -10000, gap0_std = -10000;
  GapByAvgSep( distrs, region1, region2, links,  gap0, gap0_std );
  
  int gap1 = gap0, gap1_std = -10000;
  // Fit the gap size by the total number of points
  GapByFitCount( distrs, region1, region2, links, cr, gap1, gap1_std );

  //vec< normal_distribution > tmp;
  //for( size_t iLib = 0; iLib < nLib; iLib++ ){
  //  int gap2 = gap0, gap2_std = -10000;
  //  if ( x1x2s_lib[iLib].size() < 5 ) continue;
  //  GapByFitSep( distrs[iLib], region1, region2, Mean(x1x2s_lib[iLib]), Max(x1x2s_lib[iLib]),  gap2, gap2_std );
  //  double std0 = distrs[iLib].SD() * sqrt(2.0) / sqrt( x1x2s_lib[iLib].size() ); // rough estimation of error
  //  if ( gap2 > -10000 )
  //    tmp.push( normal_distribution( gap2, std0 ) );
  //}
  //normal_distribution tmp2 = CombineNormalDistributions( tmp );
  //int gap2 = tmp2.mu_;
  //int gap2_std = tmp2.sigma_;
  
  int gap2 = gap0, gap2_std = -10000;
  GapByFitAllCombine( distrs, region1, region2, cutoffs, cr, links, gap2, gap2_std );
  cout << "grep= " << gap1 << " " << gap0 << " " << gap2 << endl;
  gap = gap2;

  return;
}


// Estimate the gap size by fitting the pair separation distributions to the observed links, and the coverage 
// links[libId][pairId] = pair(x1,x2),  x1s[libId] and x2s[libId] are the coverage of over the contigs
template <class T>
void FitGap3( const vec<T>& distrs, const int len1, const int len2, const vec< vec< pair<int, int> > > & links, 
    const vec< vec<int> >& x1s, const vec< vec<int> >& x2s, const vec<int> cutoffs,
    int& gap, int& std, bool verbose = false )
{
  unsigned nLib = distrs.size();
  // Sanity check
  ForceAssertEq( nLib, links.size() );
  ForceAssertEq( nLib, x1s.size() );
  ForceAssertEq( nLib, x2s.size() );
  // Calculate the average coverage for the libraries. Also identify the 
  // repetative region, for which the links has been filtered out ( low coverage ).
  unsigned int cl1 = Min( 20000, len1 ); // donot sample more than 10000 bases
  unsigned int cl2 = Min( 20000, len2 ); 
  if (verbose) cout << "cl1= " << cl1 << endl;
  if (verbose) cout << "cl2= " << cl2 << endl;
  if (verbose) for( size_t iLib=0;iLib<nLib; iLib++ ) cout << "x1s.size()= " << x1s[iLib].size() << endl;
  if (verbose) for( size_t iLib=0;iLib<nLib; iLib++ ) cout << "x2s.size()= " << x2s[iLib].size() << endl;
  vec<int> total_cr1( cl1, 0 ); // totoal coverage
  vec<int> total_cr2( cl2, 0 ); 
  // Note that the numbering will be inverted for contig1 coverage x1s
  // Also any coverage for p < 96 is zero for using invarant pair separation
  for( size_t iLib = 0; iLib < nLib; iLib++ ) {
    for( size_t p = 96; p < x1s[iLib].size() && p < cl1; p++ ) 
      total_cr1[p] += x1s[iLib][p];
    for( size_t p = 96; p < x2s[iLib].size() && p < cl2; p++ ) 
      total_cr2[p] += x2s[iLib][p];
  }
  // find the region with positive coverage
  vec< pair<int, int> > region1;
  vec< pair<int, int> > region2;
  unsigned int rep_cutoff = 90;
  FindValidRegion( total_cr1, region1, rep_cutoff );
  FindValidRegion( total_cr2, region2, rep_cutoff );

  if (  region1.empty() ||  region2.empty() ) return;
  // output
  {
    for( size_t i = 0; i< region1.size(); i++)
      if (verbose) cout << "rep1= " << region1[i].first << " - " << region1[i].second << endl;
    for( size_t i = 0; i< region2.size(); i++)
      if (verbose) cout << "rep2= " << region2[i].first << " - " << region2[i].second << endl;
  }

  // find the average coverage over the selected regions
  vec<double> avg_cov1( nLib, 0 );
  vec<double> avg_cov2( nLib, 0 );
  for( size_t iLib = 0; iLib < nLib; iLib++ ) {
    if (  x1s[iLib].empty() ) continue; // somtime coverage is zero, so no data stored
    int sum = 0, counter = 0;
    for( size_t i = 0; i < region1.size(); i++ ){
      for ( int p = region1[i].first; p < region1[i].second; p++) {
	sum += x1s[iLib][ p ];
	counter++;
      }
    }
    if ( counter > 0 )
      avg_cov1[iLib] = double(sum) / counter;
  }
  for( size_t iLib = 0; iLib < nLib; iLib++ ) {
    if (  x2s[iLib].empty() ) continue; // somtime coverage is zero, so no data stored
    int sum = 0, counter = 0;
    for( size_t i = 0; i < region2.size(); i++ ){
      for ( int p = region2[i].first; p < region2[i].second; p++) {
	sum += x2s[iLib][p];
	counter++;
      }
    }
    if ( counter > 0 )
      avg_cov2[iLib] = double(sum) / counter;
  }
  // average coverage excluding invalid regions
  vec<double> cr( nLib, 0 );
  for( size_t iLib = 0; iLib < nLib; iLib++ ) {
    cr[iLib] = ( avg_cov1[iLib] + avg_cov2[iLib] ) / 4.0;  // also consider 
    if (verbose) cout << "cr"<<iLib << " " <<  avg_cov1[iLib] << "," << avg_cov2[iLib] << endl;
  }

  // initial guess of the gap by calculating the average
  int gap0 = -10000, gap0_std = -10000;
  GapByAvgSep( distrs, region1, region2, links,  gap0, gap0_std );
  
  int gap1 = gap0, gap1_std = -10000;
  // Fit the gap size by the total number of points
  GapByFitCount( distrs, region1, region2, links, cr, gap1, gap1_std );

  //vec< normal_distribution > tmp;
  //for( size_t iLib = 0; iLib < nLib; iLib++ ){
  //  int gap2 = gap0, gap2_std = -10000;
  //  if ( x1x2s_lib[iLib].size() < 5 ) continue;
  //  GapByFitSep( distrs[iLib], region1, region2, Mean(x1x2s_lib[iLib]), Max(x1x2s_lib[iLib]),  gap2, gap2_std );
  //  double std0 = distrs[iLib].SD() * sqrt(2.0) / sqrt( x1x2s_lib[iLib].size() ); // rough estimation of error
  //  if ( gap2 > -10000 )
  //    tmp.push( normal_distribution( gap2, std0 ) );
  //}
  //normal_distribution tmp2 = CombineNormalDistributions( tmp );
  //int gap2 = tmp2.mu_;
  //int gap2_std = tmp2.sigma_;
  
  int gap2 = gap0, gap2_std = -10000;
  GapByFitAllCombine( distrs, region1, region2, cutoffs, cr, links, gap2, gap2_std );
  if (verbose) cout << "grep= " << gap1 << " " << gap0 << " " << gap2 << endl;

  gap = gap2;
  return;
}



/// Make correction to the pair separations based on their distribution over the gaps
template <class T>
void PairDistCorrectionOld(const vec<T>& libDist, const PairsManager& pairs, const vec<alignlet>& aligns, const vec<int>& index, 
   vec<int>& seps, vec<int>& sds, bool verbose=false)
{

  int MIN_N_LINKS = 5; // minimal number of links required to apply the correction
  cout << Date() << ": Start sorting" << endl;
  map<int, map<int, map<int, vec<longlong> > > > maps;  // maps[libID][index1][index2] returns a vector of pairIDs
  size_t nLibs = pairs.nLibraries();
  ForceAssertEq(nLibs, libDist.size());
  int counter3=0;
  for( size_t i=0; i < pairs.nPairs(); i++ ) 
  {
    int aid1 = index[pairs.ID1(i)];
    int aid2 = index[pairs.ID2(i)];
    int libID = pairs.libraryID(i);
    if( aid1 < 0 || aid2 < 0 ) continue;
    int cg1 = aligns[aid1].TargetId();
    int cg2 = aligns[aid2].TargetId();
    if ( cg1 == cg2 ) continue;
    if ( cg1 > cg2 ) swap(cg1,cg2), swap(aid1,aid2);

    // map to cg1 to cg2 , where cg1 < cg2
    int index1,index2;
    if (aligns[aid1].Fw1()) 
      index1 = 2*cg1;
    else
      index1 = 2*cg1 +1;
    if (aligns[aid2].Fw1()) 
      index2 = 2*cg2;
    else
      index2 = 2*cg2 +1;
    maps[libID][index1][index2].push(i);
    counter3++;
  }
  cout << Date() << ": End sorting" << endl;
  if (verbose) cout << "            number of scaffold links= " << counter3 << endl;

  cout << Date() << ": Start converting " << endl;
  int counter2 = 0;
  for( map<int,map<int,map<int,vec<longlong> > > >::iterator it0 = maps.begin(); it0 != maps.end();
      it0++ )
  {
    int libId  = it0->first;
    int libSep = pairs.getLibrarySep(libId);
    int libSD = pairs.getLibrarySD(libId);

    for( map<int,map<int,vec<longlong> > >:: iterator it1 = it0->second.begin(); it1 !=  it0->second.end();
	it1++ )
    {
      size_t i = it1->first;
      for( map<int,vec<longlong> >:: iterator it2 = it1->second.begin(); it2 !=  it1->second.end();
	  it2++ )
      {
	size_t j = it2->first;
	unsigned int nLinks = it2->second.size();
	if (nLinks < 1 ) continue;
	int cg1 = i/2;
	int cg2 = j/2;
	bool fw1 = ( i % 2 == 0 );
	bool fw2 = ( j % 2 == 0 );
	int len1 = -1, len2 = -1;
	vec<int> gap0s;
	vec< pair<int, int> > links;
	for(size_t k=0; k < nLinks; k++) 
	{
	  longlong pairId = maps[libId][i][j][k];

	  int aid1 = index[pairs.ID1(pairId)];
	  int aid2 = index[pairs.ID2(pairId)];
	  int cg1Target = aligns[aid1].TargetId();
	  int cg2Target = aligns[aid2].TargetId();
	  ForceAssertNe( cg1Target, cg2Target );
	  if ( cg1Target > cg2Target ) swap(aid1, aid2), swap(cg1Target, cg2Target);

	  ForceAssertEq( cg1, cg1Target);
	  ForceAssertEq( cg2, cg2Target);
	  const alignlet& a1 = aligns[aid1];
	  const alignlet& a2 = aligns[aid2];
	  int len1a = a1.TargetLength();
	  int len2a = a2.TargetLength();
	  if ( len1 < 0 ) len1 = len1a;
	  if ( len2 < 0 ) len2 = len2a;
	  ForceAssertEq(len1, len1a);
	  ForceAssertEq(len2, len2a);

	  int sep = seps[pairId]; // the original separation
	  int dev = sds[pairId];  // the original deviation
	  int dist_to_end1, dist_to_end2;
	  // End-to-end distance between the reads
	  if (fw1) 
	    dist_to_end1 = len1 - a1.Pos2();
	  else
	    dist_to_end1 = a1.pos2();
	  if (fw2) 
	    dist_to_end2 = len2 - a2.Pos2();
	  else 
	    dist_to_end2 = a2.pos2();
	  //// Now use the distance from the invariant end of the jump reads
	  //if (fw1) 
	  //  dist_to_end1 = len1 - a1.pos2();
	  //else
	  //  dist_to_end1 = a1.Pos2();
	  //if (fw2) 
	  //  dist_to_end2 = len2 - a2.pos2();
	  //else 
	  //  dist_to_end2 = a2.Pos2();
	  if ( dist_to_end1 + dist_to_end2 > libSep + libSD * 5 + 1000) continue;
	  int gap0 = sep  - dist_to_end1 - dist_to_end2;
	  gap0s.push(gap0);
	  links.push( make_pair( dist_to_end1, dist_to_end2) );
	}
	if(gap0s.isize() < MIN_N_LINKS) continue;
	// initial guess of gap size and standard deviation
	int gap0 = Mean(gap0s);
	int gap0new  = -1000000, std0new = -1000000;
	MostProbableGap( libDist[libId], len1, len2, links, gap0new, std0new, verbose );
	//gap0new = Prob3(libDist[libId], gap0, len1, len2, libSep, libSD);
	int shift = gap0new - gap0;
	// !!!! hard-coded parameter
	if ( abs(shift) < libSD/3) continue;
	for(size_t k=0; k < nLinks; k++) 
	{
	  longlong pairId = maps[libId][i][j][k];
	  seps[pairId] = pairs.sep(pairId) + shift;
	  //sds[pairId] = std0new;
	  counter2++;
	}
	if (verbose) 
	  cout << "c"<< cg1 << "," << cg2 << ( fw1 ? '+':'-') <<  (fw2 ? '+' : '-')
	       << " lib "<< libId << ", " << nLinks << " links, shift= " << shift << endl;
      }
    }
  }
  if (verbose) cout << "Total pair separation adjustments= " << counter2 << endl;
  cout << Date() << ": End converting " << endl;
}


// Make correction to the pair separations based on their distribution over the gaps
// Note the we use PairsManager only to get the pairing information. The input pairing searations are pass 
// in seps and sds vectors. In this way, we have better control of the behavior of the function.
template <class T>
void PairDistCorrection(const vec<T>& libDist, const PairsManager& pairs, const vec<alignlet>& aligns, const vec<int>& index, 
   vec<int>& seps, vec<int>& sds, bool verbose=false)
{

  int MIN_N_LINKS = 5; // minimal number of links required to apply the correction
  cout << Date() << ": Start sorting" << endl;
  size_t nLib = pairs.nLibraries();
  ForceAssertEq(nLib, libDist.size());

  vec< int > cutoffs(nLib); // cutoff for dist_to_end1 + dist_to_end2;
  for( size_t libId = 0; libId < nLib; libId++ ) {
    int libSep = pairs.getLibrarySep(libId);
    int libSD = pairs.getLibrarySD(libId);
    cutoffs[libId] = libSep + libSD * 5 + 2000;
  }

  // xs[index1][libId][pos] coverage counter
  map< int, vec< vec< int > > > xs;  
  // maps[index1][index2][libID] returns an array of pairs of links
  typedef vec< pair<int, int> > Pairs;
  map<int, map<int, vec< vec< pair<int, int> > > > > link_maps;  
  map<int, map<int, vec< vec< longlong > > > > pair_maps;  
  int counter3=0;
  for( size_t iPair=0; iPair < pairs.nPairs(); iPair++ ) 
  {
    int aid1 = index[pairs.ID1(iPair)];
    int aid2 = index[pairs.ID2(iPair)];
    int libID = pairs.libraryID(iPair);
    // assign the coverage 
    vec <int> aids(2, -1);
    aids[0] = aid1; aids[1] = aid2;
    vec <int> index(2, -1); // 2*cg or 2*cg + 1 depending on direction
    vec <int> dist_to_end(2, -1);
    vec <int> tigs(2, -1);
    // pairs sorting
    for( size_t kk = 0; kk < 2; kk++ ){
      int aid = aids[kk];
      if( aid < 0) continue;
      const alignlet & a = aligns[aid];
      int cg = a.TargetId();
      tigs[kk] = cg;
      xs[2*cg].resize(nLib);
      xs[2*cg+1].resize(nLib);
      int len = a.TargetLength();
      int cl = Min( len, 20000 ); // maximum 20000 base for coverage sampling
      int p = a.Fw1() ? a.pos2() : a.Pos2() - 1;
      xs[2*cg][libID].resize(cl,0);
      xs[2*cg+1][libID].resize(cl,0);
      if ( a.Fw1() ) { // fw alignment on contig
      //          ------->
      // __________________________   
      //          <---- x1 ------->
	index[kk] = 2 * cg;
	dist_to_end[kk] = len - 1 - a.pos2();
	int dist_to_begin = a.pos2();
	if ( dist_to_end[kk] < cl ) {
	  xs[2*cg][libID][dist_to_end[kk]]++;
	}
	if ( dist_to_begin < cl ) {
	  xs[2*cg+1][libID].resize(cl,0);
	  xs[2*cg+1][libID][dist_to_begin]++;
	}
      } else { // rc alignment on contig
      //          <-------
      // __________________________   
      // <-- x1---------->
	index[kk] = 2 * cg + 1;
	dist_to_end[kk] = a.Pos2() - 1;
	int dist_to_begin = len - a.Pos2();
	if ( dist_to_end[kk] <  cl ) {
	  xs[2*cg+1][libID][dist_to_end[kk]]++;
	}
	if ( dist_to_begin < cl ) {
	  xs[2*cg][libID][dist_to_begin]++;
	}
      }
    }
    if ( tigs[0] < 0 || tigs[1] < 0 || tigs[0] == tigs[1] ) continue;
    if ( dist_to_end[0] + dist_to_end[1] > cutoffs[libID] ) continue;
    // save links from index1 to index2, where index1 < index2
    if ( index[0] < index[1] ) {
      link_maps[index[0]][index[1]].resize(nLib);
      link_maps[index[0]][index[1]][libID].push( 
	  make_pair(dist_to_end[0], dist_to_end[1]) );
      pair_maps[index[0]][index[1]].resize(nLib);
      pair_maps[index[0]][index[1]][libID].push(iPair);
    }else{
      link_maps[index[1]][index[0]].resize(nLib);
      link_maps[index[1]][index[0]][libID].push( 
	  make_pair(dist_to_end[1], dist_to_end[0]) );
      pair_maps[index[1]][index[0]].resize(nLib);
      pair_maps[index[1]][index[0]][libID].push(iPair);
    }
    counter3++;
  }

  cout << Date() << ": End sorting" << endl;
  if (verbose) cout << "            number of scaffold links= " << counter3 << endl;

  cout << Date() << ": Start converting " << endl;
  int counter2 = 0;
  for( map<int,map<int, vec<Pairs > > >::iterator it0 = link_maps.begin(); it0 != link_maps.end();
      it0++ )
  {
    int i = it0->first;
    for( map<int,vec< Pairs > >:: iterator it1 = it0->second.begin(); it1 !=  it0->second.end();
	it1++ )
    {
      size_t j = it1->first;
      int cg1 = i/2;
      int cg2 = j/2;
      //cout << "i= " << i << " ";
      //cout << "j= " << j << endl;
      //cout << "cg1= " << cg1 << " ";
      //cout << "cg2= " << cg2 << endl;
      bool fw1 = ( i % 2 == 0 );
      bool fw2 = ( j % 2 == 0 );
      vec< vec< pair<int,int> > >& all_links = it1->second;
      vec< vec< longlong > >& pair_ids = pair_maps[i][j];
      int total_links = 0;
      for ( size_t iLib = 0; iLib < nLib; iLib++ ) total_links += all_links[iLib].size();
      if( total_links < MIN_N_LINKS ) continue;
      //cout << "total_links= " << total_links << endl;
      unsigned int len1 = 0, len2 = 0; // contig length are obtained from the alignment to avoid extra arguments passing in
      for ( size_t iLib = 0; iLib < nLib; iLib++ ) {
	if ( xs[i][iLib].size() > len1 ) len1 = xs[i][iLib].size();
	if ( xs[j][iLib].size() > len2 ) len2 = xs[j][iLib].size();
      }
      if ( len1 == 0 || len2 == 0 ) continue; // nothing you can do with zero coverage
      // initial guess
      vec< vec<int> > gap0_all( nLib );
      for ( size_t iLib = 0; iLib < nLib; iLib++ ) 
	for ( size_t k=0; k < all_links[ iLib ].size(); k++ ) {
	  gap0_all[ iLib ].push( seps[ pair_ids[iLib][k] ] - all_links[iLib][k].first - all_links[iLib][k].second );
	}
      int gap_eval  = -10000, std_eval = -10000;
      bool vb = false;
      FitGap3( libDist, len1, len2, all_links, xs[i], xs[j], cutoffs, gap_eval, std_eval, vb);
      if ( gap_eval >  -10000 ) {
	for ( size_t iLib = 0; iLib < nLib; iLib++ ) {
	  int nlinks =  gap0_all[iLib].size();
	  if ( nlinks < 1 ) continue;
	  int gap0 = Mean( gap0_all[iLib]); // the original gap estimatin
	  int shift = gap_eval - gap0;
	  int libSD = pairs.getLibrarySD(iLib);
	  if ( abs(shift) < libSD / 3) continue;
	  for(size_t k=0; k < pair_ids[iLib].size(); k++) 
	  {
	    longlong pairId = pair_ids[iLib][k];
	    seps[pairId] += shift;
	    counter2++;
	  }
	  if (verbose) 
	  cout << "c"<< cg1 << "," << cg2 << ( fw1 ? '+':'-') <<  (fw2 ? '+' : '-')
	    << " lib "<< iLib << ", " << nlinks << " links, shift= " << shift << endl;
	}
      }
    }
  }
  if (verbose) cout << "Total pair separation adjustments= " << counter2 << endl;
  cout << "Total pair separation adjustments= " << counter2 << endl;
  cout << Date() << ": End converting " << endl;
}


// perform pair distance correction using the distribution obtained from IntDist.
// Note that the invariant separations are used in IntDist, while current MakeScaffold
// codes requires end-to-end separation. Therefore, read lengths are needed for the 
// conversion. 
//
void PairDistCorrectionFromIntDistOld(const String reads_head,  const PairsManager& pairs, const vec<alignlet>& aligns, const vec<int>& index, 
    vec<int>&seps, vec<int>& sds, bool verbose) ;

void PairDistCorrectionFromIntDist(const String reads_head,  const PairsManager& pairs, const vec<alignlet>& aligns, const vec<int>& index, 
    vec<int>&seps, vec<int>& sds, bool verbose) ;

// perform pair distance correction using the given distribution 
void PairDistCorrectionFromFile(const vec<String>& pairs_dist_files, const PairsManager& pairs, const vec<alignlet>& aligns, const vec<int>& index, 
    vec<int>&seps, vec<int>& sds, bool verbose) ;



#endif
