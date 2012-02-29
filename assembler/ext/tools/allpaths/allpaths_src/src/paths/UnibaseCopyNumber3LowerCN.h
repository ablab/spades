///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// UnibaseCopyNumber3Core: predict copy numbers.  Warning: this creates files
// in a way which it shouldn't.  This feature will presumably be eliminated in
// a subsequent version.

#ifndef UNIBASE_COPY_NUMBER3_LOWER_CN
#define UNIBASE_COPY_NUMBER3_LOWER_CN

#include "Basevector.h"
#include "CoreTools.h"

#include "paths/RemodelGapTools.h"
#include "paths/Ulink.h"
#include "paths/PairDistFitting.h"
#include "util/SearchFastb2Core.h"
#include "paths/simulation/Placement.h"
#include "paths/UnipathFixerTools.h"
#include "math/Functions.h"
#include "ParallelVecUtilities.h"
#include "PairsManager.h"


Bool cmp_ridx( const segalign& a1, const segalign& a2 ){    
  return a1.rid < a2.rid;    
}

// !! requires:   ParallelSort(JSEGS, cmp_ridx);
void MakeUlinksX( const String& run_dir, const String& JUMP_READS, const int K, 
     const vec<int>& to_rc, const vec<segalign>& SUGS, const vecbasevector& unibases, 
     const vec<int>& innie_sep, const vec<int>& innie_dev, 
     const vec<double>& innie_percent, const double min_innie_percent,
     const int min_kmers, const int max_devs, vec<ulink_with_uids>& ulinks );

void GapStatsAlt( vec<int> gap, vec<int> gapdev, int& gap_ave, int& gapdev_ave );


// Wrapper of the original code in main() that find all links between unibases
// and determine the most probable gap size and standard deviation between them
void GetUnibaseSeps( 
    // input
    const String & run_dir,
    const String & JUMP_READS,
    const int K,
    const vec<int> & to_rc,
    const vec<segalign> & JSEGS,
    const vecbasevector & unibases,
    const vec<int> & innie_sep,
    const vec<int> & innie_dev,
    const vec<double> & innie_percent,
    const vec<double> & CN_raw,
    // output
    vec<ulink_with_uids> & condensed_links,
    vec<double> & cn_from,
    vec<double> & cn_to,
    vec< vec< triple<ho_interval,int,int> > > & cov_from, 
    vec< vec< triple<ho_interval,int,int> > > & cov_to,
    vec<int> & most_links_from,
    vec<int> & most_links_to
    );

// Define upper cutoff.  To do this, walk backwards from the end of the 
// distribution, sampling 100 points at a time.  Taking into account
// the visibility, estimate the noise rate and error bars for it, assuming
// that the 100 points are all noise.  Now keep walking until we think
// the signal:noise ratio is at least 10.  This defines the upper cutff.
void FindDistributionCutoff( 
          // input
          const map<int,int>& hist, 
          const vec<int64_t>&  VISIBLE, 
          // output
          int & cutoff 
          ) ;


void GetDistribution( 
          //input
          const map<int,int>& hist, const vec<int>& lens_input, 
          //output
          vec<double>& distr, 
          int& start,           // the distribution density at (start) is distr[0] 
          // parameters
          bool SMOOTH = true
          );

void UnibaseGap( 
          // input
          const vec< ProbFuncIntDist >& distrs, 
          const vec< vec< pair<int,int> > >& links,
          int len1, int len2,  
          // output
          int & gap, int & std 
    );

// Experimental code that find all links between unibases
// and determine the most probable gap size and standard deviation between them
void GetUnibaseSepsExp( 
    // input
    const String & run_dir,
    const String & JUMP_READS,
    const int K,
    const vec<int> & to_rc,
    const vec<segalign> & SUGS,
    const vecbasevector & unibases,
    const vec<int> & innie_sep,
    const vec<int> & innie_dev,
    const vec<double> & innie_percent,
    const vec<double> & CN_raw,
    // output
    vec<ulink_with_uids> & condensed_links,
    vec<double> & cn_from,
    vec<double> & cn_to,
    vec< vec< triple<ho_interval,int,int> > > & cov_from, 
    vec< vec< triple<ho_interval,int,int> > > & cov_to,
    vec<int> & most_links_from,
    vec<int> & most_links_to
    );


// Wrapper to call functions in RemodelGapsTools to align reads to unibases
// and estimate the gaps
void UnibaseGapRemodeling ( 
          // input
          const String& run_dir,
          const String& JUMP_READS,
          const vecbasevector& unibases,
          const int K,
          const vec<int> & to_rc,
          const vec<double> & CN_raw,
          // output
          vec<ulink_with_uids> & condensed_links,
          vec<double> & cn_from,
          vec<double> & cn_to,
          vec< vec< triple<ho_interval,int,int> > > & cov_from, 
          vec< vec< triple<ho_interval,int,int> > > & cov_to,
          vec<int> & most_links_from,
          vec<int> & most_links_to,
          // parameters
          bool REGAP,  // only re-evaluate the gaps that are already in condensed_links
          int VERBOSITY
          );


#endif
