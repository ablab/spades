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

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS
#include <omp.h>

#include "Basevector.h"
#include "CoreTools.h"

#include "ParallelVecUtilities.h"
#include "paths/UnipathFixerTools.h"
#include "paths/Ulink.h"
#include "math/HoInterval.h"
#include "paths/PairDistFitting.h"
#include "paths/UnibaseCopyNumber3LowerCN.h"
#include "paths/RemodelGapTools.h"

#include "paths/LinkingPairs.h"
#include "math/Array.h"

void MakeUlinksX( const String& run_dir, const String& JUMP_READS, const int K, 
     const vec<int>& to_rc, const vec<segalign>& SUGS, const vecbasevector& unibases, 
     const vec<int>& innie_sep, const vec<int>& innie_dev, 
     const vec<double>& innie_percent, const double min_innie_percent,
     const int min_kmers, const int max_devs, vec<ulink_with_uids>& ulinks )
{    cout << Date( ) << ": set up for link computation" << endl;
     String head = run_dir + "/" + JUMP_READS;
     vecbasevector reads( head + ".fastb" );
     uint64_t nreads = reads.size( );
     vec<int> readlen( reads.size( ) );
     for ( size_t i = 0; i < reads.size( ); i++ )
          readlen[i] = reads[i].size( );
     Destroy(reads);
     PairsManager pairs( head + ".pairs" );
     vec<size_t> S_START_RID(nreads+1);
     {    size_t SEG_POS = 0;
          for ( uint64_t rid = 0; rid <= nreads; rid++ )
          {   while( SEG_POS < SUGS.size( ) && SUGS[SEG_POS].rid < rid ) ++SEG_POS;
              S_START_RID[rid] = SEG_POS;    }    }
     cout << Date( ) << ": compute links from "
         << ToStringAddCommas( pairs.nPairs( ) ) << " pairs" << endl;
     #pragma omp parallel for
     for ( size_t i = 0; i < pairs.nPairs( ); i++ ){    
          longlong id1 = pairs.ID1(i), id2 = pairs.ID2(i);
    
     // Require that jump reads have only one placement.  This may be
     // too stringent.
    
     if ( S_START_RID[id1+1] - S_START_RID[id1] > 1 ) continue;
     if ( S_START_RID[id2+1] - S_START_RID[id2] > 1 ) continue;
    
     for (size_t j1 = S_START_RID[id1]; j1 < S_START_RID[id1+1]; j1++)
     {    for (size_t j2 = S_START_RID[id2]; j2 < S_START_RID[id2+1]; j2++)
          {    int u1 = SUGS[j1].u, u2 = to_rc[ SUGS[j2].u ];
	       if ( u1 == u2 || u1 == to_rc[u2] ) continue;
	       if ( unibases[u1].isize( ) - K + 1 < min_kmers ) continue;
	       if ( unibases[u2].isize( ) - K + 1 < min_kmers ) continue;
	       int nu1 = unibases[u1].size( ), nu2 = unibases[u2].size( );
	       int nr1 = readlen[id1], nr2 = readlen[id2];
	       int upos1 = SUGS[j1].upos, upos2 = SUGS[j2].upos;
	       int lib = pairs.libraryID(i);
	       // Library -> separation/deviation
	       // No need to add the read length since the reads are already flipped before alignment
	       int sep_in = innie_sep[lib] //+ nr1 + nr2
	            - ( nu1 - upos1 + SUGS[j1].rpos + nu2 
	            - upos2 + SUGS[j2].rpos );
	       int dev_in = innie_dev[lib];
	       int sep_out = pairs.getLibrarySep(pairs.libraryID(i))
	            - ( upos1 - SUGS[j1].rpos
	            + upos2 - SUGS[j2].rpos ) ; //+ nr1 + nr2;
	       int dev_out = pairs.getLibrarySD(pairs.libraryID(i));
	
	       // Add this link to the ulinks list - in both fw
	       // and rc form, as both an innie and an outie.
	
               #pragma omp critical
	       {    if ( innie_percent[lib] >= min_innie_percent )
                    {    ulinks.push( u1, u2, sep_in, dev_in, 
			      Min( nu1, upos1 + nr1 ), 
			      Max( 0, nu2 - upos2 - nr2 ), 1 );
	                 ulinks.push( to_rc[u2], to_rc[u1], sep_in, 
			      dev_in, Min( nu2, upos2 + nr2 ), 
			      Max( 0, nu1 - upos1 - nr1 ), 1 );    }
	            ulinks.push( u2, u1, sep_out, dev_out,
		         Min( nu2, nu2 - upos2 ), Max(0, upos1), 1 );
	            ulinks.push( to_rc[u1], to_rc[u2], sep_out, 
		         dev_out, Min( nu1, nu1 - upos1 ),
		         Max(0, upos2), 1 );    }    }    }    }
     cout << Date( ) << ": sorting " << ToStringAddCommas( ulinks.size( ) ) 
          << " unipath links" << endl;
     ParallelSort(ulinks);    }

void GapStatsAlt( vec<int> gap, vec<int> gapdev, int& gap_ave, int& gapdev_ave ){
     // If there are less than six gaps, we directly compute their mean.
     // Otherwise, we attempt to remove outliers, as follows.  We sort 
     // the gaps and extract the middle half.  From this middle half, we compute the 
     // mean and standard deviation.  Then we select those gaps lying withing 5 
     // standard deviations of the mean, and form their mean.  

     vec<normal_distribution> S;
     if ( gap.size( ) >= 6 ){    
       vec<int> mid_gaps;
       SortSync( gap, gapdev );
       for ( unsigned int i = gap.size( )/4; i < 3*(1+gap.size( ))/4; i++ )
	 mid_gaps.push_back( gap[i] );
       float sum1 = 0, sum2 = 0;
       for ( unsigned int i = 0; i < mid_gaps.size( ); i++ ){    
	 sum1 += mid_gaps[i];
	 sum2 += float(mid_gaps[i]) * float(mid_gaps[i]);    
       }
       float n = mid_gaps.size( );
       float mean = sum1/n;
       float sd = sqrt(sum2/n - mean * mean);
       float start = mean - 5 * sd, stop = mean + 5 * sd;
       for ( unsigned int i = 0; i < gap.size( ); i++ ){    
	 if ( start <= gap[i] && gap[i] <= stop )
	   S.push( gap[i], gapdev[i] );    
       }    
     }
     else{    
       for ( int l = 0; l < gap.isize( ); l++ )
	 S.push( gap[l], gapdev[l] );    
     }
     normal_distribution s = CombineNormalDistributions(S);
     gap_ave = int(round(s.mu_));
     gapdev_ave = int(round(s.sigma_));    
}


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
    )
{
      // Heuristics.
  
      const double max_devs = 6.0;
      const int min_kmers = 120;
      const int min_links_initial = 2;
      const double min_innie_percent = 5.0;
  
      // For each pair of unipaths that are connected by two or more links, predict
      // their order and separation.  Both orders may be possible.  Note that a more 
      // wholistic approach may be needed.
  
      vec<ulink_with_uids> ulinks;
      MakeUlinksX( run_dir, JUMP_READS, K, to_rc, JSEGS, unibases, innie_sep, 
           innie_dev, innie_percent, min_innie_percent, min_kmers, 
           max_devs, ulinks );
  
      // Make the list of putative unipath links into a set of condensed links.
      // Here we use GapStats to merge links, and we apply a min_links_initial 
      // threshold.  A condensed links is accepted if its predicted overlap is no 
      // more than K-1, plus slop.
      
      cout << Date( ) << ": condense links" << endl;
      int dev_mult = 3;
      for ( size_t i = 0; i < ulinks.size( ); i++ )
      {    size_t j;
           int u1 = ulinks[i].u1, u2 = ulinks[i].u2;
           for ( j = i + 1; j < ulinks.size( ); j++ )
           if ( ulinks[j].u1 != u1 || ulinks[j].u2 != u2 ) break;
           int nlinks = j - i;
           if ( nlinks >= min_links_initial )
           {    vec<int> seps, devs;
                int sep, dev, start1 = unibases[u1].size( ), stop2 = 0;
                for ( size_t m = i; m < j; m++ )
                {    const ulink_with_uids& f = ulinks[m];
	             seps.push_back(f.sep), devs.push_back(f.dev);    
	             start1 = Min( start1, f.start1 );
	             stop2 = Max( stop2, f.stop2 );    }
                GapStatsAlt( seps, devs, sep, dev );
		// Because of the biased distribution of pairs over the gap,
		// the gap size calculated from simple averaging can sometimes be 
		// off by as much as the distribution width, even if there are
		// hundreds of links between them. We need to keep such gaps for
		// more accurate remodeling in the following stages.
                if ( sep + max_devs * dev >= -(K-1) ||
                     sep + Mean(devs) * 2.5 >= -(K-1) && nlinks >= 5 )
                {    condensed_links.push( u1, u2, sep, dev, start1, stop2, j-i );
	             int Start2 = sep + dev_mult * dev;
	             int Stop2 = sep + unibases[u2].isize( ) - dev_mult * dev;
	             if ( Start2 < Stop2 ) 
                     {    ho_interval h2( Start2, Stop2 );
	                  cov_from[u1].push( h2, nlinks, u2 );    }
	             int Start1 = -sep - unibases[u1].isize( ) + dev_mult * dev;
	             int Stop1 = -sep - dev_mult * dev;
	             if ( Start1 < Stop1 ) 
                     {    ho_interval h1( Start1, Stop1 );
	                  cov_to[u2].push( h1, nlinks, u1 );    }
	             if ( nlinks > most_links_from[u1] )
                     {    most_links_from[u1] = j-i;
	                  cn_from[u1] = CN_raw[u2];    }
	             if ( nlinks > most_links_to[u2] )
                     {    most_links_to[u2] = j-i;
	                  cn_to[u2] = CN_raw[u1];    }    }    }
           i = j - 1;    }

}

// Define upper cutoff.  To do this, walk backwards from the end of the 
// distribution, sampling 100 points at a time.  Taking into account
// the visibility, estimate the noise rate and error bars for it, assuming
// that the 100 points are all noise.  Now keep walking until we think
// the signal:noise ratio is at least 10.  This defines the upper cutff.
void FindDistributionCutoff( const map<int,int>& hist, 
     const vec<int64_t>&  VISIBLE, // input
     int & cutoff // output
          ) 
{
     cutoff = hist.rbegin()->first + 1;  // default guess

     const int bin_count = 100;
     const double err_mult = 2.0;
     const double max_noise_frac = 0.1;
     // hist_sum[x] is number of links whos separation is less or equal to x
     map<int,int>  hist_sum;
     int sum = 0;
     for( map<int,int>::const_iterator it = hist.begin(); it != hist.end(); ++it ) {
          sum += it->second;
          hist_sum[ sum ] = it->first;
     }
     cout << "sum= " << sum<< endl;
     if ( sum < bin_count )  return;
     // start
     int N = sum;
     double min_freq_high = 1000000000;
     double mean_visable = Mean(VISIBLE);
     for ( int rx = 0; N - (rx+1)*bin_count >= 0; rx++ ) {    
          double freq = 0, effective_size = 0;
          int high = hist_sum.lower_bound( N - (rx)*bin_count - 1 )->second;
          for ( int j = N - (rx+1)*bin_count; j < N - (rx)*bin_count; j++ ) {    
               int p = hist_sum.lower_bound( j )->second, q = high;
               freq += mean_visable / double( VISIBLE[abs(p)] );
               effective_size += double(VISIBLE[abs(q)]) / double( VISIBLE[abs(p)] );    
          }
          int low = hist_sum.lower_bound( N - (rx+1)*bin_count )->second;
          if ( high == low ) {
               cutoff = high;
               break;
          }
          freq /= double( high - low );
          double freq_err = freq / sqrt(effective_size);
          double freq_low = Max( 0.0, freq - err_mult * freq_err );
          double freq_high = freq + err_mult * freq_err;
          min_freq_high = Min( min_freq_high, freq_high );
          double delta = freq_low / min_freq_high;
          //PRINT5( low, high, freq_low, freq_high, delta );
          if ( 1.0 / delta <= max_noise_frac ) {    
               cutoff  = high;
               break;    
          }    
     }    
}


void GetDistribution( const map<int,int>& hist, const vec<int>& lens_input, //input
          vec<double>& distr, 
          int& start, // output
          bool SMOOTH // parameters
          )
{
     // Define visibility at given distances.
     // Both positive and negative sides are considered
     int max_tig = Max ( abs(hist.rbegin()->first) + 1, abs(hist.begin()->first) + 1);
     vec<int64_t> VISIBLE( max_tig, 0 );
     {
          vec<int> lens( lens_input );
          Sort( lens );
          int ntigs = lens.size();
          vec<int> sum_len( ntigs, 0 ); // sum_len[j] = sum_{i, where i >= j} ( len(i) )
          partial_sum( lens.rbegin(), lens.rend(), sum_len.rbegin() ); 
          int j = 0;
          for ( int i = 0; i < max_tig; i++ ) {
               // find the contig j, so that lens[j] > i
               if ( j >= ntigs ) break;
               if ( i > lens[j] ) {
                    for ( ; j < ntigs; j++ ) 
                         if ( lens[j] > i ) break;
                    if ( j >= ntigs ) break;
               }
               VISIBLE[i] = sum_len[j] - ( ntigs - j ) * i;
          }
     }
     // upper bound of the distribution
     int cutoff_ub;
     FindDistributionCutoff(hist, VISIBLE, cutoff_ub );
     cout << "  cutoff_ub= " << cutoff_ub << endl;
     // lower bound of the distribution
     int cutoff_lb;
     {
          map <int, int> hist_neg; // negative side
          map<int,int>::const_reverse_iterator it = hist.rbegin();
          for ( ; it != hist.rend(); it++ )
               hist_neg[ -it->first ] = it->second;
          int ub;
          FindDistributionCutoff(hist_neg, VISIBLE, ub);
          cutoff_lb = - ub + 1;
     }
     cout << "  cutoff_lb= " << cutoff_lb << endl;
     // populate the vector , make corrections, and re-normalize
     // distr is an array holds the distribution density in [ cutoff_lb, cutoff_ub )
     start = cutoff_lb;
     distr.assign( cutoff_ub - cutoff_lb, 0);
     for( map<int,int>::const_iterator it = hist.begin(); it != hist.end(); ++it ) {
         if ( it->first < cutoff_lb) continue;
         if ( it->first >= cutoff_ub) break;
         distr[it->first - cutoff_lb] += it->second;
     }
     for ( int i = 0; i < distr.isize(); i++ ) {
          int x = abs( i + cutoff_lb );
          ForceAssertGe( VISIBLE[x], 0 );
          distr[i] /= double(VISIBLE[x]);
     }
     // smoothing
     if ( SMOOTH ) { SmoothArrayGaussian( distr, 50 ); }

     double total = Sum( distr );
     for ( int i = 0; i < distr.isize(); i++ ) 
          distr[i] /= total;
}

void UnibaseGap( const vec< ProbFuncIntDist >& distrs, const vec< vec< pair<int,int> > >& links,
         int len1, int len2,  // input
    int & gap, int & std // output
    )
{
     int nlibs2 = distrs.size();
     ForceAssertEq( nlibs2, links.isize() );
     vec< normal_distribution > combine;
     for ( int libid = 0; libid < nlibs2; libid++ ) {
          int gap_i = -1, std_i = -1;
          if ( ! links[libid].empty() )
               MostProbableGap( distrs[libid], len1, len2, links[libid], gap_i, std_i);
          if ( std_i > 0 )
               combine.push( gap_i, std_i );
     }
     if ( combine.empty() ) return;
     normal_distribution nd = CombineNormalDistributions( combine);
     gap = nd.mu_;
     std = nd.sigma_;
}


void UnibaseGap2( const vec< ProbFuncIntDist >& distrs, 
         int len1, int len2,  
         const vec< vec< pair<int,int> > >& links12,
         const LibCov_t & starts1,
         const LibCov_t & stops2,
         const vec< vec< pair<int,int> > >& links21,
         const LibCov_t & starts2,
         const LibCov_t & stops1,
         // end of input
    int & gap, int & std // output
    )
{
     int nlibs = distrs.size();
     ForceAssertEq( nlibs, links12.isize() );
     ForceAssertEq( nlibs, links21.isize() );
     ForceAssertEq( nlibs, starts1.isize() );
     ForceAssertEq( nlibs, starts2.isize() );
     ForceAssertEq( nlibs, stops1.isize() );
     ForceAssertEq( nlibs, stops2.isize() );
}


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
    )
{
     String head = run_dir + "/" + JUMP_READS;
     PairsManager pairs( head + ".pairs" );
     int nlibs = pairs.nLibraries();
     ForceAssertEq( nlibs, innie_sep.isize() );
     ForceAssertEq( nlibs, innie_dev.isize() );
     ForceAssertEq( nlibs, innie_percent.isize() );

     // Heuristics.
     const double max_devs = 6.0;
     const int min_kmers = 120;
     const int min_links_initial = 2;
     const double min_innie_percent = 5.0;

     // For each pair of unipaths that are connected by two or more links, predict
     // their order and separation.  Both orders may be possible.  Note that a more 
     // wholistic approach may be needed.
     cout << Date( ) << ": set up for link computation" << endl;
     vecbasevector reads( head + ".fastb" );
     uint64_t nreads = reads.size( );
     vec<size_t> S_START_RID(nreads+1);
     {    size_t SEG_POS = 0;
          for ( uint64_t rid = 0; rid <= nreads; rid++ )
          {   while( SEG_POS < SUGS.size( ) && SUGS[SEG_POS].rid < rid ) ++SEG_POS;
              S_START_RID[rid] = SEG_POS;    }    }
     cout << Date( ) << ": compute links from "
         << ToStringAddCommas( pairs.nPairs( ) ) << " pairs" << endl;
     
     // Only pick unibases that are long enough
     vec< bool > unibase_pick( unibases.size(), false ); 
     for ( size_t u = 0; u < unibases.size(); u++ )
          if ( unibases[u].isize( ) - K + 1 >= min_kmers )
               unibase_pick[u] = true;

     // prepare the alignments for each read. must be uniquely aligned
     vec< pair<int,int> > aligns( nreads, make_pair(-1,-1) );
     for ( size_t i = 0; i < pairs.nPairs( ); i++ ) {    
          longlong id1 = pairs.ID1(i), id2 = pairs.ID2(i);
          int libid = pairs.libraryID(i);
          // Force uniq alignment !
          if ( S_START_RID[id1+1] - S_START_RID[id1] == 1 ) {
               size_t j1 = S_START_RID[id1];
               int x1 = SUGS[j1].upos - SUGS[j1].rpos; 
               int u1 = SUGS[j1].u;
	       int nu1 = unibases[u1].size( );
	       if ( unibase_pick[u1] && x1 >= 0 && x1 < nu1 )
                         aligns[ id1 ] = make_pair( u1, x1 );
          }
          if ( S_START_RID[id2+1] - S_START_RID[id2] == 1 ) {
               size_t j2 = S_START_RID[id2];
               int x2 = SUGS[j2].upos - SUGS[j2].rpos; 
               int u2 = SUGS[j2].u;
	       int nu2 = unibases[u2].size( );
	       if ( unibase_pick[u2] && x2 >= 0 && x2 < nu2 )
                    aligns[ id2 ] = make_pair( u2, x2 );
          }
     }
     // Gather the linking infomation
     vec< int > unibase_lens( unibases.size(), -1); 
     for ( size_t u = 0; u < unibases.size(); u++ )
          unibase_lens[u] = unibases[u].size();
     LinkingPairs linking( nlibs, unibase_lens );
     GatherLinks( pairs, aligns, to_rc, linking );

     // The distribution models
     vec< ProbFuncIntDist > distrs( nlibs);
     for ( int libid = 0; libid < nlibs; libid++ ) {
          const map<int,int>& dist = linking.GetDists(libid);
         // cout << "histx" << libid << " ";
         // for( map<int,int>::const_iterator it = dist.begin(); it != dist.end(); it++ )
         //      cout << it->first << " ";
         // cout << endl;
         // cout << "histy" << libid << " ";
         // for( map<int,int>::const_iterator it = dist.begin(); it != dist.end(); it++ )
         //      cout << it->second << " ";
         // cout << endl;
          vec<int> lens;
          for ( int i = 0; i < (int) unibases.size(); i++ )
               if ( unibase_pick[i] ) lens.push_back( unibases[i].size() );
          Sort(lens);
          vec<double> distr_array;
          int cutoff_lb ;
          GetDistribution( dist, lens, distr_array, cutoff_lb,  true);
          //cout << scientific; 
          //for ( int i = 0; i < distr_array.isize(); i++ ) {
          //     if ( i % 10 != 0 ) continue;
          //     int x = i + cutoff_lb;
          //     cout << "distr" << libid << " " << x << " " << distr_array[i] << endl;
          //}
	  distrs[libid].FromArray( distr_array, cutoff_lb );
     }
     cout << Date( ) << ": condense links" << endl;
     int dev_mult = 3;
     
     vec< pair<int,int> > all_linked = linking.GetAllLinked( );
     for ( int i = 0; i < all_linked.isize(); i++ ) {
          // We get all the linked pairs where there are read pairs that start from (rc align)
          // u1 and stop at (fw align) u2. Two possible arrangements :
          // (1)
          //                  -----------------u1-------------- (gap) ------u2------
          //                                         <---@ r1             @-->r2
          // (2)
          // ----u2---- (gap) -----------------u1--------------
          //  @-->r2              <---@ r1           
          // Case (2) can be treated as case (1) with a large negative gap.
          //
          //
          // We should also check pairs which start from u2 and end at u1.
          // (1)
          //                  -----------------u1-------------- (gap) ------u2------
          //                                         @---> r1             <--@r2
          // (2)
          // ----u2---- (gap) -----------------u1--------------
          //  <--@r2              @---> r1           
          int u1 = all_linked[i].first;
          int u2 = all_linked[i].second;
          const LibLink_t & ulinks12= linking.GetLinks( u1, u2 ); // starts from u1, ends at u2
          const LibLink_t & ulinks21= linking.GetLinks( u2, u1 ); // starts from u2, ends at u1
          int nlinks = 0;
          for ( int libid = 0; libid < nlibs; libid++ ) 
               nlinks += ulinks12[libid].size() + ulinks21[libid].size();
          if ( nlinks < min_links_initial ) continue;

          int sep=-1, dev=-1;
          // Below we still obtain the sep and dev without using distribution
          // new module can be inserted here considering pair distribution corrections
          {  
               vec <int> seps,  devs;
               for ( int libid = 0; libid < nlibs; libid++ ) {
                    for ( int k = 0; k < ulinks12[libid].isize(); k++ ) {
                         seps.push_back( pairs.getLibrarySep(libid) - ulinks12[libid][k].first - 
                                   ulinks12[libid][k].second ); 
                         devs.push_back(  pairs.getLibrarySD( libid )); }
               }
               GapStatsAlt( seps, devs, sep, dev );
          }

          int dev_temp = -1;
          //UnibaseGap( distrs, ulinks, unibases[u1].size(), unibases[u2].size(), sep, dev_temp);
          UnibaseGap2( distrs, unibases[u1].size(), unibases[u2].size(), 
                    ulinks12, linking.GetStarts(u1), linking.GetStops(u2), 
                    ulinks21, linking.GetStarts(u2), linking.GetStops(u1), // end of input
                    sep, dev_temp);

          if ( sep + max_devs * dev <  -(K-1) ) continue;
          int start1 = 0, stop2 = 0; // those numbers are actually not used afterwards


          condensed_links.push( u1, u2, sep, dev, start1, stop2, nlinks );
          int Start2 = sep + dev_mult * dev;
          int Stop2 = sep + unibases[u2].isize( ) - dev_mult * dev;
          if ( Start2 < Stop2 ) 
          {    ho_interval h2( Start2, Stop2 );
               cov_from[u1].push( h2, nlinks, u2 );    }
          int Start1 = -sep - unibases[u1].isize( ) + dev_mult * dev;
          int Stop1 = -sep - dev_mult * dev;
          if ( Start1 < Stop1 ) 
          {    ho_interval h1( Start1, Stop1 );
               cov_to[u2].push( h1, nlinks, u1 );  }
          if ( nlinks > most_links_from[u1] )
          {    most_links_from[u1] = nlinks;
               cn_from[u1] = CN_raw[u2];    }
          if ( nlinks > most_links_to[u2] )
          {    most_links_to[u2] = nlinks;
               cn_to[u2] = CN_raw[u1];    }    
     }
}

void FlipAlign( int m, 
          const vecbasevector& tigs,
          const PairsManager& pairs,
          const vec< unsigned short int >& base_lens,
          const vec< vec<longlong> > & aligns_index,
          vec<Bool> & placed_fw, 
          vec<Bool> & placed_rc,
          vec< pair<int,int> > & placement
     )
{
     const vec<longlong>& pids_on_tig = aligns_index[m];
     for ( size_t i = 0; i < pids_on_tig.size(); i++ ) {
          longlong pid = pids_on_tig[i];
          longlong rid1 = pairs.ID1(pid);
          longlong rid2 = pairs.ID2(pid);
          if ( ( placed_fw[rid1] || placed_rc[rid1] ) && placement[rid1].first == m ) {
               swap( placed_fw[rid1], placed_rc[rid1] );
               placement[rid1].second = tigs[m].size() - base_lens[rid1] - placement[rid1].second;
          }
          if ( ( placed_fw[rid2] || placed_rc[rid2] ) && placement[rid2].first  == m ) {
               swap( placed_fw[rid2], placed_rc[rid2] );
               placement[rid2].second = tigs[m].size() - base_lens[rid2] - placement[rid2].second;
          }
     }
}

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
          )
{
     // Heuristics.
     int dev_mult = 3;
     const double max_devs = 6.0;
     const int min_kmers = 120;
     const int min_links_initial = 2;   // need at least 5 links for gap estimation
     
     // Save the existing gap data
     map< pair<int,int>, ulink_with_uids> input_gaps;
     if ( REGAP ) {
          for ( int i = 0; i < condensed_links.isize(); i++ )
               input_gaps[ make_pair( condensed_links[i].u1, condensed_links[i].u2 ) ]
                    = condensed_links[i];
     }
     int nuni = unibases.size();

     // get rid of the rc of each unibases
     // and creat mapping

     vecbasevector unibases_norc;   // unibases_norc[m] is a non-redundant set of unibases
     vec< int > m2u;        
     vec< int > u2m(unibases.size(), -1 );        
     for ( int i = 0; i < (int) unibases.size(); i++ ){
          if ( to_rc[i] < i ) continue;
          if ( unibases[i].isize() - K + 1 < min_kmers  ) continue;
          unibases_norc.push_back( unibases[i] );
          m2u.push_back( i );
          u2m[i] = unibases_norc.size() - 1;
     }

     // Build a map of kmers in the contigs.

     cout << Date( ) << ": building kmer maps for all unibases" << endl;
     const int KK = 40;
     vec< triple<kmer<KK>,int,int> > kmers_plus;
     MakeKmerLookup( unibases, kmers_plus );
     cout << Date( ) << ": copying" << endl;
     vec< kmer<KK> > kmers( kmers_plus.size( ) );
     for ( size_t i = 0; i < kmers.size( ); i++ )
          kmers[i] = kmers_plus[i].first;

     // Count number of libraries.

     vec<int> libs_to_use;
     int nlibs = 0;
     GetLibCount( run_dir, JUMP_READS, nlibs, VERBOSITY, True );

     if ( VERBOSITY >= 1)
          cout << "nlibs= " << nlibs << endl;

     // Align reads.

     cout << Date( ) << ": loading reads" << endl;
     vecbasevector bases( run_dir + "/" + JUMP_READS + ".fastb" );
     PairsManager pairs( run_dir +  "/" + JUMP_READS + ".pairs" );
     vec<Bool> placed_fw, placed_rc;
     vec< pair<int,int> > placement;
     vec< vec<longlong> > aligns_index;
     AlignReads( 40, 20, true, unibases_norc, bases, pairs, placed_fw, placed_rc, placement, aligns_index, True);

     // Define distributions for libraries.

     vec<Bool> lfail;
     vec<int> max_dist;
     vec< vec<int> > DISTS;
     vec< vec<double> > X;
     vec< unsigned short int > base_lens( bases.size() );
     for ( int i = 0; i < base_lens.isize(); i++ ) base_lens[i] = bases[i].size();
     const Bool HISTOGRAMS = False;
     DefineDistributions( nlibs, 0, unibases_norc, base_lens, pairs, placed_fw,
          placed_rc, placement, aligns_index, True, HISTOGRAMS, lfail, 
          max_dist, DISTS, X );

     map< pair<int,int>, int > nlinks_unibases;
     // Which pairs of unibases are possibly connected?
     for ( size_t pid = 0; pid < pairs.nPairs( ); pid++ ) {    
          int64_t id1 = pairs.ID1(pid), id2 = pairs.ID2(pid);
          if ( ! (placed_fw[id1] || placed_rc[id1]) ) continue;
          if ( ! (placed_fw[id2] || placed_rc[id2]) ) continue;
          int m1 = placement[id1].first;
          int m2 = placement[id2].first;
          if ( m1 == m2 ) continue;
          // What unibases are the reads forward aligned?
          int u1 = placed_fw[id1] ? m2u[m1] : to_rc[ m2u[m1] ] ;
          int u2 = placed_fw[id2] ? m2u[m2] : to_rc[ m2u[m2] ] ;
          // possible gaps
          nlinks_unibases[ make_pair(u1, to_rc[u2]) ]++;
          nlinks_unibases[ make_pair(u2, to_rc[u1]) ]++;
     }

     // The gaps we want to remodel
     vec< pair<int,int> > all_gaps;
     if ( REGAP ) 
          for( map< pair<int,int>, ulink_with_uids> :: iterator it = input_gaps.begin();
                    it != input_gaps.end(); it++ )
               all_gaps.push_back( it->first );
     else {
          for (map< pair<int,int>, int >:: iterator it = nlinks_unibases.begin(); 
                    it != nlinks_unibases.end(); it++ ) {
               int nlinks = it->second;
               if ( nlinks >= min_links_initial ) all_gaps.push_back( it->first );
          }
     }

     // Becase gap from u1 to u2 is the same as gap from u2_rc to u1_rc,
     // potentially we need to check only one of them.
     // Also potentially the gap remodel function may give slightly different
     // values for gap( u1, u2 ) and gap( u2_rc, u1_rc ), resulting
     // in different copy numbers for u1 and u1_rc, which should be prevented.
     map< pair<int,int>, bool > mark_pick;
     vec <Bool> to_erase( all_gaps.size(), False );
     for ( int i = 0; i < all_gaps.isize(); i++ ) {
          int u1 = all_gaps[i].first;
          int u2 = all_gaps[i].second;
          if ( u2m[u1] >= 0 && u2m[u2] >=0 ) mark_pick[ all_gaps[i] ] = true;
          else if ( u2m[u1] < 0 && u2m[u2] < 0 ) to_erase[i] = true;
          else if ( mark_pick.find( make_pair( to_rc[u2], to_rc[u1] ) )
                    == mark_pick.end() )
               mark_pick[ all_gaps[i] ] = true;
          else
               to_erase[i] = true;
     }
     EraseIf( all_gaps, to_erase );

     // Parallization of the gap estimation.
     // We need to pay special attention because we some times need to flip the unibases 
     // to test alternative unibase linking orientation, we need change the corresponding
     // alignment to use the gap-remodel modules which assume only forward orientation.
     // Because we may need to flip the alignments of a pair of unibases, we need to 
     // make sure other threads will work on other unibases.
     
     vec <ulink_with_uids> condensed_links_new;
     cout << Date( ) << ": Estimating sizes for " << all_gaps.size() << " gaps."  << endl;
     // We first separate the gaps in different groups so that any two gaps in the same 
     // group don't share unibases.
     while ( all_gaps.size() > 0 ) {
          vec< pair<int,int> >  gaps_in_group;
          set<int>  tigs_in_group;
          vec< Bool > to_del( all_gaps.size(), false);
          for ( int i = 0; i < all_gaps.isize(); i++ ) {
               int u1 = all_gaps[i].first;
               int u2 = all_gaps[i].second;
               int m1 = ( u2m[u1] >= 0 ) ? u2m[u1] : u2m[ to_rc[u1] ];
               int m2 = ( u2m[u2] >= 0 ) ? u2m[u2] : u2m[ to_rc[u2] ];
               if ( tigs_in_group.find(m1) != tigs_in_group.end()
                 || tigs_in_group.find(m2) != tigs_in_group.end() ) 
                         continue;
               gaps_in_group.push_back( all_gaps[i] );
               tigs_in_group.insert( m1 );
               tigs_in_group.insert( m2 );
               to_del[i] = true;
          }
          EraseIf(all_gaps, to_del);
          if ( VERBOSITY >= 4) {
               cout << "group : ";
               for ( int j = 0; j < gaps_in_group.isize(); j++ ) 
                    cout <<  gaps_in_group[j].first << "<->" << gaps_in_group[j].second << " ";
               cout << endl;
               cout << "gaps left : " << all_gaps.size() << endl;
          }

     #pragma omp parallel for
     for ( int j = 0; j < gaps_in_group.isize(); j++ ) {
          pair<int,int> upair = gaps_in_group[j];
          int u1 = upair.first;
          int u2 = upair.second;
          int nlinks = nlinks_unibases[ upair ];
          bool flip1 = u2m[u1] < 0, flip2 = u2m[u2] < 0; 
          ForceAssert( ! (flip1 && flip2) );
          int m1 = flip1 ? u2m[ to_rc[u1] ] : u2m[u1];
          int m2 = flip2 ? u2m[ to_rc[u2] ] : u2m[u2];
          basevector unibase1( unibases[u1] );
          basevector unibase2( unibases[u2] );
          vec<int> accepted_overlaps;
          int max_overlap;

          ComputeOverlaps( unibase1, unibase2, u1, max_dist,
                    kmers_plus, kmers, 0 , accepted_overlaps, max_overlap );

          // Remodel the gap if there are enough links
          ulink_with_uids ulinks;
          int gap = 0;
          double dev = 0;
          Bool gap_predicted = false;
          if ( nlinks >= min_links_initial ) {
               if ( flip1 ) FlipAlign(m1, unibases_norc, pairs, base_lens, aligns_index, placed_fw, placed_rc, placement);
               if ( flip2 ) FlipAlign(m2, unibases_norc, pairs, base_lens, aligns_index, placed_fw, placed_rc, placement);
               PredictGap( "jump", kmers_plus, kmers, max_overlap, nlibs, 0, libs_to_use, lfail, unibases_norc,
                         DISTS, X, bases, pairs, base_lens, aligns_index, placed_fw, placed_rc, 
                         placement,  gap_id( gap_id::BETWEEN, m1, m2 ), 1, gap_predicted, gap, dev );
               if ( flip1 ) FlipAlign(m1, unibases_norc, pairs, base_lens, aligns_index, placed_fw, placed_rc, placement);
               if ( flip2 ) FlipAlign(m2, unibases_norc, pairs, base_lens, aligns_index, placed_fw, placed_rc, placement);
               if (  gap_predicted ) {
                    if ( accepted_overlaps.solo( ) )
                    {    int ogap = -accepted_overlaps[0];
                         if ( Abs( gap - ogap ) <= 3.0 * dev ) gap = ogap;    
                    }
               }
          }

          if ( VERBOSITY >= 4)
	    cout << "regapping " << u1 << "_" << u2 << "( "  << to_rc[u2] << "_" << to_rc[u1] << " )" 
	      " nlinks= " << nlinks << " newgap= " << 
	      ( gap_predicted ? ToString(gap) + " +/- " + ToString(dev) : String("failed") ) << endl;

	  // Initial gap.
	  ulinks = input_gaps[ upair ];

	  // Try to regap.
	  if ( REGAP ) {
	    if ( gap_predicted ) {
	      // only re-assign gap and dev
	      ulinks.sep = gap;
	      ulinks.dev = dev;
	    }
          }
	  else {
	    if (  gap_predicted ) {
	      int start1 = 0, stop2 = 0;  // those numbers are useless at the moment
	      ulinks = ulink_with_uids( u1, u2, gap, dev, start1, stop2, nlinks );
	    }
          }

	  // Unibases overlap - a suspect case? Use original gap.
          if ( ulinks.sep + max_devs * ulinks.dev <  -(K-1) ) {
	    if ( VERBOSITY >= 4 ) 
	      cout << "gap " << u1 << "_" << u2 << " " << "gap= " << ulinks.sep << " +/ " << ulinks.dev << " "
		   << " predicted gap size implies an overlap >= K-1 " << endl;
	  }
	  
          #pragma omp critical
          {  
               condensed_links_new.push( ulinks ); 
               // the gap fro u2* to u1* is the same
               ulink_with_uids ulinks_rc = ulinks;
               ulinks_rc.u1 = to_rc[u2];
               ulinks_rc.u2 = to_rc[u1];
               condensed_links_new.push( ulinks_rc ); 
          }
     }
     }
     
     // Assign the from- and to- data
     most_links_from.assign( nuni, 0 );
     most_links_to.assign( nuni, 0 );
     cov_from.assign(nuni,vec<triple<ho_interval, int, int> >() );
     cov_to.assign(nuni,  vec<triple<ho_interval, int, int> >() );
     cn_from.assign(nuni, 0); 
     cn_to.assign(nuni, 0);
     for ( int i = 0; i < condensed_links_new.isize(); i++ ) {
          const ulink_with_uids& ulinks = condensed_links_new[i];
          int nlinks = ulinks.nlinks;
          int u1 = ulinks.u1;
          int u2 = ulinks.u2;
          int gap = ulinks.sep;
          int dev = ulinks.dev;
          int Start2 = gap + dev_mult * dev;
          int Stop2 = gap + unibases[u2].isize( ) - dev_mult * dev;
          if ( Start2 < Stop2 ) 
          {    ho_interval h2( Start2, Stop2 );
               cov_from[u1].push( h2, nlinks, u2 );    
          }
          int Start1 = -gap - unibases[u1].isize( ) + dev_mult * dev;
          int Stop1 = -gap - dev_mult * dev;
          if ( Start1 < Stop1 ) 
          {    ho_interval h1( Start1, Stop1 );
               cov_to[u2].push( h1, nlinks, u1 );    
          }
          if ( nlinks > most_links_from[u1] )
          {    most_links_from[u1] = nlinks;
               cn_from[u1] = CN_raw[u2];    
          }
          if ( nlinks > most_links_to[u2] )
          {    most_links_to[u2] = nlinks;
               cn_to[u2] = CN_raw[u1];    
          }    
     }
     swap( condensed_links, condensed_links_new );
}

#endif
