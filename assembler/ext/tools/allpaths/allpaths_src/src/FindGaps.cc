// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

#include "AnAssemblyClass.h"
#include "AnAssemblyClassUtil.h"
#include "CoreTools.h"
#include "FindGaps.h"

// FindGaps: given a supercontig, with estimated gaps, recompute them (using
// values from RecomputeGaps to bootstrap).  Also compute standard deviations.

void FindGaps( assembly& A, const vec<int>& simple_reads_orig_index, int s, 
     vec<int>& gaps_dev, int start_contig, int stop_contig )
{
     RecomputeGaps( A, s, False, False );

     super &sup1 = A.supers[s];
     int n = sup1.mtig.size( );
     gaps_dev.resize_and_set(n-1, 0);
     static vec<Bool> gap_computed, bad_dev;
     gap_computed.resize_and_set( n-1, False );
     bad_dev.resize_and_set( n-1, False );

     // At first, in computing the size of the gap, use only the two contigs on 
     // either side of it.  Then use progressively more.

     for ( int tigs_to_use = Min(2, n-1); tigs_to_use < n; tigs_to_use++ )
     {    for ( int j = 1; j < n; j++ )
          {    if ( gap_computed[j-1] || bad_dev[j-1] ) continue;

               // Find the relevant links which cross the gap.

               static vec<int> loc_indices;
               loc_indices.clear( );

               vec<int> &mtig = A.supers[s].mtig;
               
               int begin_pos = Max( 0, j - tigs_to_use + 1 );
               int end_pos = Min( j + tigs_to_use - 1, (int) mtig.size() );

               for ( int k = begin_pos; k < j; k++ )
               {    int m1 = mtig[k];
                    const vec<int>& m1reads = A.reads_orig_index[m1];
                    for ( unsigned int l = 0; l < m1reads.size( ); l++ )
                    {    int loc_idx = m1reads[l];
                         const read_location& r1 = A.reads_orig[ loc_idx ];
                         if ( r1.OrientationOnContig( ) == ReverseOr ) continue;
                         int id1 = r1.ReadId( );
                         if ( A.pairs_index[id1] >= 0 )
                         {    const read_pairing& p = A.pairs[ A.pairs_index[id1] ];
                              int id2 = p.Partner(id1);
                              int si = simple_reads_orig_index[id2];
                              if ( si >= 0 )
                              {    const read_location& r2 = A.reads_orig[si];
                                   if ( r2.OrientationOnContig( ) != ReverseOr ) 
                                        continue;
                                   int m2 = r2.Contig( );
                                   vec<int>::iterator found_contig = 
                                        find( mtig.begin() + j, 
                                             mtig.begin() + end_pos, m2 );
                                   if ( found_contig != mtig.begin() + end_pos )
                                        loc_indices.push_back( 
                                             m1reads[l] );    }    }    }    }
               if ( loc_indices.empty( ) ) continue;

               // Set up supercontigs.

               int n1 = Min( tigs_to_use, j ), n2 = Min( tigs_to_use, n - j );
               static super sup1a, sup1b;
               static vec<int> deva, devb;
               sup1a.mtig.resize(n1);
               sup1a.gap.resize(n1-1);
               deva.resize(n1-1);
               for ( int k = 0; k < n1; k++ )
               {    sup1a.mtig[k] = sup1.mtig[j-n1+k];
                    if ( k < n1 - 1 ) 
                    {    sup1a.gap[k] = sup1.gap[j-n1+k];
                         deva[k] = gaps_dev[j-n1+k];    }    }
               sup1b.mtig.resize(n2);
               sup1b.gap.resize(n2-1);
               devb.resize(n2-1);
               for ( int k = 0; k < n2; k++ )
               {    sup1b.mtig[k] = sup1.mtig[j+k];
                    if ( k < n2 - 1 ) 
                    {    sup1b.gap[k] = sup1.gap[j+k];
                         devb[k] = gaps_dev[j+k];    }    }

               // Compute gaps.  Force short links to be used.  Check for
               // completely outrageous values in deva and devb first.  Such 
               // values have been observed to arise in totally screwed up
               // supercontigs, but could conceivably also occur in good supers.

               int gap = -1, gap_dev = -1;
               for ( int x = 0; x < (int) deva.size( ); x++ )
                    if ( deva[x] > 1000000 ) bad_dev[j-1] = True;
               for ( int x = 0; x < (int) devb.size( ); x++ )
                    if ( devb[x] > 1000000 ) bad_dev[j-1] = True;
               if ( !bad_dev[j-1] )
                    gap = ExpectedGapSuperFull( A, sup1a, deva, sup1b, devb, 
                         loc_indices, gap_dev, False, False, False, True,
                         simple_reads_orig_index );

               if ( gap_dev > 0 )
               {    
                    // If we're not between start_contig and stop_contig, use the
                    // original gap.

                    if ( start_contig <= j-1 && j <= stop_contig ) 
                         sup1.gap[j-1] = gap;
                    else gap = sup1.gap[j-1];

                    gaps_dev[j-1] = gap_dev;

                    gap_computed[j-1] = True;    }    }    }

     // If all else fails, use the original gaps, and a guess for the
     // standard deviation.

     for ( int j = 1; j < n; j++ )
     {    if ( !gap_computed[j-1] )
          {    sup1.gap[j-1] = sup1.gap[j-1];
               gaps_dev[j-1] = Abs(sup1.gap[j-1]/4);    }    }    }


void GapStats( vec<int> gap, vec<int> gapdev, int& gap_ave, int& gapdev_ave )
{
     // If there are less than six gaps, we directly compute their mean.
     // Otherwise, we attempt to remove outliers, as follows.  We sort 
     // the gaps and extract the middle half.  From this middle half, we compute the 
     // mean and standard deviation.  Then we select those gaps lying withing 5 
     // standard deviations of the mean, and form their mean.  

     if ( gap.size( ) >= 6 )
     {    vec<int> mid_gaps;
          Sort(gap);
          for ( unsigned int i = gap.size( )/4; i < 3*(1+gap.size( ))/4; i++ )
               mid_gaps.push_back( gap[i] );
          float sum1 = 0, sum2 = 0;
          for ( unsigned int i = 0; i < mid_gaps.size( ); i++ )
          {    sum1 += mid_gaps[i];
               sum2 += float(mid_gaps[i]) * float(mid_gaps[i]);    }
          float n = mid_gaps.size( );
          float mean = sum1/n;
          float sd = sqrt(sum2/n - mean * mean);
          float start = mean - 5 * sd, stop = mean + 5 * sd;
          longlong gap_count = 0, gap_total = 0;
          Sort(gap);
          float x = 0;
          for ( unsigned int i = 0; i < gap.size( ); i++ )
          {    if ( start <= gap[i] && gap[i] <= stop )
               {    ++gap_count;
                    gap_total += gap[i];    
                    x +=  float(gapdev[i]) * float(gapdev[i]);    }    }
          if ( gap_count > 0 )
          {    gap_ave = gap_total / gap_count;
               gapdev_ave = int(floor( sqrt(x) / gapdev.size( ) ));
               return;    }    }
     float x = 0;
     longlong total = BigSum(gap);
     gap_ave = total / longlong(gap.size( ));
     for ( unsigned int i = 0; i < gapdev.size( ); i++ )
          x += float(gapdev[i]) * float(gapdev[i]);
     gapdev_ave = int(floor( sqrt(x) / gapdev.size( ) ));    }

void Expect( assembly& a, int m1, int m2, int& gap_ave, int& gapdev_ave,
     Bool strong, Bool rc1, Bool rc2 )
{
     // Find the links from m1 to m2.

     static vec<int> gap, gapdev;
     gap.clear( );
     gapdev.clear( );
     vec<int>& m1reads = a.reads_orig_index[m1];
     for ( unsigned int i = 0; i < m1reads.size( ); i++ )
     {    read_location& r1 = a.reads_orig[ m1reads[i] ];
          int id1 = r1.ReadId( );
          if ( a.pairs_index[id1] >= 0 )
          {    read_pairing& p = a.pairs[ a.pairs_index[id1] ];
               int id2 = p.Partner(id1);
               int ri = a.simple_reads_orig_index[id2];
               if ( ri < 0 ) continue;
               read_location& r2 = a.reads_orig[ri];
               if ( r2.Contig( ) != m2 ) continue;

               if ( ( r1.OrientationOnContig( ) + r2.OrientationOnContig( )
                    + rc1 + rc2 ) % 2 == 1 )
               {    if ( !strong || 
                         ( r1.OrientationOnContig( ) ^ rc1 == ForwardOr &&
                           r2.OrientationOnContig( ) ^ rc2 == ReverseOr ) )
                    {    int end1, end2;
                         if ( !rc1 ) 
                              end1 = r1.LengthOfContig( ) - r1.StopOnContig( );
                         else end1 = r1.StartOnContig( );
                         if ( !rc2 ) end2 = r2.StartOnContig( );
                         else end2 = r2.LengthOfContig( ) - r2.StopOnContig( );
                         gap.push_back( p.sep - end1 - end2 );
                         gapdev.push_back( p.sd );    }    }    }    }

     if ( gap.size( ) == 0 ) // XXX
     {    cout << "No pairings found from " << m1 << " to " << m2 // XXX
               << ".\n"; // XXX
          exit(1);    } // XXX
     ForceAssert( gap.size( ) > 0 );
     GapStats( gap, gapdev, gap_ave, gapdev_ave );    }


void Expect2( assembly& a, int m1, int m2, int& gap_ave, int& gapdev_ave,
     Bool rc1, Bool rc2, Bool exit_if_no_pairings )
{
     ForceAssert( rc1 == False );

     // Find the links from m1 to m2.

     static vec<int> gap, gapdev;
     gap.clear( );
     gapdev.clear( );
     vec<int>& m1reads = a.reads_orig_index[m1];
     for ( unsigned int i = 0; i < m1reads.size( ); i++ )
     {    const read_location& r1 = a.reads_orig[ m1reads[i] ];
          int id1 = r1.ReadId( );
          if ( a.pairs_index[id1] >= 0 )
          {    read_pairing& p = a.pairs[ a.pairs_index[id1] ];
               int id2 = p.Partner(id1);
               int ri = a.simple_reads_orig_index[id2];
               if ( ri < 0 ) continue;
               read_location& r2 = a.reads_orig[ri];
               if ( r2.Contig( ) != m2 ) continue;

               // We have a pairing p between contigs m1 and m2 
               // (or their reverse complements).

               if ( !rc2 )
               {    if ( r1.OrientationOnContig( ) != r2.OrientationOnContig( ) )
                    {    
                         if ( r1.OrientationOnContig( ) == ForwardOr )
                         {    int end1 = r1.LengthOfContig( ) - r1.StopOnContig( );
                              int end2 = r2.StartOnContig( );
                              gap.push_back( p.sep - end1 - end2 );
                              gapdev.push_back( p.sd );    }
     
                         else // contigs out of order!
                         {
                              // First compute the gap you would expect if the
                              // contigs were in the right order.
     
                              int end1 = r2.LengthOfContig( ) - r2.StopOnContig( );
                              int end2 = r1.StartOnContig( );
                              int g = p.sep - end1 - end2;
     
                              // Now adjust.
     
			      g = - g - a.mtig[m1].size( ) - a.mtig[m2].size( );
			      gap.push_back(g);
			      gapdev.push_back( p.sd );  }  }    
               else
               {    if ( r1.OrientationOnContig( ) == r2.OrientationOnContig( ) )
                    {    
                         if ( r1.OrientationOnContig( ) == ForwardOr )
                         {    int end1 = r1.LengthOfContig( ) - r1.StopOnContig( );
                              int end2 = r2.LengthOfContig( ) - r2.StopOnContig( );
                              gap.push_back( p.sep - end1 - end2 );
                              gapdev.push_back( p.sd );    }
     
                         else // contigs out of order!
                         {
                              // First compute the gap you would expect if the
                              // contigs were in the right order.
     
                              int end1 = r2.LengthOfContig( ) - r2.StopOnContig( );
                              int end2 = r1.LengthOfContig( ) - r1.StopOnContig( );
                              int g = p.sep - end1 - end2;
     
                              // Now adjust.
     
                              g = - g - a.mtig[m1].size( ) - a.mtig[m2].size( );
                              gap.push_back(g);
                              gapdev.push_back( p.sd );    }    }    }    }    }    }

     
     if ( gap.size( ) == 0 ) // XXX
     {   
       // cout << "No pairings found from " << m1 << " to " << m2 << ".\n"; // XXX
     
       if ( exit_if_no_pairings )
	 exit(1);
       else
       {
	 gap_ave=-1;
	 gapdev_ave=-1;
       }
     } // XXX
     else
     {
       ForceAssert( gap.size( ) > 0 );
       GapStats( gap, gapdev, gap_ave, gapdev_ave );   
     }
}


int ExpectedGap( assembly& a, int m1, int m2 )
{    int gap_ave, gapdev_ave;
     Expect( a, m1, m2, gap_ave, gapdev_ave );
     return gap_ave;    }

int ExpectedGap2( assembly& a, int m1, int m2 )
{    int gap_ave, gapdev_ave;
     Expect2( a, m1, m2, gap_ave, gapdev_ave );
     return gap_ave;    }
     

int ExpectedGapSuper( assembly& a, int m1, int m2, int& gap_dev, 
     Bool rc1, Bool rc2, Bool strong )
{    
     int s1 = a.mtigs_to_supers[m1], s2 = a.mtigs_to_supers[m2];
     super &sup1 = a.supers[s1], &sup2 = a.supers[s2];
     int n1 = sup1.mtig.size( ), n2 = sup2.mtig.size( );

     // Compute the lengths of the mtigs.

     vec<int> len1(n1), len2(n2);
     for ( int i = 0; i < n1; i++ )
          len1[i] = a.mtig[ sup1.mtig[i] ].size( );
     for ( int i = 0; i < n2; i++ )
          len2[i] = a.mtig[ sup2.mtig[i] ].size( );

     // Compute the starting positions of the mtigs on the supercontigs.

     vec<int> pos1(n1), pos2(n2);
     pos1[0] = 0;
     pos2[0] = 0;
     for ( int i = 1; i < n1; i++ )
          pos1[i] = pos1[i-1] + len1[i-1] + sup1.gap[i-1];
     for ( int i = 1; i < n2; i++ )
          pos2[i] = pos2[i-1] + len2[i-1] + sup2.gap[i-1];

     // Compute the indices of the mtigs in the supercontigs.

     unsigned int i1, i2;
     for ( i1 = 0; i1 < sup1.mtig.size( ); i1++ )
          if ( sup1.mtig[i1] == m1 ) break;
     for ( i2 = 0; i2 < sup2.mtig.size( ); i2++ )
          if ( sup2.mtig[i2] == m2 ) break;

     // Compute the distance from the end of m1 to the end of s1,
     // and the distance from the end of m2 to the end of s2.

     int to_end1 = pos1.back( ) + len1.back( ) - pos1[i1] - len1[i1];
     int to_end2 = pos2.back( ) + len2.back( ) - pos2[i2] - len2[i2];

     // Compute the distance from the beginning of s1 to the beginning of m1,
     // and the distance from the beginning of s2 to the beginning of m2.

     int to_begin1 = pos1[i1], to_begin2 = pos2[i2];

     // Compute the distance from the end of m1 to the end of s1,
     // and the distance from the beginning of s2 to the beginning of m2,
     // taking into account orientations.

     int dist1 = (!rc1) ? to_end1 : to_begin1;
     int dist2 = rc2 ? to_end2 : to_begin2;

     int gap_ave;
     Expect( a, m1, m2, gap_ave, gap_dev, strong, rc1, rc2 );
     return gap_ave - dist1 - dist2;    }

int ExpectedGapSuperFull( assembly& a, const super& sup1, const vec<int>& dev1,
     const super& sup2, const vec<int>& dev2, const vec<int>& read_indices, 
     int& gap_dev_ave, Bool rc1, Bool rc2, Bool verbose, Bool force_short,
     const vec<int>& simple_reads_orig_index, int min_gap, int max_gap )
{    
     if (verbose) cout << "\nentering ExpectedGapSuperFull\n";
     const vec<int> &m1 = sup1.mtig, &m2 = sup2.mtig;
     int n1 = m1.size( ), n2 = m2.size( );

     // Create indices to the contigs in sup1 and sup2.

     static vec<int> cpos1a, cpos1b, cpos2a, cpos2b;
     {    cpos1a.resize(n1);
          cpos1b.resize(n1);
          cpos2a.resize(n2);
          cpos2b.resize(n2);
          static vec< pair<int, int> > contig_pos1, contig_pos2;
          contig_pos1.resize(n1);
          contig_pos2.resize(n2);
          for ( int i = 0; i < n1; i++ )
               contig_pos1[i] = make_pair( m1[i], i );
          for ( int i = 0; i < n2; i++ )
               contig_pos2[i] = make_pair( m2[i], i );
          Sort(contig_pos1), Sort(contig_pos2);
          for ( int i = 0; i < n1; i++ )
          {    cpos1a[i] = contig_pos1[i].first;
               cpos1b[i] = contig_pos1[i].second;    }
          for ( int i = 0; i < n2; i++ )
          {    cpos2a[i] = contig_pos2[i].first;
               cpos2b[i] = contig_pos2[i].second;    }    }

     // Compute the lengths of the mtigs.

     static vec<int> len1, len2;
     len1.resize(n1);
     len2.resize(n2);
     for ( int i = 0; i < n1; i++ )
          len1[i] = a.mtig[ m1[i] ].size( );
     for ( int i = 0; i < n2; i++ )
          len2[i] = a.mtig[ m2[i] ].size( );

     // Compute the starting positions of the mtigs on the supercontigs.

     static vec<int> pos1, pos2;
     pos1.resize(n1);
     pos2.resize(n2);
     pos1[0] = 0;
     pos2[0] = 0;
     for ( int i = 1; i < n1; i++ )
          pos1[i] = pos1[i-1] + len1[i-1] + sup1.gap[i-1];
     for ( int i = 1; i < n2; i++ )
          pos2[i] = pos2[i-1] + len2[i-1] + sup2.gap[i-1];

     // Compute the distance from the end of m1[i] to the end of sup1,
     // and the distance from the end of m2[i] to the end of sup2.

     static vec<int> to_end1, to_end2;
     to_end1.resize(n1);
     to_end2.resize(n2);
     for ( int i1 = 0; i1 < n1; i1++ )
          to_end1[i1] = pos1.back( ) + len1.back( ) - pos1[i1] - len1[i1];
     for ( int i2 = 0; i2 < n2; i2++ )
          to_end2[i2] = pos2.back( ) + len2.back( ) - pos2[i2] - len2[i2];

     // Compute the distance from the beginning of sup1 to the beginning of m1[i1],
     // and the distance from the beginning of sup2 to the beginning of m2[i2].

     static vec<int> to_begin1, to_begin2;
     to_begin1.resize(n1);
     to_begin2.resize(n2);
     for ( int i1 = 0; i1 < n1; i1++ )
          to_begin1[i1] = pos1[i1];
     for ( int i2 = 0; i2 < n2; i2++ )
          to_begin2[i2] = pos2[i2];

     // Compute the distance from the end of m1[i] to the end of sup1,
     // and the distance from the beginning of sup2 to the beginning of m2[i],
     // taking into account orientations.

     static vec<int> dist1, dist2;
     dist1.resize(n1);
     dist2.resize(n2);
     for ( int i1 = 0; i1 < n1; i1++ )
          dist1[i1] = (!rc1) ? to_end1[i1] : to_begin1[i1];
     for ( int i2 = 0; i2 < n2; i2++ )
          dist2[i2] = rc2 ? to_end2[i2] : to_begin2[i2];

     // Compute deviations for these.

     static vec<int> dist1_dev, dist2_dev;
     dist1_dev.resize_and_set(n1, 0);
     dist2_dev.resize_and_set(n2, 0);
     if ( !dev1.empty( ) )
     {    static vec<int> pos1_dev, pos2_dev;
          pos1_dev.resize(n1);
          pos2_dev.resize(n2);
          pos1_dev[0] = 0;
          pos2_dev[0] = 0;
          for ( int i = 1; i < n1; i++ )
               pos1_dev[i] = pos1_dev[i-1] + dev1[i-1];
          for ( int i = 1; i < n2; i++ )
               pos2_dev[i] = pos2_dev[i-1] + dev2[i-1];
          static vec<int> to_end1_dev, to_end2_dev;
          to_end1_dev.resize(n1);
          to_end2_dev.resize(n2);
          for ( int i1 = 0; i1 < n1; i1++ )
               to_end1_dev[i1] = pos1_dev.back( ) - pos1_dev[i1];
          for ( int i2 = 0; i2 < n2; i2++ )
               to_end2_dev[i2] = pos2_dev.back( ) - pos2_dev[i2];
          static vec<int> to_begin1_dev, to_begin2_dev;
          to_begin1_dev.resize(n1);
          to_begin2_dev.resize(n2);
          for ( int i1 = 0; i1 < n1; i1++ )
               to_begin1_dev[i1] = pos1_dev[i1];
          for ( int i2 = 0; i2 < n2; i2++ )
               to_begin2_dev[i2] = pos2_dev[i2];
          dist1_dev.resize(n1);
          dist2_dev.resize(n2);
          for ( int i1 = 0; i1 < n1; i1++ )
               dist1_dev[i1] = (!rc1) ? to_end1_dev[i1] : to_begin1_dev[i1];
          for ( int i2 = 0; i2 < n2; i2++ )
               dist2_dev[i2] = rc2 ? to_end2_dev[i2] : to_begin2_dev[i2];    }

     // Start computing gaps.

     static vec<int> gap, gapdev;
     gap.clear( );
     gapdev.clear( );

     for ( int j1 = 0; j1 < (int) read_indices.size( ); j1++ )
     {    const read_location& r1 = a.reads_orig[ read_indices[j1] ];
          int i1 = cpos1b[ BinPosition( cpos1a, r1.Contig( ) ) ];
          int id1 = r1.ReadId( );
          if ( a.pairs_index[id1] >= 0 )
          {    const read_pairing& p = a.pairs[ a.pairs_index[id1] ];
               int id2 = p.Partner(id1);
               int ri = simple_reads_orig_index[id2];
               if ( ri < 0 ) continue;
               const read_location& r2 = a.reads_orig[ri];
               if ( !BinMember( cpos2a, r2.Contig( ) ) ) continue;
               int i2 = cpos2b[ BinPosition( cpos2a, r2.Contig( ) ) ];

               // We have a pairing p between contigs m1[i1] and m2[i2] 
               // (or their reverse complements).

               int m1_length = r1.LengthOfContig( );
               int m2_length = r2.LengthOfContig( );
               orientation r1_orient = r1.OrientationOnContig( );
               orientation r2_orient = r2.OrientationOnContig( );

               if ( !rc1 )
               {    if ( !rc2 )
                    {    if ( r1_orient != r2_orient )
                         {    
                              if ( r1_orient == ForwardOr )
                              {    int end1 = m1_length - r1.StopOnContig( );
                                   int end2 = r2.StartOnContig( );
                                   gap.push_back( p.sep - end1 - end2 
                                        - dist1[i1] - dist2[i2] );
                                   if (verbose) 
                                   {    PRINT4( 1, id1, id2, gap.back( ) );
                                        PRINT4( end1, end2, dist1[i1], dist2[i2] );
                                             }
                                   gapdev.push_back( 
                                        p.sd + dist1_dev[i1] + dist2_dev[i2] );    }
     
                              else // contigs out of order!
                              {
                                   // First compute the gap you would expect if the
                                   // contigs were in the right order.
          
                                   int end1 = m2_length - r2.StopOnContig( );
                                   int end2 = r1.StartOnContig( );
                                   int g = p.sep - end1 - end2;
          
                                   // Now adjust.
               
                                   g = - g - a.mtig[ m1[i1] ].size( ) 
                                        - a.mtig[ m2[i2] ].size( );
                                   gap.push_back( g - dist1[i1] - dist2[i2] );
                                   if (verbose) 
                                   {    PRINT4( 2, id1, id2, gap.back( ) );
                                        PRINT4( end1, end2, dist1[i1], dist2[i2] );
                                             }
                                   gapdev.push_back(
                                        p.sd + dist1_dev[i1] + dist2_dev[i2] );    
                                             }    }    }
                    else
                    {    if ( r1_orient == r2_orient )
                         {    
                              if ( r1_orient == ForwardOr )
                              {    int end1 = m1_length - r1.StopOnContig( );
                                   int end2 = m2_length - r2.StopOnContig( );
                                   gap.push_back( p.sep - end1 - end2 
                                        - dist1[i1] - dist2[i2] );
                                   if (verbose) 
                                   {    PRINT4( 3, id1, id2, gap.back( ) );
                                        PRINT4( end1, end2, dist1[i1], dist2[i2] );
                                             }
                                   gapdev.push_back(
                                        p.sd + dist1_dev[i1] + dist2_dev[i2] );    }
                    
                              else // contigs out of order!
                              {
                                   // First compute the gap you would expect if 
                                   // the contigs were in the right order.
               
                                   int end1 = m2_length - r2.StopOnContig( );
                                   int end2 = m1_length - r1.StopOnContig( );
                                   int g = p.sep - end1 - end2;
               
                                   // Now adjust.
                    
                                   g = - g - a.mtig[ m1[i1] ].size( ) 
                                        - a.mtig[ m2[i2] ].size( );
                                   gap.push_back( g - dist1[i1] - dist2[i2] );
                                   if (verbose) 
                                   {    PRINT4( 4, id1, id2, gap.back( ) );
                                        PRINT4( end1, end2, dist1[i1], dist2[i2] );
                                             }
                                   gapdev.push_back( 
                                        p.sd + dist1_dev[i1] + dist2_dev[i2] );    
                                             }    }    }    }
               else
               {    if ( !rc2 )
                    {    if ( r1_orient == r2_orient )
                         {    
                              if ( r1_orient != ForwardOr )
                              {    int end1 = r1.StartOnContig( );
                                   int end2 = r2.StartOnContig( );
                                   gap.push_back( p.sep - end1 - end2 
                                        - dist1[i1] - dist2[i2] );
                                   if (verbose) 
                                   {    PRINT4( 5, id1, id2, gap.back( ) );
                                        PRINT4( end1, end2, dist1[i1], dist2[i2] );
                                             }
                                   gapdev.push_back( 
                                        p.sd + dist1_dev[i1] + dist2_dev[i2] );    }
     
                              else // contigs out of order!
                              {
                                   // First compute the gap you would expect if the
                                   // contigs were in the right order.
          
                                   int end1 = m2_length - r2.StopOnContig( );
                                   int end2 = m1_length - r1.StopOnContig( );
                                   int g = p.sep - end1 - end2;
          
                                   // Now adjust.
               
                                   g = - g - a.mtig[ m1[i1] ].size( ) 
                                        - a.mtig[ m2[i2] ].size( );
                                   gap.push_back( g - dist1[i1] - dist2[i2] );
                                   if (verbose) 
                                   {    PRINT4( 6, id1, id2, gap.back( ) );
                                        PRINT4( end1, end2, dist1[i1], dist2[i2] );
                                             }
                                   gapdev.push_back(
                                        p.sd + dist1_dev[i1] + dist2_dev[i2] );    
                                         }    }    }
                    else
                    {    if ( r1_orient != r2_orient )
                         {    
                              if ( r1_orient != ForwardOr )
                              {    int end1 = r1.StartOnContig( );
                                   int end2 = m2_length - r2.StopOnContig( );
                                   gap.push_back( p.sep - end1 - end2 
                                        - dist1[i1] - dist2[i2] );
                                   if (verbose) 
                                   {    PRINT4( 7, id1, id2, gap.back( ) );
                                        PRINT4( end1, end2, dist1[i1], dist2[i2] );
                                             }
                                   gapdev.push_back( 
                                        p.sd + dist1_dev[i1] + dist2_dev[i2] );    }
                    
                              else // contigs out of order!
                              {
                                   // First compute the gap you would expect if 
                                   // the contigs were in the right order.
               
                                   int end1 = m2_length - r2.StopOnContig( );
                                   int end2 = r1.StartOnContig( );
                                   int g = p.sep - end1 - end2;
               
                                   // Now adjust.
                    
                                   g = - g - a.mtig[ m1[i1] ].size( ) 
                                        - a.mtig[ m2[i2] ].size( );
                                   gap.push_back( g - dist1[i1] - dist2[i2] );
                                   if (verbose) 
                                   {    PRINT4( 8, id1, id2, gap.back( ) );
                                        PRINT4( end1, end2, dist1[i1], dist2[i2] );
                                             }
                                   gapdev.push_back( 
                                        p.sd + dist1_dev[i1] + dist2_dev[i2] );    
                                             }    }    }    }    }    }

     if ( gap.size( ) == 0 )
     {    // cout << "No pairings found between supercontigs.\n";
          gap_dev_ave = -1;
          return 0;    }

     // Filter out gaps outside [min_gap, max_gap] if there are any inside.

     Bool have_good_gap = False, have_bad_gap = False;
     for ( int i = 0; i < gap.isize( ); i++ )
     {    if ( min_gap <= gap[i] && gap[i] <= max_gap ) have_good_gap = True;
          if ( min_gap > gap[i] || gap[i] > max_gap ) have_bad_gap = True;    }
     if ( have_good_gap && have_bad_gap )
     {    static vec<int> gap_new, gapdev_new;
          gap_new.clear( ), gapdev_new.clear( );
          for ( int i = 0; i < gap.isize( ); i++ )
          {    if ( min_gap <= gap[i] && gap[i] <= max_gap )
               {    gap_new.push_back( gap[i] );
                    gapdev_new.push_back( gapdev[i] );    }    }
          gap = gap_new;
          gapdev = gapdev_new;    }

     // The following heuristic is designed so that if credible short links are
     // available, they will be used exclusively.  In most cases, this makes the
     // answer more accurate.
     //
     // Compute the mean and standard deviation (gap_ave, gap_dev_ave) if all links 
     // are used.  Then let d_min be the minimum s.d. appearing amongst the links.
     // Let d = Max( 2000, 3 * d_min ).  Now consider only those links having s.d.
     // <= d.  Compute gap_ave_sub, gap_dev_ave_sub.  If |gap_ave - gap_ave_sub|
     // <= 3 * gap_dev_ave, return the sub values instead.

     int gap_ave;
     GapStats( gap, gapdev, gap_ave, gap_dev_ave );
     int d = Max( 2000, 3 * Min(gapdev) );
     static vec<int> gap_sub, gapdev_sub;
     gap_sub.clear( );
     gapdev_sub.clear( );
     for ( int i = 0; i < (int) gap.size( ); i++ )
          if ( gapdev[i] <= d )
          {    gap_sub.push_back( gap[i] );
               gapdev_sub.push_back( gapdev[i] );    }
     int gap_ave_sub, gap_dev_ave_sub;
     GapStats( gap_sub, gapdev_sub, gap_ave_sub, gap_dev_ave_sub );
     if ( ( force_short && gap_sub.size( ) >= 2 )
              || Abs( gap_ave - gap_ave_sub ) <= 3 * gap_dev_ave )
     {    gap_ave = gap_ave_sub;
          gap_dev_ave = gap_dev_ave_sub;    }
     return gap_ave;    }
