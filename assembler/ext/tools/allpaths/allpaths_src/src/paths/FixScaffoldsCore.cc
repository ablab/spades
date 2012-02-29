///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include <omp.h>

#include "Fastavector.h"
#include "Charvector.h"
#include "PairsManager.h"
#include "ParseSet.h"
#include "Superb.h"
#include "VecTemplate.h"
#include "lookup/LookAlign.h"
#include "math/HoInterval.h"
#include "math/NStatsTools.h"
#include "paths/Alignlet.h"
#include "paths/ContigsManager.h"
#include "paths/FixScaffoldsCore.h"
#include "paths/ScaffoldsUtils.h"
// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

class blink_fw {

     public:

     int start1; // start on t1
     int stop1_low, stop1_high; // where read2 would land if it did
     int t2;
     int stop2;  // stop on t2

     friend Bool operator<( const blink_fw& b1, const blink_fw& b2 )
     {    if ( b1.t2 < b2.t2 ) return True;
          if ( b1.t2 > b2.t2 ) return False;
          if ( b1.start1 < b2.start1 ) return True;
          return False;    }

};

class blink_rc {

     public:

     int stop1; // stop on t1
     int start1_low, start1_high; // where read2 would land if it did
     int t2;
     int start2;  // start on t2

     friend Bool operator<( const blink_rc& b1, const blink_rc& b2 )
     {    if ( b1.t2 < b2.t2 ) return True;
          if ( b1.t2 > b2.t2 ) return False;
          if ( b1.stop1 < b2.stop1 ) return True;
          return False;    }

};

void FixScaffoldsCore( const vec<int> &trace_ids,
		       const PairsManager &pairs,
		       const int MIN_REACH_AWAY,
		       vec<fastavector> &contigs,
		       vec<superb> &scaffolds,
		       vec<alignlet> &aligns0,
		       vec<int> &aligns0_index,
		       vec<int> &aligns0_index_unfilt,
		       ostream &log,
		       bool VERBOSE ){
  
  vec<alignlet> ualigns0( 0 );
  vec< vec<int> > ualigns0_index( 0 );
  FixScaffoldsCore( trace_ids, pairs, MIN_REACH_AWAY, contigs, scaffolds,
		    aligns0, aligns0_index, aligns0_index_unfilt,
		    ualigns0, ualigns0_index, log, VERBOSE );

}


void FixScaffoldsCore( const vec<int> &trace_ids,
		       const PairsManager &pairs,
		       const int MIN_REACH_AWAY,
		       vec<fastavector> &contigs,
		       vec<superb> &scaffolds,
		       vec<alignlet> &aligns0,
		       vec<int> &aligns0_index,
		       vec<int> &aligns0_index_unfilt,
		       vec<alignlet> &ualigns0,
		       vec< vec<int> > &ualigns0_index,
		       ostream &log,
		       bool VERBOSE )
{
     // Define heuristic constants.

     const int min_dist_from_end = 1000;
     const double dev_mult = 4.0;
     const int blink_dev_mult = 1;
     const int min_scaffold = 1000;
     const int trim_back = 80;
     const int min_contig = 1000;
     const int min_spread = 10;

     
     //check alignments
     for ( size_t i = 0; i < ualigns0_index.size(); i++ )
       for ( size_t j = 0; j < ualigns0_index[i].size(); j++ )
	 if ( ualigns0_index[i][j] >= 0 ){
	   size_t ai = ualigns0_index[i][j];
	   ForceAssertLe( ualigns0[ ai ].Pos2(), (int)contigs[ ualigns0[ ai ].TargetId() ].size() );
	   ForceAssertGe( ualigns0[ ai ].pos2(), 0 );
	 }


     
     int lib_id = -1;
     int lib_size = 0;
     for (size_t ii=0; ii<pairs.nLibraries( ); ii++) {
       if ( lib_id < 0 || pairs.getLibrarySep( ii ) > lib_size ) {
	 lib_id = ii;
	 lib_size = pairs.getLibrarySep( ii );
       }
     }
     const int max_gap = lib_size + ( 3 * pairs.getLibrarySD( lib_id ) );
     
     // Load.
     
     // !!!! BE AWARE that ntigs calculated maybe smaller than the number of contigs, 
     // because in some cases the scaffolds do not include all contigs. To be safe,
     // let's include every contigs here.
     int ntigs = contigs.size();
     //int ntigs = 0;
     //for (int ii=0; ii<scaffolds.isize( ); ii++)
     //  for (int jj=0; jj<scaffolds[ii].Ntigs( ); jj++)
     //    ntigs = Max( ntigs, scaffolds[ii].Tig( jj ) );
     //ntigs++;

     vec<int> tig_sizes( ntigs, 0 );
     for (int ii=0; ii<scaffolds.isize( ); ii++)
       for (int jj=0; jj<scaffolds[ii].Ntigs( ); jj++)
	 tig_sizes[ scaffolds[ii].Tig( jj ) ] = scaffolds[ii].Len( jj );
     
     size_t nreads = aligns0_index.size( );

     // Index scaffolds.

     log << Date( ) << ": indexing scaffolds" << endl;
     vec<int> to_scaffold(ntigs, -1), to_scaffold_pos(ntigs, -1);
     for ( int i = 0; i < scaffolds.isize( ); i++ ) {
       const superb& s = scaffolds[i];
       for ( int j = 0; j < s.Ntigs( ); j++ ) {
	 ForceAssertEq( to_scaffold[ s.Tig(j) ], -1 );
	 to_scaffold[ s.Tig(j) ] = i;
	 to_scaffold_pos[ s.Tig(j) ] = j;
       }
     }

	 
     // Find pairs that land logically on a scaffold, and mark the physically
     // covered parts.

     log << Date( ) << ": finding coverage" << endl;
     size_t nbatches = 10 * omp_get_max_threads( );
     vec< vec< vec<ho_interval> > > covi(nbatches);
     for ( size_t i = 0; i < nbatches; i++ )
          covi[i].resize(ntigs);
     #pragma omp parallel for
     for ( size_t bi = 0; bi < nbatches; bi++ ){    
       size_t start = ( bi * nreads ) / nbatches;
       size_t stop = ( (bi+1) * nreads ) / nbatches;
       for ( size_t id1 = start; id1 < stop; id1++ ){    
	 if ( !pairs.isPaired(id1) ) continue;
	 int id2 = pairs.getPartnerID(id1);
	 if ( aligns0_index[id1] < 0 || aligns0_index[id2] < 0 ) continue;
	 alignlet la1 = aligns0[ aligns0_index[id1] ];
	 alignlet la2 = aligns0[ aligns0_index[id2] ];
	 int t1 = la1.TargetId( ), t2 = la2.TargetId( );
	 if ( to_scaffold[t1] < 0 || to_scaffold[t2] < 0 ) continue;
	 int sid = to_scaffold[t1];
	 if ( to_scaffold[t2] != sid || !la1.Fw1( ) || la2.Fw1( ) ) continue;
	 const superb& s = scaffolds[sid];
	 int p1 = to_scaffold_pos[t1], p2 = to_scaffold_pos[t2];
	 if ( !( p1 <= p2 ) ) continue;
	 int start = la1.pos2( ), stop = la2.Pos2( );
	 if ( p1 == p2 && !( start <= stop ) ) continue;
	 int psep = pairs.sep( pairs.getPairID(id1) );
	 int pdev = pairs.sd( pairs.getPairID(id1) );
	 int sep, dev = pdev;
	 if ( p1 == p2 ) sep = la2.pos2( ) - la1.Pos2( );
	 else{    
	   sep = s.SubSuperLength( p1, p2 ) - s.Len(p1) - s.Len(p2)
	     + la2.pos2( ) + tig_sizes[t1] - la1.Pos2( );
	   dev += s.SubSuperLengthDev( p1, p2 );    
	 }
	 if ( double( Abs( sep - psep ) ) > dev_mult * double(dev) ) continue;
	 start += trim_back;
	 stop -= trim_back;
	 if ( p1 == p2 ){
	   if ( start <= stop ) covi[bi][t1].push( start, stop );
	 }else{    
	   if ( start <= tig_sizes[t1] ) 
	     covi[bi][t1].push( start, tig_sizes[t1] );
	   for ( int j = p1 + 1; j < p2; j++ )
	     covi[bi][ s.Tig(j) ].push( 0, tig_sizes[ s.Tig(j) ] );
	   if ( stop >= 0 ) covi[bi][t2].push( 0, stop );    
	 }    
       }    
     }
     log << Date( ) << ": combining coverage" << endl;
     vec< vec<ho_interval> > cov(ntigs);
     #pragma omp parallel for
     for ( int t = 0; t < ntigs; t++ ){    
       for ( size_t bi = 0; bi < nbatches; bi++ )
	 cov[t].append( covi[bi][t] );    
     }

     // Find uncovered parts of each edge.

     log << Date( ) << ": finding uncovered parts" << endl;
     for ( int i = 0; i < trace_ids.isize( ); i++ ){    
       int t = trace_ids[i];
       Sort( cov[t] );
       log << "\ncoverage of " << t << endl;
       for ( int j = 0; j < cov[t].isize( ); j++ )
	 log << cov[t][j] << "\n";    
     }
     vec< vec<ho_interval> > un(ntigs);
     vec< vec<ho_interval> > covregs(ntigs);
     #pragma omp parallel for
     for ( int t = 0; t < ntigs; t++ ){
       Uncovered( tig_sizes[t], cov[t], un[t] );
       Uncovered( tig_sizes[t], un[t], covregs[t] );
     }
     

     // Find out if uncovered regions correspond to high copy numer unipath alignments
     // as found in ualigns0

     vec< vec<ho_interval> > hcns(ntigs);
     for ( size_t i = 0; i < ualigns0_index.size(); i++ )
       for ( size_t j = 0; j < ualigns0_index[i].size(); j++ ){
	 if ( ualigns0_index[i][j] < 0 ) continue;
	 int ai = ualigns0_index[i][j];
	 
	 int tid = ualigns0[ai].TargetId();
	 if ( ualigns0[ai].Pos2() - ualigns0[ai].pos2() < 200 ) continue; // not attempting breakes 
	                                                                  // around short highCN unibases
	 hcns[ tid ].push_back( ho_interval(ualigns0[ai].pos2(), ualigns0[ai].Pos2()) );
     }
     for ( int t = 0; t < ntigs; t++ )
       Sort( hcns[t] );
     for ( int t = 0; t < ntigs; t++ )
       for ( int i = 0; i < hcns[t].isize() -1; i++ ){
	 if ( hcns[t][i] == hcns[t][i+1] ){
	   hcns[t].erase( hcns[t].begin() + i + 1 );
	   i--;
	 } 
       }


     // initialize indicator values for candidate slide-split regions

     vec< vec<Bool> > split_indicator( hcns.size() );
     for ( size_t t = 0; t < hcns.size(); t++ ){
       split_indicator[t].resize( hcns[t].size() );
       for ( size_t i = 0; i < hcns[t].size(); i++ )
	 split_indicator[t][i] = True;
     }

     // finding spanned hcn intervals

     #pragma omp parallel for
     for ( int t = 0; t < ntigs; t++ ){
       for ( size_t si = 0; si < split_indicator[t].size(); si++ )
	 if ( split_indicator[t][si] ){
	   for ( size_t ci = 0; ci < cov[t].size(); ci++ ){
	     int cstart = cov[t][ci].Start() - trim_back >= 0 ? cov[t][ci].Start() - trim_back : 0;
	     int cstop  = cov[t][ci].Stop() + trim_back > (int)contigs[t].size() ? (int)contigs[t].size() : 
	       cov[t][ci].Stop() + trim_back;
	     ho_interval hint( cstart, cstop );
	     if ( Overlap( hcns[t][si], hint ) >= hcns[t][si].Stop() - hcns[t][si].Start() )
	       split_indicator[t][si] = False;
	   }
	 }
     }
     

     // filtering not spanned intervals
     
     vec< vec<ho_interval> > slideSplits( ntigs );
     for ( int t = 0; t < ntigs; t++ ){
       for ( size_t hi = 0; hi < hcns[t].size(); hi++ ){
	 if ( ! split_indicator[t][hi] ) continue;
	 // require covered region to the left and right of condidate slide split region
	 if ( hcns[t][hi].Stop() >= (int)contigs[t].size() && hcns[t][hi].Start() <= 0 ){
	   split_indicator[t][hi] = False;
	   continue;
	 }
	 Bool goodLeft = False, goodRight = False;
	 for ( size_t cj = 0; cj < covregs[t].size(); cj++ ){
	   if ( covregs[t][cj].Start() < hcns[t][hi].Start() )
	     goodLeft = True;
	   if ( covregs[t][cj].Stop() > hcns[t][hi].Stop() )
	     goodRight = True;
	 }
	 
	 if ( ! goodLeft || ! goodRight )
	   split_indicator[t][hi] = False;
       }
     }
     
     for ( int t = 0; t < ntigs; t++ )
       for ( size_t hi = 0; hi < hcns[t].size(); hi++ )
	 if ( split_indicator[t][hi] )
	   slideSplits[t].push_back( hcns[t][hi] );
     if ( VERBOSE )
       for ( int t = 0; t < ntigs; t++ ){
	 if ( slideSplits[t].size() > 0 ){
	   cout << "slide splits for contig=" << t << ":";
	   for ( int i = 0; i < slideSplits[t].isize() -1; i++ )
	     cout << slideSplits[t][i] << " ";
	   cout << endl;
	 }
       }


     // Find suspicious gaps between contigs.  Not really handling tiny contigs
     // correctly.  This is compensating for a defect in scaffolding.

     log << Date( ) << ": finding suspicious gaps" << endl;
     vec< vec<int> > suspicious_gaps( scaffolds.size( ) );
     #pragma omp parallel for
     for ( int i = 0; i < scaffolds.isize( ); i++ )
     {    const superb& s = scaffolds[i];
          for ( int j = 0; j < s.Ntigs( ); j++ )
          {    int t = s.Tig(j);
               if ( tig_sizes[t] < min_contig ) continue;
               if ( j < s.Ntigs( ) - 1 )
               {    if ( tig_sizes[ s.Tig(j+1) ] < min_contig ) continue;
                    vec<int> to_end;
                    for ( int l = 0; l < cov[t].isize( ); l++ )
                    {    if ( cov[t][l].Stop( ) == tig_sizes[t] )
                              to_end.push_back( cov[t][l].Start( ) );    }
                    Sort(to_end);
                    if ( to_end.empty( ) || to_end.back( ) - to_end.front( ) 
                         < min_spread )
                    {    log << "suspicious gap: " << t << " --> "
                              << s.Tig(j+1) << endl;    
                         suspicious_gaps[i].push_back(j);    }    }
               if ( j > 0 )
               {    if ( tig_sizes[ s.Tig(j-1) ] < min_contig ) continue;
                    vec<int> to_end;
                    for ( int l = 0; l < cov[t].isize( ); l++ )
                    {    if ( cov[t][l].Start( ) == 0 )
                              to_end.push_back( cov[t][l].Stop( ) );    }
                    Sort(to_end);
                    if ( to_end.empty( ) || to_end.back( ) - to_end.front( ) 
                         < min_spread )
                    {    log << "suspicious gap: " << s.Tig(j-1) << " --> "
                              << t << endl;    
                         suspicious_gaps[i].push_back(j-1);    }    }    }    }
     #pragma omp parallel for
     for ( int s = 0; s < scaffolds.isize( ); s++ )
          UniqueSort( suspicious_gaps[s] );

     // Find links that suggest misassambly.

     log << Date( ) << ": finding evil links" << endl;
     vec< vec< vec<blink_fw> > > blink_fwsi(nbatches);
     vec< vec< vec<blink_rc> > > blink_rcsi(nbatches);
     for ( size_t j = 0; j < nbatches; j++ )
     {    blink_fwsi[j].resize(ntigs), blink_rcsi[j].resize(ntigs);    }
     #pragma omp parallel for
     for ( size_t bi = 0; bi < nbatches; bi++ )
     {    size_t start = ( bi * nreads ) / nbatches;
          size_t stop = ( (bi+1) * nreads ) / nbatches;
          for ( size_t id1 = start; id1 < stop; id1++ )
          {    if ( !pairs.isPaired(id1) ) continue;
               int id2 = pairs.getPartnerID(id1);
               if ( aligns0_index[id1] < 0 || aligns0_index[id2] < 0 ) continue;
               alignlet la1 = aligns0[ aligns0_index[id1] ];
               alignlet la2 = aligns0[ aligns0_index[id2] ];
               int t1 = la1.TargetId( ), t2 = la2.TargetId( );
	       ForceAssertLt( t1, ntigs );
	       ForceAssertLt( t2, ntigs );
               if ( to_scaffold[t1] < 0 || to_scaffold[t2] < 0 ) continue;
               int sid1 = to_scaffold[t1], sid2 = to_scaffold[t2];
               if ( sid1 == sid2 ) continue; // not really right
               const superb &s1 = scaffolds[sid1], &s2 = scaffolds[sid2];
               int p1 = to_scaffold_pos[t1], p2 = to_scaffold_pos[t2];
               int read_length2 = la2.Pos2( ) - la2.pos2( );
               int psep = pairs.sep( pairs.getPairID(id1) );
               int pdev = pairs.sd( pairs.getPairID(id1) );
     
               if ( la1.Fw1( ) )
               {    
                    // Ignore links that go off the end of the scaffold.
     
                    if ( la1.Pos2() + psep + int(round(dev_mult*pdev)) + read_length2
                         > s1.SubSuperLength( p1, s1.Ntigs( ) - 1 ) )
                    {    continue;    }
     
                    blink_fw b;
                    b.start1 = la1.pos2( );
                    b.stop1_low 
                         = b.start1 + psep - blink_dev_mult * pdev + read_length2;
                    b.stop1_high 
                         = b.start1 + psep + blink_dev_mult * pdev + read_length2;
                    b.t2 = t2;
                    b.stop2 = la2.Pos2( );
                    blink_fwsi[bi][t1].push_back(b);    }

               else
               {    
                    // Ignore links that go off the end of the scaffold.
     
                    if ( la1.pos2() - psep - int(round(dev_mult*pdev)) - read_length2
                         < 0 ) 
                    {    continue;    }
     
                    blink_rc b;
                    b.stop1 = la1.Pos2( );
                    b.start1_low 
                         = b.stop1 - psep - blink_dev_mult * pdev - read_length2;
                    b.start1_high 
                         = b.stop1 - psep + blink_dev_mult*pdev - read_length2;
                    b.t2 = t2;
                    b.start2 = la2.pos2( );
                    blink_rcsi[bi][t1].push_back(b);    }    }    }
     log << Date( ) << ": combining evil links" << endl;
     vec< vec<blink_fw> > blink_fws(ntigs);
     vec< vec<blink_rc> > blink_rcs(ntigs);
     #pragma omp parallel for
     for ( int t = 0; t < ntigs; t++ )
     {    for ( size_t bi = 0; bi < nbatches; bi++ )
          {    blink_fws[t].append( blink_fwsi[bi][t] );
               blink_rcs[t].append( blink_rcsi[bi][t] );    }    }

     // Display uncovered parts of scaffolds and misassembly evidence associated
     // with them.  Find places where the scaffolds should be broken.

     log << Date( ) << ": displaying uncovered parts" << endl;
     log << "\nuncovered parts of edges:\n";
     vec< vec<ho_interval> > splits(ntigs);
     for ( int i = 0; i < scaffolds.isize( ); i++ )
     {    const superb& s = scaffolds[i];
          if ( s.FullLength( ) < min_scaffold ) continue;
          log << "\nscaffold " << i << "\n";
          for ( int l = 0; l < s.Ntigs( ); l++)
          {    int t = s.Tig(l);
               log << t << "[l=" << tig_sizes[t] << "]:";
               for ( int j = 0; j < un[t].isize( ); j++ )
               {    const ho_interval& u = un[t][j];
                    if ( u.Stop( ) + s.SubSuperLength( 0, l ) - s.Len(l) 
                         < min_dist_from_end )
                    {    continue;    }
                    if ( -u.Start( ) + s.SubSuperLength( l, s.Ntigs( ) - 1 ) 
                         < min_dist_from_end )
                    {    continue;    }
                    log << " " << u;    

                    // Search for supporting links.
     
                    log << " {";
                    vec<blink_fw> fws;
                    for ( int l = 0; l < blink_fws[t].isize( ); l++ )
                    {    const blink_fw& b = blink_fws[t][l];
                         if ( !( b.start1 < u.Start( ) ) ) continue;
                         if ( !( b.stop1_high > u.Stop( ) ) ) continue;
                         fws.push_back(b);    }
                    Sort(fws);
                    Bool support = False;
                    for ( int l = 0; l < fws.isize( ); l++ )
                    {    int m;
                         for ( m = l + 1; m < fws.isize( ); m++ )
                              if ( fws[m].t2 != fws[l].t2 ) break;
                         int count = 0;
                         for ( int v = l; v < m; v++ )
                         {    if ( v == l || fws[v].stop2 != fws[v-1].stop2 )
                                   ++count;    }
                         if ( count < MIN_REACH_AWAY ) continue;
                         support = True;
                         log << " " << fws[l].t2 << "<";
                         for ( int v = l; v < m; v++ )
                         {    log << "(" << fws[v].start1 << ","
                                   << fws[v].stop2 << ")";    }
                         log << ">";
                         l = m - 1;    }
                    log << " ;";
                    vec<blink_rc> rcs;
                    for ( int l = 0; l < blink_rcs[t].isize( ); l++ )
                    {    const blink_rc& b = blink_rcs[t][l];
                         if ( !( b.stop1 > u.Stop( ) ) ) continue;
                         if ( !( b.start1_low < u.Start( ) ) ) continue;
                         rcs.push_back(b);    }
                    Sort(rcs);
                    for ( int l = 0; l < rcs.isize( ); l++ )
                    {    int m;
                         for ( m = l + 1; m < rcs.isize( ); m++ )
                              if ( rcs[m].t2 != rcs[l].t2 ) break;
                         int count = 0;
                         for ( int v = l; v < m; v++ )
                         {    if ( v == l || rcs[v].start2 != rcs[v-1].start2 )
                                   ++count;    }
                         if ( count < MIN_REACH_AWAY ) continue;
                         support = True;
                         log << " " << rcs[l].t2 << "<";
                         for ( int v = l; v < m; v++ )
                         {    log << "(" << rcs[v].stop1 << ","
                                   << rcs[v].start2 << ")";    }
                         log << ">";
                         l = m - 1;    }
                    log << " }";    
                    if (support) splits[t].push_back(u);    }
               log << "\n";    }    }

     // Find gaps in original supers that need to be cut.

     log << Date( ) << ": breaking scaffold gaps" << endl;
     vec< vec<int> > gap_cutters( scaffolds.size( ) );
     
     // Here we would cut at suspicious gaps, but we've turned this
     // code off.
     // gap_cutters = suspicious_gaps;
     
     // Breaking a contig causes also the super to be broken.

     {
       vec<int> to_super( ntigs, -1 );
       vec<int> to_superpos( ntigs, -1 );
       for (int ii=0; ii<scaffolds.isize( ); ii++) {
	 for (int jj=0; jj<scaffolds[ii].Ntigs( ); jj++) {
	   to_super[ scaffolds[ii].Tig( jj ) ] = ii;
	   to_superpos[ scaffolds[ii].Tig( jj ) ] = jj;
	 }
       }
 
       for (int ii=0; ii<ntigs; ii++) {
	 if ( splits[ii].size( ) > 0 ) {
	   int super_id = to_super[ii];
	   int super_pos = to_superpos[ii];
	   if ( super_pos < scaffolds[super_id].Ntigs( ) - 1 )
	     gap_cutters[super_id].push_back( super_pos );
	 }
       }
       
       #pragma omp parallel for
       for (int ii=0; ii<scaffolds.isize( ); ii++)
	 UniqueSort( gap_cutters[ii] );
     }
     
     log << "\nBROKEN GAPS SUMMARY:\n" << endl;

     vec<superb> new_scaffolds;
     for (int old_sid=0; old_sid<scaffolds.isize( ); old_sid++) {
       if ( gap_cutters[old_sid].size( ) < 1 ) {
	 new_scaffolds.push_back( scaffolds[old_sid] );
	 continue;
       }
       
       log << "breaking super_" << old_sid << " ("
	   << scaffolds[old_sid].Ntigs( )
	   << " contigs): ";

       for (int chunk=0; chunk<gap_cutters[old_sid].isize( )+1; chunk++) {
	 int begin_to_use 
	   = chunk == 0 
	   ? 0 
	   : gap_cutters[old_sid][chunk-1] + 1;
	 int end_to_use 
	   = chunk == gap_cutters[old_sid].isize( ) 
	   ? scaffolds[old_sid].Ntigs( )
	   : gap_cutters[old_sid][chunk] + 1;

	 vec<int> to_use;
	 to_use.reserve( end_to_use - begin_to_use );
	 for (int ii=begin_to_use; ii<end_to_use; ii++)
	   to_use.push_back( ii );

	 log << "[" << begin_to_use << "-" << end_to_use << ")";
	 new_scaffolds.push_back( scaffolds[old_sid].SubSuper( to_use ) );
       }
       log << "\n";
       
     }
     log << endl;
     
     // Break contigs (this may add new singleton supers to new_scaffolds).

     log << Date( ) << ": breaking contigs" << endl;
     //ContigsManager manager( contigs, aligns0, aligns0_index, ualigns0, ualigns0_index );
     ContigsManager manager( contigs, aligns0, aligns0_index_unfilt, ualigns0, ualigns0_index );
     //update the fitered index if both unfiltered and filtered are provided
     if ( aligns0_index_unfilt != aligns0_index )
     {
       size_t n_events = 0;
       for (size_t ii=0; ii<aligns0_index_unfilt.size( ); ii++) {
	 if ( aligns0_index_unfilt[ii] >= 0 && aligns0_index[ii] >=0 ) {
	   ForceAssertEq(aligns0_index_unfilt[ii],  aligns0_index[ii]);
	 }
	 if ( aligns0_index_unfilt[ii] < 0 && aligns0_index[ii] >=0 ) {
	   aligns0_index[ii] = aligns0_index_unfilt[ii];
	   n_events++;
	 }
       }
       cout << "n_events= " << n_events << endl;
     }

     log << "BROKEN CONTIGS SUMMARY:\n" << endl;

     vec<int> to_super( ntigs, -1 );
     vec<int> to_pos( ntigs, -1 );
     for (int ii=0; ii<new_scaffolds.isize( ); ii++) {
       for (int jj=0; jj<new_scaffolds[ii].Ntigs( ); jj++) {
	 to_super[ new_scaffolds[ii].Tig( jj ) ] = ii;
	 to_pos[ new_scaffolds[ii].Tig( jj ) ] = jj;
       }
     }
     
     for (int contig_id=0; contig_id<ntigs; contig_id++) {
       if ( splits[contig_id].size( ) < 1 && slideSplits[contig_id].size() < 1 ) continue;

       vec<ho_interval> to_split = splits[contig_id];
       vec<ho_interval> to_slideSplit = slideSplits[contig_id];
       vec<ho_interval> to_split_all = to_split; to_split_all.append( to_slideSplit );
       vec<Bool> is_slide( to_split_all.size(), False );
       for ( size_t i = to_split.size(); i < to_split_all.size(); i++ )
	 is_slide[i] = True;

       SortSync( to_split_all, is_slide );
       for ( size_t i = 0; i < to_split_all.size() -1; i++ )
	 for ( size_t j = i + 1; j < to_split_all.size(); j++ )
	   if ( Overlap( to_split_all[i], to_split_all[j] ) > 0 ){
	     to_split_all.erase( to_split_all.begin() + j);
	     is_slide.erase( is_slide.begin() + j);
	     j--;
	   }

       ForceAssertEq( to_split_all.size(), is_slide.size() );
       
       while( to_split_all.size( ) > 0 ) {
	 ForceAssertEq( to_split_all.size(), is_slide.size() );
	 int cglen_orig = manager.Contigs( )[contig_id].size( );
	 ho_interval win = to_split_all.back( );
	 //ForceAssertLe( win.Stop( ), cglen_orig );
	 //This shouldn't happen.
	 if (win.Stop() > cglen_orig) {
	   log << "ERROR: ";
	   PRINT5(contig_id, cglen_orig, win.Start(), win.Stop(), ToStringBool( is_slide.back() ) );
	   ForceAssertLe( win.Stop( ), cglen_orig );
	 }

	 if ( ! is_slide.back() ){
	   if ( win.Stop( ) < cglen_orig ) {
	     int pos = win.Stop( );
	     log << "breaking tig " << contig_id << " at pos " << pos << "\n";
	     vec<size_t> new_ids = manager.SplitContig( contig_id, pos );
	     int new_contig = new_ids[1];
	     int new_cglen = manager.Contigs( ).back( ).size( );
	     superb newsup;
	     newsup.PlaceFirstTig( new_contig, new_cglen );
	     new_scaffolds.push_back( newsup );
	   }
	   log << "removing tail [" << win.Start( )
	       << "," << win.Stop( )
	       << ") from contig " << contig_id
	       << "\n";
	   manager.CutTail( contig_id, win.Start( ) );
	   int chopped = cglen_orig - win.Start( );
	   int sid = to_super[contig_id];
	   int cpos = to_pos[contig_id];
	   new_scaffolds[sid].SetLen( cpos, cglen_orig - chopped );
	   if ( cpos < new_scaffolds[sid].Ntigs( ) - 1 ) {
	     int gaplen_orig = new_scaffolds[sid].Gap( cpos );
	     new_scaffolds[sid].SetGap( cpos, gaplen_orig + chopped );
	   }
	 }else{
	   int pos1 = win.Start(), pos2 = win.Stop();
	   log << "slide-breaking tig " << contig_id << " at interval (" << pos1 << "," << pos2 << ")\n";
	   vec<size_t> new_ids = manager.SlideSplitContig( contig_id, pos1, pos2 );
	   int new_contig = new_ids[1];
	   int new_cglen = manager.Contigs( ).back( ).size( );
	   superb newsup;
	   newsup.PlaceFirstTig( new_contig, new_cglen );
	   new_scaffolds.push_back( newsup );
	   
	   int chopped = cglen_orig - pos2;
	   int sid = to_super[contig_id];
	   int cpos = to_pos[contig_id];
	   new_scaffolds[sid].SetLen( cpos, cglen_orig - chopped );
	   if ( cpos < new_scaffolds[sid].Ntigs( ) - 1 ) {
	     int gaplen_orig = new_scaffolds[sid].Gap( cpos );
	     new_scaffolds[sid].SetGap( cpos, gaplen_orig + chopped );
	   }
	   
	 }
	 to_split_all.resize( to_split_all.size( ) - 1 );
	 is_slide.resize( is_slide.size( ) -1 );
       }
       
     }
     log << endl;

     swap( scaffolds, new_scaffolds );

     // Cut large gaps.
     
     log << Date( ) << ": cutting large gaps" << endl;

     int n_cuts = 0;
     vec<superb> scaffolds_cut;
     for (int super_id=0; super_id<(int)scaffolds.size( ); super_id++) {
       const superb &sup = scaffolds[super_id];
       int begin = 0;
       while ( begin < sup.Ntigs( ) ) {
	 int end = -1;
	 for (int pos=begin; pos<sup.Ntigs( ); pos++) {
	   if ( pos == sup.Ntigs( ) - 1 ) {
	     end = sup.Ntigs( );
	     break;
	   }
	   int gap = sup.Gap( pos );
	   if ( gap < - max_gap || gap > max_gap ) {
	     log << "  cutting s" << super_id << " between "
		 << pos << " (c" << sup.Tig( pos ) << ") and "
		 << pos+1 << " (c" << sup.Tig( pos+1 ) << "): gap = "
		 << gap << " +/- " << sup.Dev( pos ) << "\n";
	     n_cuts++;
	     end = pos+1;
	     break;
	   }
	 }
	 ForceAssert( end > -1 );
	 superb newsuper;
	 newsuper.SetNtigs( end - begin );
	 for (int pos=begin; pos<end; pos++) {
	   int newpos = pos - begin;
	   newsuper.SetTig( newpos, sup.Tig( pos ) );
	   newsuper.SetLen( newpos, sup.Len( pos ) );
	   if ( pos < end - 1 ) {
	     newsuper.SetGap( newpos, sup.Gap( pos ) );
	     newsuper.SetDev( newpos, sup.Dev( pos ) );
	   }
	 }
	 scaffolds_cut.push_back( newsuper );
	 begin = end;
       }
     }
     if ( n_cuts < 1 )
       log << "no large gaps found.\n" << endl;
     else {
       log << n_cuts << " large gaps found (and cut)\n" << endl;
       swap( scaffolds, scaffolds_cut );
     }
     
     if ( VERBOSE ) {
       log << "\nStats after fixing:\n" << endl;
       ReportScaffoldsN50( scaffolds, log );
     }


     //check alignments
     //check alignments
     for ( size_t i = 0; i < ualigns0_index.size(); i++ )
       for ( size_t j = 0; j < ualigns0_index[i].size(); j++ )
	 if ( ualigns0_index[i][j] >= 0 ){
	   size_t ai = ualigns0_index[i][j];
	   ForceAssertLe( ualigns0[ ai ].Pos2(), (int)contigs[ ualigns0[ ai ].TargetId() ].size() );
	   ForceAssertGe( ualigns0[ ai ].pos2(), 0 );
	 }

     // Done.
     
     log << Date() << ": done" << endl;
     
}
