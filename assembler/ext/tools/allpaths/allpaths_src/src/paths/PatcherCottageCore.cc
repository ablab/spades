//////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "Intvector.h"
#include "CoreTools.h"
#include "ParseSet.h"
#include "Qualvector.h"
#include "VecUtilities.h"
#include "lookup/LookAlign.h"
#include "math/Functions.h"
#include "paths/GetNexts.h"
#include "paths/PatcherCottageCore.h"
#include "paths/UnipathFixerTools.h"
#include "random/Random.h"
#include "random/SequenceEntropy.h"
#include "system/SharedMem.h"
// MakeDepend: cflags OMP_FLAGS

void PatcherCottageCore( basevector L, basevector R, const int sep,
     const int dev, vecbasevector& reads, vecqualvector& quals,
     vec< pair<int,int> >& pairs, String& report,      
     const PatcherCottage_LogParams& log_params, const String& data_dir, 
     vec<int>& start_stop, const int K, const int MAX_READS, const int MAX_JOINS,
     const int MIN_OVERLAP_END, const String& ROOT, size_t p,
     const int u1, const int u2, const Bool LR_special )
{
     double clock = WallClockTime( );

     // Define alignment procedure.

     String temp_filep = ROOT + ToString(p) + ".cottage.fastb";
     String QLT = "QueryLookupTable" + ARG(K, 12) + ARG(MM, 12) + ARG(MC, 0.15)
          + ARG(NH, True) + ARG(VISUAL, True) + ARG(SMITH_WAT, True) 
          + ARG(QUIET, True) + ARG(L, data_dir + "/../genome.lookup")
          + ARG(SEQS, temp_filep);

     // Define heuristic constants.

          const int Q = 15;
          int QP = Q + 1;
          const int seed = 12;
          const int min_align_length = 40;
          const double max_align_errors = 0.2;
          const int end_len = 200;
          const int min_score = 15;
          const double score_prox = 2.0;

          // Initialize outputs.

          vec<basevector> bridge, L_bridge, bridge_R, L_bridge_R, L_bridge_R_abbr;
          start_stop.clear( );
          vecbasevector reads_new;
          vecqualvector quals_new;
          vec<xalign> xaligns0;
          vec<int> xaligns0_index;
          vec<xalign> baligns;
          vec<int> baligns_index;
          vec<Bool> reads_rc;
          int min_qdiff;
          int nreads = reads.size( ), naligns = 0;
          if ( nreads > MAX_READS )
          {    reads.clear( );
               return;    }

          // Print aligned reads.

          if ( log_params.log_align_reads 
               && IsRegularFile( data_dir + "/../genome.lookup" ) )
          {    report += "ALIGNMENT OF READS\n";
               reads.WriteAll(temp_filep);
               String temp_fileq = ROOT + ToString(p) + ".cottage.qualb";
               quals.WriteAll(temp_fileq);
               report += AllOfOutput1( QLT + ARG(QUALS, temp_fileq) );
               Remove(temp_filep), Remove(temp_fileq);    }

          min_qdiff = 0;
          CombinePairs( reads, quals, pairs, score_prox, min_score, min_qdiff,
               log_params.verbosity, report );

          // Scope so alignments don't stay in memory.

          reads_new = reads;
          quals_new = quals;
          {    
               // Build alignments between the reads.

               vec<xalign> xaligns;
               xalign::reads = &reads;
               BuildAlignments( reads, seed, 0, max_align_errors, xaligns0 );
               naligns = xaligns0.size( );

               // Correct errors in the reads.  First we correct the errors, then
               // we identify and delete any alignments that have 3 or more errors,
               // then we error correct again.  The goal of this procedure is to 
               // reduce the number of false alignments that are used in error 
               // correction.

               CorrectErrors( reads, quals, xaligns0, min_align_length, reads_new, 
                    quals_new, log_params.log_correct, report );
               xalign::reads = &reads_new;
               const int max_errors_corr = 3;
               vec<Bool> xdel( xaligns0.size( ), False );
               for ( int i = 0; i < xaligns0.isize( ); i++ )
               {    const xalign& a = xaligns0[i];
                    if ( a.Errs( ) > max_errors_corr ) xdel[i] = True;    }
               EraseIf( xaligns0, xdel );
               xalign::reads = &reads;
               CorrectErrors( reads, quals, xaligns0, min_align_length, reads_new, 
                    quals_new, log_params.log_correct, report );

               xalign::reads = &reads_new;
               if ( log_params.verbosity >= 2 )
               {    report += "\nREADS AFTER CORRECTION:\n";
                    {    ostringstream rout;
                         for ( int i = 0; i < nreads; i++ )
                              reads_new[i].Print( rout, i );
                         report += rout.str( );    }
                    report += "ALIGNMENTS OF ERROR CORRECTED READS\n";
                    {    ostringstream xout;
                         for ( int i = 0; i < xaligns0.isize( ); i++ )
                         {    const xalign& a = xaligns0[i];
                              static int last_id1(-1);
                              if ( a.id1 != last_id1 ) xout << "\n";
                              last_id1 = a.id1;
                              xout << "id1 = " << a.id1 << ", id2 = " << a.id2;
                              xout << ", pos = " << a.pos1 << "/" << a.pos2
                                   << ", ext = " << a.LeftExt1( ) << "/"
                                   << a.RightExt2( ) << ", over+ = " 
                                   << a.OverlapPlus( ) << "\n";    }
                         report += xout.str( );    }    }
               xaligns0_index.resize(nreads+1);
               {    size_t pos = 0;
                    for ( int id = 0; id <= nreads; id++ )
                    {    while( pos < xaligns0.size( ) && xaligns0[pos].id1 < id ) 
                              ++pos;
                         xaligns0_index[id] = pos;    }    }    }

          {

          double exp_clock = WallClockTime( ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

          vecbasevector reads_new_untrimmed(reads_new);

          // Add two pseudo-reads: the last 200 bases of L, and the first 200 bases
          // of R, or less, if L or R are shorter than 200.

          if (LR_special)
          {    L = reads_new[0], R = reads_new[1];    }
          int Ltail = Min( end_len, L.isize( ) ); 
          int Rhead = Min( end_len, R.isize( ) );
          reads_new.push_back( basevector( L, L.isize( ) - Ltail, Ltail ) );
          reads_new.push_back( basevector( R, 0, Rhead ) );
     
          // Find all 16-mers S in the reads for which there is a path from the 
          // first 16-mer of the end of L to the last 16-mer of the beginning of R.
          // Then trim reads from the left and the right until all their 16-mers are
          // in S.

          int FUDGE = 150; // TO REMOVE LATER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          int max_path = sep + 4*dev + Ltail + Rhead - QP + 1 + FUDGE;
          vec<basevector> mers;
          for ( int i = 0; i < nreads+2; i++ )
          {    for ( int j = 0; j <= reads_new[i].isize( ) - QP; j++ )
                    mers.push( reads_new[i], j, QP );    }
          UniqueSort(mers);
          vec< vec<int> > nexts, backs;
          GetNexts( QP, mers, nexts );
          int M = mers.size( );
          for ( int i = 0; i < M; i++ )
               mers[i].ReverseComplement( );
          GetNexts( QP, mers, backs );
          for ( int i = 0; i < M; i++ )
               mers[i].ReverseComplement( );
          vec<int> starts, stops;
          if ( reads_new[nreads].isize( ) >= QP )
          {    starts.push_back( BinPosition( 
                    mers, basevector( reads_new[nreads], 0, QP ) ) );    }
          if ( Rhead >= QP )
          {    stops.push_back( BinPosition( 
                    mers, basevector( reads_new[nreads+1], Rhead - QP, QP ) ) );    }
          vec<int> from_starts( M, False ), to_stops( M, False );
          while( starts.nonempty( ) )
          {    int s = starts.back( );
               starts.pop_back( );
               if ( !from_starts[s] )
               {    from_starts[s] = True;
                    for ( int j = 0; j < nexts[s].isize( ); j++ )
                         starts.push_back( nexts[s][j] );    }    }
          while( stops.nonempty( ) )
          {    int s = stops.back( );
               stops.pop_back( );
               if ( !to_stops[s] )
               {    to_stops[s] = True;
                    for ( int j = 0; j < backs[s].isize( ); j++ )
                         stops.push_back( backs[s][j] );    }    }
          vec<basevector> betweenmers;
          for ( int i = 0; i < M; i++ )
               if ( from_starts[i] && to_stops[i] ) betweenmers.push_back( mers[i] );
          vec<int> left_trim( nreads, 0 ), right_trim( nreads, 0 );
          for ( size_t i = 0; i < reads_new.size( ); i++ )
          {    int j1, j2;
               for ( j1 = 0; j1 <= reads_new[i].isize( ) - QP; j1++ )
               {    basevector b( reads_new[i], j1, QP );
                    if ( BinMember( betweenmers, b ) ) break;    }
               for ( j2 = reads_new[i].isize( ) - QP; j2 >= 0; j2-- )
               {    basevector b( reads_new[i], j2, QP );
                    if ( BinMember( betweenmers, b ) ) break;    }
               if ( j2 >= j1 )
               {    if ( (int) i < nreads )
                    {    if ( i % 2 == 0 )
                              pairs[i/2].first += reads_new[i].isize( ) - (j2+QP);
                         else pairs[i/2].first += j1;
                         quals_new[i].SetToSubOf( quals_new[i], j1, j2 - j1 + QP ); 
                         left_trim[i] = j1;
                         right_trim[i] = reads[i].isize( ) - (j2+QP);    }
                    reads_new[i].SetToSubOf( reads_new[i], j1, j2 - j1 + QP );    }
               else 
               {    reads_new[i].resize(0);
                    if ( (int) i < nreads ) quals_new[i].resize(0);    }    }

          // Chimera detection, V3.  For each trimmed read C, find all untrimmed
          // reads X that align to the left end of C, and all untrimmed reads that
          // align to the right end of C.  Suppose that:
          // 1. There are no reads that extend reads in X to the right and 
          //    reads in Y to the left, except C.
          // 2. There are reads that extend X to the right, which would reach
          //    the end of C -- or --
          // 3. There are reads that extend Y to the left, which would reach 
          //    the beginning of C.
          // 4. X and Y are both nonempty.
          // Then we call C a chimera.

          xalign::reads = &reads_new_untrimmed;
          vec<xalign> yaligns;
          BuildAlignments( reads_new_untrimmed, Q, Q, 1, yaligns );
          vec<int> yaligns_index(nreads+1);
          {    size_t pos = 0;
               for ( int id = 0; id <= nreads; id++ )
               {    while( pos < yaligns.size( ) && yaligns[pos].id1 < id ) ++pos;
                    yaligns_index[id] = pos;    }    }
          for ( int id1 = 0; id1 < nreads; id1++ )
          {    const basevector& r1 = reads_new[id1];
               int len1 = r1.size( );
               if ( len1 == 0 ) continue;
               vec<int> X, Y, Xright, Yleft, XPos1, Ypos1;
               for ( int j = yaligns_index[id1]; j < yaligns_index[id1+1]; j++ )
               {    const xalign& a = yaligns[j];
                    int id2 = a.id2, pos1 = a.pos1, pos2 = a.pos2;
                    const basevector& r2 = reads_new_untrimmed[id2];
                    int len2 = r2.size( );
                    if ( left_trim[id1] > 0 )
                    {    if ( pos1 == 0 ) pos2 += left_trim[id1];
                         else pos1 -= left_trim[id1];    }
                    int offset = pos1 - pos2;
                    int overlap = IntervalOverlap( 0, len1, offset, offset + len2 );
                    if ( overlap < Q ) continue;
                    if ( offset < 0 )
                    {    int errs = 0;
                         for ( int l = 0; l < Q; l++ )
                              if ( r1[l] != r2[l-offset] ) ++errs;
                         if ( errs == 0 )
                         {    X.push_back(id2);
                              XPos1.push_back( offset + len2 );    }    }
                    if ( offset + len2 - len1 > 0 )
                    {    int errs = 0;
                         for ( int l = 1; l <= Q; l++ )
                              if ( r1[len1-l] != r2[len1-l-offset] ) ++errs;
                         if ( errs == 0 )
                         {    Y.push_back(id2);
                              Ypos1.push_back(offset);    }    }    }
               if ( X.empty( ) || Y.empty( ) ) continue;
               Bool condition2 = False, condition3 = False;
               for ( int i = 0; i < X.isize( ); i++ )
               {    int idx = X[i];
                    for ( int j = xaligns0_index[idx]; 
                         j < xaligns0_index[idx+1]; j++ )
                    {    const xalign& a = xaligns0[j];
                         if ( a.RightExt2( ) > 0 && a.id2 != id1 ) 
                         {    Xright.push_back(a.id2);    
                              if ( XPos1[i] + a.RightExt2( ) >= len1 )
                                   condition2 = True;    }    }    }
               for ( int i = 0; i < Y.isize( ); i++ )
               {    int idy = Y[i];
                    for ( int j = xaligns0_index[idy]; 
                         j < xaligns0_index[idy+1]; j++ )
                    {    const xalign& a = xaligns0[j];
                         if ( a.LeftExt2( ) > 0 && a.id2 != id1 ) 
                         {    Yleft.push_back(a.id2);    
                              if ( Ypos1[i] - a.LeftExt2( ) <= 0 )
                              {    condition3 = True;    }    }    }    }
               if ( condition2 || condition3 )
               {    UniqueSort(Xright), UniqueSort(Yleft);
               vec<int> M = Intersection( Xright, Yleft );
               if ( M.empty( ) ) 
                    {    if ( log_params.verbosity >= 2 )
                         {    report += "newer method, see putative chimera " 
                                   + ToString(id1) + "\n";    }
                         reads_new[id1].resize(0);    
                         quals_new[id1].resize(0);    }    }    }

          // Cap the size of the pseudo-reads to try to avoid junk on their ends.
          // There is probably a better way to get rid of this junk but I haven't
          // figured out what it is.  Also punt if the ends end up being too small.

          const int mpr = 100;
          if ( reads_new[nreads].isize( ) > mpr ) reads_new[nreads].resize(mpr);
          if ( reads_new[nreads+1].isize( ) > mpr ) 
          {    reads_new[nreads+1].SetToSubOf( reads_new[nreads+1],
                    reads_new[nreads+1].isize( ) - mpr, mpr );    }
          if ( reads_new[nreads].isize( ) < mpr 
               || reads_new[nreads+1].isize( ) < mpr )
          {    reads_new[nreads].resize(0);
               reads_new[nreads+1].resize(0);    }

          // Combine pairs and track changes.

          vec<Bool> readsizes1(nreads), readsizes2(nreads);
          for ( int i = 0; i < nreads; i++ )
               readsizes1[i] = reads_new[i].size( );
          min_qdiff = 0;
          CombinePairs( reads_new, quals_new, pairs, score_prox, min_score, 
               min_qdiff, log_params.verbosity, report );
          for ( int i = 0; i < nreads; i++ )
               readsizes2[i] = reads_new[i].size( );

          // Find all perfect overlaps of length 15 or more between the error
          // corrected reads.

          xalign::reads = &reads_new;
          vec<xalign> xaligns;
          BuildAlignments( reads_new, Q, Q, 0, xaligns );
          vec<int> xaligns_index(nreads+3);
          {    size_t pos = 0;
               for ( int id = 0; id <= nreads+2; id++ )
               {    while( pos < xaligns.size( ) && xaligns[pos].id1 < id ) ++pos;
                    xaligns_index[id] = pos;    }    }

          // Identify illegal trimming events.  (PRELIMINARY VERSION)
          // Let X denote a read, Xt its trimmed version.  Suppose that Xt
          // does not overlap Rt.  Suppose that X has been trimmed on the right,
          // and that for some Y, we have a perfect overlap like this
          //
          //             ==========Xt-------->   (-- is trimmed part)
          //                     ======Yt========>
          // 
          // Then we declare the trimming of X to be illegal, and delete the read.

         for ( int id = 0; id < nreads; id++ )
          {    if ( right_trim[id] == 0 ) continue;
               // Don't mess with pairs that have been combined.
               if ( readsizes2[id] != readsizes1[id] ) continue;
               Bool hits_R = False;
               for ( int j = xaligns_index[id]; j < xaligns_index[id+1]; j++ )
               {    if ( xaligns[j].id2 == nreads + 1 )
                    {    hits_R = True;
                         break;    }    }
               if (hits_R) continue;
               xalign::reads = &reads;
               for ( int j = xaligns0_index[id]; j < xaligns0_index[id+1]; j++ )
               {    const xalign& a = xaligns0[j];
                    if ( a.Errs( ) > 0 ) continue;
                    if ( a.Overlap( ) < right_trim[id] + 10 ) continue;
                    if ( a.RightExt2( ) == 0 ) continue;
                    reads_new[id].resize(0);
                    quals_new[id].resize(0);
                    if ( log_params.verbosity >= 2 )
                         report += "invalidating read " + ToString(id) + "\n";
                    left_trim[id] = right_trim[id] = 0;
                    break;    }
               xalign::reads = &reads_new;    }

          // Print trimmed reads.

          if ( log_params.verbosity >= 1 )
          {    int total_bases = 0;
               for ( size_t i = 0; i < reads_new.size( ); i++ )
                    total_bases += reads_new[i].size( );
               report += "\ntotal bases in trimmed reads: " 
                    + ToString(total_bases) + "\n\n";
               report += "\nREADS AFTER TRIMMING (LAST TWO ARE L AND R):\n";
               ostringstream rout;
               for ( int i = 0; i < nreads+2; i++ )
                    reads_new[i].Print( rout, i );
               report += rout.str( );    }

          // Identify and mark those alignments x --> xi that are to be tried as 
          // right extensions.  We apply the following tests:
          //
          // 1. Read xi must extend x to the right by at least one base.
          //
          // 2. If the overlap sequence is <= 20 bases, its dinucleotide entropy
          // must be at least one.
          //
          // 3. The extending read must itself have a right extension.
          // (commented out, and in any case should not be applied to read
          // nreads+1, the beginning of R)

          vec<Bool> has_right( nreads+2, False );
          for ( int i = 0; i < xaligns.isize( ); i++ )
          {    const xalign& a = xaligns[i];
               if ( a.RightExt2( ) > 0 ) has_right[a.id1] = True;    }
          vec<Bool> good_right( xaligns.size( ), False );
          for ( int i = 0; i < xaligns.isize( ); i++ )
          {    const xalign& a = xaligns[i];
               int id1 = a.id1, id2 = a.id2, pos1 = a.pos1, pos2 = a.pos2;
               if ( reads_new[id1].size( ) == 0 ) continue;
               int len1 = reads_new[id1].size( ), len2 = reads_new[id2].size( );
               if ( a.RightExt2( ) <= 0 ) continue;
               int over = Min( len1 - pos1, len2 - pos2 );
               if ( over <= 20 )
               {    basevector b;
                    b.SetToSubOf( reads_new[ pos2 == 0 ? id2 : id1 ], 0, over );
                    if ( DinukeEntropy(b) < 1 ) continue;    }
               good_right[i] = True;    }
          {    if ( log_params.verbosity >= 2 )
               {    ostringstream xout;
                    for ( int i = 0; i < xaligns.isize( ); i++ )
                    {    const xalign& a = xaligns[i];
                         static int last_id1(-1);
                         if ( a.id1 != last_id1 ) xout << "\n";
                         last_id1 = a.id1;
                         xout << "id1 = " << a.id1;
                         if ( a.id1 == nreads ) xout << "[begin]";
                         xout << ", id2 = " << a.id2;
                         if ( a.id2 == nreads+1 ) xout << "[end]";
                         xout << ", pos = " << a.pos1 << "/" << a.pos2 << ", ext = "
                              << a.LeftExt1( ) << "/" << a.RightExt2( )
                              << ", over+ = " << a.OverlapPlus( );
                         if ( good_right[i] ) xout << " (right!)";
                         xout << "\n";    }
                    report += xout.str( );    }    }

          // Start looking for closures.

          const int max_iterations = 5000;
          int iterations = 0;
          align_chain::reads = &reads_new;
          vec<align_chain> closure_chains;
          align_chain c(nreads);
          vec<align_chain> active_chains;
          active_chains.push_back(c);
          ostringstream* pout = 0;
          if ( log_params.verbosity >= 3 ) pout = new ostringstream;
          while( active_chains.nonempty( ) )
          {    if ( ++iterations > max_iterations ) break;
               align_chain c( active_chains.back( ) );
               active_chains.pop_back( );
               if ( log_params.verbosity >= 3 )
               {    *pout << "\nexploring";
                    for ( int j = 0; j < c.ids.isize( ); j++ )
                         *pout << " " << c.ids[j];
                    *pout << endl;
                    c.Seq( ).Print( 
                         *pout, "explore_" + ToString(iterations-1) );    }
               if ( c.ids.back( ) == nreads+1 ) closure_chains.push_back(c);
               else
               {    vec<align_chain> c_next;
                    GetRights( c, xaligns, xaligns_index, good_right, c_next, nreads,
                         MIN_OVERLAP_END );
                    if ( log_params.verbosity >= 3 )
                    {    *pout << "found " << c_next.size( )
                              << " right extensions" << "\n";    }
                    for ( int j = 0; j < c_next.isize( ); j++ )
                         active_chains.push_back(c_next[j]);    }    }
          if ( log_params.verbosity >= 3 )
          {    report += pout->str( );
               delete pout;    }
          if ( log_params.verbosity >= 1 )
          {    report += "iterations = " + ToString(iterations) + "\n";
               if ( iterations > max_iterations )
                    report += "max_iterations hit\n";    }

          // Convert closures into sequences.

          if ( iterations > max_iterations ) closure_chains.clear( );
          if ( log_params.verbosity >= 1 )
          {    report += "\nExpected closure length: " + ToString(sep) + " +/- "
                    + ToString(dev) + ", max = " + ToString(max_path) + "\n";
               report += "\nCLOSURES:\n";    }
          for ( int i = 0; i < closure_chains.isize( ); i++ )
          {    const align_chain& c = closure_chains[i];
               if ( log_params.verbosity >= 1 )
               {    report += "length = " + ToString( c.Len( ) ) + "\n";
                    report += "[" + ToString(i+1) + "]";
                    for ( int j = 0; j < c.ids.isize( ); j++ )
                         report += " " + ToString(c.ids[j]);
                    report += "\n";    }
               basevector b = c.Seq( );
               if ( log_params.verbosity >= 1 )
               {    ostringstream bout;
                    b.Print( bout, i+1 );
                    report += bout.str( );    }
               if ( !Member( bridge, b ) )
               {    bridge.push_back(b);
                    L_bridge.push_back( 
                         Cat( basevector( L, 0, L.isize( ) - Ltail ), b ) );
                    bridge_R.push_back(
                         Cat( b, basevector( R, Rhead, R.isize( ) - Rhead ) ) );
                    if (log_params.log_aligns)
                    {    L_bridge_R.push_back(
                              Cat( basevector( L, 0, L.isize( ) - Ltail ), b, 
                                   basevector( R, Rhead, R.isize( ) - Rhead ) ) );
                         int Lpart = Min( 500, L.isize( ) - Ltail );
                         L_bridge_R_abbr.push_back(
                              Cat( basevector( L, L.isize( )-Ltail-Lpart, Lpart ),
                                   b, basevector( R, Rhead, 
                                   Min( 500, R.isize( ) - Rhead ) ) ) );    }
                    start_stop.push_back( L.isize( ) - Ltail, Rhead );    }    }
          if ( log_params.verbosity >= 1 )
          {    report += "unique closures: " + ToString( bridge.size( ) ) + "\n";
               report += TimeSince(exp_clock) 
                    + " used in experimental section\n";    }
          if ( bridge.isize( ) > MAX_JOINS )
          {    if ( log_params.verbosity >= 1 ) report += "MAX_JOINS hit\n";
               bridge.clear( ), L_bridge.clear( ), bridge_R.clear( );    }

          // Record results.

          reads.clear( );
          for ( int i = 0; i < bridge.isize( ); i++ )
          {    reads.push_back( bridge[i] );
               reads.push_back( L_bridge[i] );
               reads.push_back( bridge_R[i] );    }

          // Print alignments of joins.

          if ( log_params.log_aligns 
               && IsRegularFile( data_dir + "/../genome.lookup" )
               && ( bridge.nonempty( ) || log_params.log_aligns_all ) )
          {    vec<look_align> L_aligns, R_aligns, J_aligns_abbr;
               look_align la;
               report += "\n================================================="
                    + String("===================================\n\n");
     
               // Align L and R.
     
               vecbasevector LR(2);
               LR[0] = L, LR[1] = R;
               for ( int s = 0; s < 2; s++ )
               {    report += String("ALIGNMENT OF ") + 
                         ( s == 0 ? "LEFT" : "RIGHT" ) + " SEQUENCE " 
                         + ToString( s == 0 ? u1 : u2 ) + "\n";
                    LR.WriteOne( temp_filep, ( s == 0 ? 0 : 1 ) );
                    vec<String> LR_out = AllOfOutput( QLT + ARG(PARSEABLE, True) );
                    for ( int j = 0; j < LR_out.isize( ); j++ )
                    {    if ( LR_out[j].Contains( "QUERY", 0 ) )
                         {    la.ReadParseable( LR_out[j] );
                              if ( s == 0 ) L_aligns.push_back(la);
                              else R_aligns.push_back(la);    }
                         else 
                         {    if ( LR_out[j].Contains( "0", 0 ) )
                                   report += ( s == 0 ? "LEFT" : "RIGHT" )
                                        + LR_out[j].After( "0" ) + "\n";
                              else report += LR_out[j] + "\n";    }    }    }

               // Align joins.
     
               vecbasevector BJ;
               vec<String> BJ_out;
               const int MAX_JOINS_TO_PRINT = 10;
               if ( bridge.isize( ) > MAX_JOINS_TO_PRINT )
               {    report += ToString( bridge.size( ) ) 
                         + " JOINS [DETAILS OMITTED]\n";    }
               else
               {    for ( int j = 0; j < bridge.isize( ); j++ )
                         BJ.push_back( L_bridge_R[j] );
                    BJ.WriteAll(temp_filep);
                    BJ_out = AllOfOutput( QLT + ARG(PARSEABLE, True) );
                    int bid = -1;
                    vec< vec<String> > BJ_out_i( bridge.size( ) );
                    for ( int l = 0; l < BJ_out.isize( ); l++ )
                    {    if ( BJ_out[l].Contains( "QUERY", 0 ) )
                         {    la.ReadParseable( BJ_out[l] );
                              bid = la.query_id;    }
                         else if ( bid >= 0 && bid < bridge.isize( ) )
                         {    BJ_out_i[bid].push_back( BJ_out[l] );    }    }
                    for ( int j = 0; j < bridge.isize( ); j++ )
                    {    report += "ALIGNMENT OF JOIN " + ToString(j+1) 
                              + " OF SEQUENCE " + ToString(u1) + " TO " 
                              + ToString(u2) + "\n\n";
                         for ( int l = 0; l < BJ_out_i[j].isize( ); l++ )
                              report += BJ_out_i[j][l] + "\n";    }    }
               BJ.clear( );
               for ( int j = 0; j < bridge.isize( ); j++ )
                    BJ.push_back( L_bridge_R_abbr[j] );
               BJ.WriteAll(temp_filep);
               BJ_out = AllOfOutput( QLT + ARG(PARSEABLE, True) );
               for ( int l = 0; l < BJ_out.isize( ); l++ )
               {    if ( BJ_out[l].Contains( "QUERY", 0 ) )
                    {    la.ReadParseable( BJ_out[l] );
                         J_aligns_abbr.push_back(la);    }    }

               // Analyze alignments.
     
               if ( bridge.nonempty( ) ) 
               {    int min_errs_J = 1000000000;
                    for ( int i = 0; i < J_aligns_abbr.isize( ); i++ )
                         min_errs_J = Min( min_errs_J, J_aligns_abbr[i].Errors( ) );
                    if ( min_errs_J > 0 )
                    {    report += "\nNO PERFECT JOIN FOUND, MIN ERRORS = "
                              + ToString(min_errs_J) + "\n";    }
                    else report += "have perfect join\n";
                    Bool consistent = False;
                    const int max_dist = 1000;
                    for ( int iL = 0; iL < L_aligns.isize( ); iL++ )
                    {    const look_align& aL = L_aligns[iL];
                         for ( int iR = 0; iR < R_aligns.isize( ); iR++ )
                         {    const look_align& aR = R_aligns[iR];
                              if ( aL.TargetId( ) != aR.TargetId( ) ) continue;
                              if ( aL.Fw1( ) && aR.Fw1( ) )
                              {    if ( Abs( aR.pos2( ) - aL.Pos2( ) ) <= max_dist )
                                        consistent = True;    }
                              if ( aL.Rc1( ) && aR.Rc1( ) )
                              {    if ( Abs( aL.pos2( ) - aR.Pos2( ) ) <= max_dist )
                                        consistent = True;    }    }    }
                    if (consistent)
                    {    report += "\nFOUND CONSISTENT INPUTS TO JOIN OF " 
                              + ToString(u1) + " TO " + ToString(u2) 
                              + "\n";    }    }
     
               // Clean up.

               Remove(temp_filep);    }

          // Finalize report.

          if ( log_params.verbosity >= 0.1 )
          {    report += ToString(p) + ": (" + ToString(u1) + "," + ToString(u2) 
                    + "), " + ToString(nreads) + " reads, " + ToString(naligns) 
                    + " aligns, " + ToString( bridge.size( ) ) + " bridges, "
                    + TimeSince(clock) + "\n";    }    }    }
