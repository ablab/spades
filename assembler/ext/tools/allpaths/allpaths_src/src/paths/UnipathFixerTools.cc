///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>
#include <sys/wait.h>

#include "Basevector.h"
#include "FeudalMimic.h"
#include "Intvector.h"
#include "MainTools.h"
#include "PairsManager.h"
#include "ParallelVecUtilities.h"
#include "Qualvector.h"
#include "Set.h"
#include "VecTemplate.h"
#include "VecUtilities.h"
#include "math/Functions.h"
#include "paths/Ulink.h"
#include "paths/UnibaseCopyNumberCommon.h"
#include "paths/UnibaseCopyNumber2Core.h"
#include "paths/UnipathFixerTools.h"
#include "system/SysConf.h"
#include "system/WorklistMP.h"
#include "util/SearchFastb2Core.h"
#include <vector>

size_t const PCottageJoinData::HEADER;

void AlignReadsToUnipaths( const String& run_dir, const String& jump_reads,
     const String& frag_reads, const String& frag_reads_edit, const Bool USE_JUMPS,
     const int MAX_PLACEMENTS, const String& unifile,
     vec< triple<int64_t,int64_t,int> >& ALIGNS,
     vec< triple<int64_t,int64_t,int> >& JALIGNS, const String& checkpoint_head )
{
     // Set up.

     ALIGNS.clear( ), JALIGNS.clear( );
     String temp_file = run_dir + "UnipathFixer.fastb";
     String jump_reads_file = run_dir + "/" + jump_reads + ".fastb";
     String frag_reads_file = run_dir + "/" + frag_reads + ".fastb";
     String frag_reads_edit_file = run_dir + "/" + frag_reads_edit + ".fastb";
     ForceAssert( IsRegularFile(jump_reads_file) );
     ForceAssert( IsRegularFile(frag_reads_file) );
     ForceAssert( IsRegularFile(frag_reads_edit_file) );
     ForceAssertEq( MastervecFileObjectCount(frag_reads_file),
          MastervecFileObjectCount(frag_reads_edit_file) );
     if ( checkpoint_head != "" )
          Echo( "testing", checkpoint_head + ".test" ); // make sure we can write

     // Align jumps (just one pass for now).

     if (USE_JUMPS)
     {    cout << Date( ) << ": " 
               << "-----------------------------------------------------" << endl;
          cout << Date( ) << ": loading jumps" << endl;
          {    vecbasevector jreads( run_dir + "/" + jump_reads + ".fastb" );
               for ( size_t i = 0; i < jreads.size( ); i++ )
                    jreads[i].ReverseComplement( );
               cout << Date( ) << ": total jump reads = " 
                    << ToStringAddCommas( jreads.size( ) ) << endl;
               for ( size_t i = 0; i < jreads.size( ); i++ )
                    jreads[i].resize(20);
               jreads.WriteAll(temp_file);    }
          vec< triple<int64_t,int64_t,int> > jaligns;
          SearchFastb2( temp_file, unifile, 20, &jaligns, 0, MAX_PLACEMENTS );
          for ( size_t i = 0; i < jaligns.size( ); i++ )
               if ( jaligns[i].third >= 0 ) JALIGNS.push_back( jaligns[i] );
          if ( checkpoint_head != "" )
               BinaryWrite3( checkpoint_head + ".JALIGNS", JALIGNS );    }

     // Count reads.

     size_t nreads = MastervecFileObjectCount(frag_reads_file);
     cout << Date( ) << ": " 
          << "-----------------------------------------------------" << endl;
     cout << Date( ) << ": total fragment reads = " 
          << ToStringAddCommas(nreads) << endl;

     // Set up data structures to track alignments.

     vec< triple<int64_t,int64_t,int> > aligns;
     vec<Bool> aligned( nreads, False );

     // Look for perfect placements of fragment reads.

     cout << Date( ) << ": trying full length alignment" << endl;
     SearchFastb2( frag_reads_file, unifile, 80, &aligns, 0, MAX_PLACEMENTS );
     for ( size_t i = 0; i < aligns.size( ); i++ )
     {    aligned[ aligns[i].first ] = True;
          if ( aligns[i].third >= 0 ) ALIGNS.push_back( aligns[i] );    }
     cout << Date( ) << ": ";
     PRINT( Sum(aligned) );
     if ( checkpoint_head != "" )
     {    BinaryWrite3( checkpoint_head + ".ALIGNS1", ALIGNS );
          BinaryWrite3( checkpoint_head + ".aligned1", aligned );    }

     // Look for perfect placements of truncated fragment reads.

     for ( int pass = 1; pass <= 4; pass++ )
     {    int L = 0, KK = 0;
          if ( pass == 1 ) { L = 80, KK = 80; }
          if ( pass == 2 ) { L = 60, KK = 60; }
          if ( pass == 3 ) { L = 40, KK = 40; }
          if ( pass == 4 ) { L = 20, KK = 20; }
          cout << Date( ) << ": " 
               << "-----------------------------------------------------" << endl;
          cout << Date( ) << ": trying truncation to " << L << endl;
          {    vecbasevector reads(frag_reads_file);
               for ( size_t i = 0; i < nreads; i++ )
               {    if ( aligned[i] ) reads[i].resize(0);
                    if ( reads[i].isize( ) > L ) reads[i].resize(L);    }
               reads.WriteAll(temp_file);    }
          SearchFastb2( temp_file, unifile, KK, &aligns, 0, MAX_PLACEMENTS );
          for ( size_t i = 0; i < aligns.size( ); i++ )
          {    aligned[ aligns[i].first ] = True;
               if ( aligns[i].third >= 0 ) ALIGNS.push_back( aligns[i] );    }
          cout << Date( ) << ": ";
          PRINT( Sum(aligned) );
          if ( checkpoint_head != "" )
          {    BinaryWrite3( checkpoint_head + ".ALIGNS" + ToString(pass+1), 
                    ALIGNS );
               BinaryWrite3( checkpoint_head + ".aligned" + ToString(pass+1), 
                    aligned );    }    }

     // Try using the second 20 bases of the fragment reads..

     cout << Date( ) << ": " 
          << "-----------------------------------------------------" << endl;
     cout << Date( ) << ": trying the second 20 bases" << endl;
     int start = 20, len = 20;
     {    vecbasevector reads(frag_reads_file);
          for ( size_t i = 0; i < nreads; i++ )
          {    if ( aligned[i] || reads[i].isize() < start+len )
                 reads[i].resize(0);
               else
                 reads[i].SetToSubOf( reads[i], start, len );    }
          reads.WriteAll(temp_file);    }
     SearchFastb2( temp_file, unifile, 20, &aligns, 0, MAX_PLACEMENTS );
     for ( size_t i = 0; i < aligns.size( ); i++ )
     {    if ( aligns[i].third >= 0 )
          {    aligns[i].third -= start;
               aligned[ aligns[i].first ] = True;
               ALIGNS.push_back( aligns[i] );    }    }
     cout << Date( ) << ": ";
     PRINT( Sum(aligned) );
     if ( checkpoint_head != "" )
     {    BinaryWrite3( checkpoint_head + ".ALIGNS6", ALIGNS );
          BinaryWrite3( checkpoint_head + ".aligned6", aligned );    }

     // Try using the edited fragment reads.

     int L = 20, KK = 20;
     cout << Date( ) << ": " 
          << "-----------------------------------------------------" << endl;
     cout << Date( ) << ": trying truncation of edited reads to " 
          << L << endl;
     {    vecbasevector reads(frag_reads_edit_file);
          for ( size_t i = 0; i < nreads; i++ )
          {    if ( aligned[i] ) reads[i].resize(0);
               else if ( reads[i].isize( ) > L ) reads[i].resize(L);    }
          reads.WriteAll(temp_file);    }
     SearchFastb2( temp_file, unifile, KK, &aligns, 0, MAX_PLACEMENTS );
     for ( size_t i = 0; i < aligns.size( ); i++ )
     {    aligned[ aligns[i].first ] = True;
          if ( aligns[i].third >= 0 ) ALIGNS.push_back( aligns[i] );    }
     cout << Date( ) << ": ";
     PRINT( Sum(aligned) );
     if ( checkpoint_head != "" )
     {    BinaryWrite3( checkpoint_head + ".ALIGNS7", ALIGNS );
          BinaryWrite3( checkpoint_head + ".aligned7", aligned );    }
     cout << Date( ) << ": "
          << "-----------------------------------------------------" << endl;
     Remove(temp_file);    }

void ExtendAligns( 
     // INPUTS:
     const int K, const vec< triple<int64_t,int64_t,int> >& ALIGNS,
     const vecbasevector& reads, const vecqualvector& quals, 
     const vecbasevector& unibases, const vec< vec<int> > & nexts, 
     // LOGGING:
     const Bool print_alignments, const Bool print_segments, 
     // OUTPUTS:
     VecIntVec* calls, VecIntVec* qcalls, int64_t& total, int64_t& qtotal, 
     int64_t& qgood, vec<segalign>& SEGS, 
     // HEURISTICS:
     const Bool off_end_ok, const Bool get_calls )
{   
     SEGS.clear( );

     // Define batches for parallel run.

     size_t batch_size 
          = Max( (size_t) 1, ALIGNS.size( ) / ( 10 * omp_get_max_threads( ) ) );
     vec<size_t> batch_starts;
     batch_starts.push_back(0);
     while( batch_starts.back( ) < ALIGNS.size( ) )
     {    size_t b = batch_starts.back( ) + batch_size - 1;
          while( b < ALIGNS.size( ) - 1 && ALIGNS[b].first == ALIGNS[b+1].first )
               ++b;
          batch_starts.push_back( Min( ALIGNS.size( ), b+1 ) );    }
     int nbatches = batch_starts.size( ) - 1;

     // Define outputs of computations to be merged after parallel loop.

     vec< vec<segalign> > SEGS_b(nbatches);

     // Go through the batches.

     #pragma omp parallel for
     for ( size_t bb = 0; bb < batch_starts.size( ) - 1; bb++ )
     {    for ( size_t i = batch_starts[bb]; i < batch_starts[bb+1]; i++ )
          {    size_t j;
               for ( j = i + 1; j < ALIGNS.size( ); j++ )
                    if ( ALIGNS[j].first != ALIGNS[i].first ) break;
               vec< vec<int> > ALL_UNIS;
               vec<int> ALL_ERRS, ALL_UPOS, ALL_RPOS;
               size_t rid = ALIGNS[i].first;
               for ( size_t k = i; k < j; k++ )
               {    const triple<int64_t,int64_t,int>& A = ALIGNS[k];
                    int64_t uid = A.second, pos = A.third, pos0 = pos;
                    vec< vec<int> > UNIS, XUNIS;
                    vec<int> ERRS, RPOS, XERRS, XRPOS;
                    vec<int> p0;
                    p0.push_back(uid);
     
                    // Record info for first unipath that the read aligns to.
     
                    int errs = 0, l;
                    for ( l = 0; l < reads[rid].isize( ); l++ )
                    {    if ( pos >= unibases[uid].isize( ) ) break;
                         if ( pos >= 0 && reads[rid][l] != unibases[uid][pos] ) 
                              ++errs;
                         pos++;    }
                    UNIS.push_back(p0), ERRS.push_back(errs), RPOS.push_back(l);
     
                    // Extend across unipaths that follow it.
     
                    while( UNIS.nonempty( ) )
                    {    vec<int> p = UNIS.back( );
                         int errs = ERRS.back( );
                         int rpos = RPOS.back( );
                         UNIS.pop_back( ), ERRS.pop_back( ), RPOS.pop_back( );
                         int uid = p.back( );
                         if ( rpos == reads[rid].isize( ) || 
                              ( off_end_ok && nexts[uid].empty( ) ) )
                         {    XUNIS.push_back(p);
                              XERRS.push_back(errs);
                              XRPOS.push_back(rpos);
                              continue;    }
                         for ( int z = 0; z < nexts[uid].isize( ); z++ )
                         {    int un = nexts[uid][z], upos = K - 1;
                              int errsx = errs, l;
                              for ( l = rpos; l < reads[rid].isize( ); l++ )
                              {    if ( upos >= unibases[un].isize( ) ) break;
                                   if ( reads[rid][l] != unibases[un][upos] ) 
                                        ++errsx;
                                   upos++;    }
                              vec<int> p0(p);
                              p0.push_back(un);
                              UNIS.push_back(p0);
                              ERRS.push_back(errsx);
                              RPOS.push_back(l);    }    }    
                    for ( int w = 0; w < XUNIS.isize( ); w++ )
                    {    ALL_UNIS.push_back( XUNIS[w] );
                         ALL_ERRS.push_back( XERRS[w] );
                         ALL_UPOS.push_back(pos0);    
                         ALL_RPOS.push_back( XRPOS[w] );    }    }

               // Roughly, keep only the alignments that have the minimum number of
               // errors amongst alignments of given length.

               /*
               vec<Bool> to_delete( ALL_UNIS.size( ), False );
               ReverseSortSync( ALL_RPOS, ALL_ERRS, ALL_UNIS, ALL_UPOS );
               int best_errs = 1000000000;
               for ( size_t l1 = 0; l1 < ALL_RPOS.size( ); l1++ )
               {    size_t l2;
                    for ( l2 = l1 + 1; l2 < ALL_RPOS.size( ); l2++ )
                         if ( ALL_RPOS[l2] != ALL_RPOS[l1] ) break;
                    int min_errs = 1000000000;
                    for ( size_t j = l1; j < l2; j++ )
                         min_errs = Min( min_errs, ALL_ERRS[j] );
                    for ( size_t j = l1; j < l2; j++ )
                    {    if ( ALL_ERRS[j] > Min( min_errs, best_errs ) )
                              to_delete[j] = True;    }
                    best_errs = Min( best_errs, min_errs );
                    l1 = l2 - 1;    }
               EraseIf( ALL_UNIS, to_delete ), EraseIf( ALL_ERRS, to_delete );
               EraseIf( ALL_UPOS, to_delete ), EraseIf( ALL_RPOS, to_delete );
               */

               // Define segments.  First member is start position on read.  Second 
               // member is unibase id.  Third member is start position on unibase.

               vec< triple<int,int,int> > segs;

               // Build segments from the alignments.

               if ( ALL_ERRS.empty( ) )
               {    i = j - 1;
                    continue;    }
               int max_rpos = Max(ALL_RPOS), min_errs = 1000000000;
               for ( size_t w = 0; w < ALL_ERRS.size( ); w++ )
               {    if ( ALL_RPOS[w] == max_rpos )
                         min_errs = Min( min_errs, ALL_ERRS[w] );    }
               for ( size_t w = 0; w < ALL_UNIS.size( ); w++ )
               {    if ( ALL_ERRS[w] > min_errs ) continue;
                    if (print_alignments)
                    {    cout << "read " << rid << ", start pos = " << ALL_UPOS[w]
                              << ", unipaths = [";
                         for ( int z = 0; z < ALL_UNIS[w].isize( ); z++ )
                         {    if ( z > 0 ) cout << ",";
                              cout << ALL_UNIS[w][z];    }
                         cout << "], errs = " << ALL_ERRS[w] << endl;    }
                    int rpos = 0, upos = ALL_UPOS[w];
                    for ( int z = 0; z < ALL_UNIS[w].isize( ); z++ )
                    {    int u = ALL_UNIS[w][z];
                         {    int rposx(rpos), uposx(upos);
                              if ( rposx < 0 ) 
                              {    uposx -= rposx;
                                   rposx = 0;    }
                              if (print_segments)
                              {    cout << "read " << rid << ", starting at " 
                                        << rposx << " maps to unipath " << u
                                        << " starting at " << uposx << endl;    }
                              segs.push( rposx, u, uposx );    }
                         rpos += unibases[u].isize( ) - upos - ( K - 1 );
                         upos = 0;    }    }
               UniqueSort(segs);

               // Look for further unibase extensions on end.

               for ( size_t w = 0; w < segs.size( ); w++ )
               {    int rpos = segs[w].first, u = segs[w].second; 
                    int upos = segs[w].third;
                    int dist_to_end = // gap between read end and unibase end
                         unibases[u].isize( ) - ( reads[rid].isize() - rpos + upos );
                    int rposn = reads[rid].isize( ) + dist_to_end - (K-1);
                    if ( dist_to_end >= 0 && dist_to_end <= K-1 && rposn >= 0 )
                    {    for ( int q = 0; q < nexts[u].isize( ); q++ )
                         {    int un = nexts[u][q], uposn = 0;
                              triple<int,int,int> seg( rposn, un, uposn );
                              if ( !Member( segs, seg ) )
                                   segs.push_back(seg);    }    }    }
               for ( size_t w = 0; w < segs.size( ); w++ )
               {    SEGS_b[bb].push( True, rid, segs[w].first, segs[w].second, 
                         segs[w].third );    }
               i = j - 1;    }    }

     // Merge results of batches.

     for ( int i = 0; i < nbatches; i++ )
     {    SEGS.append( SEGS_b[i] );
          Destroy( SEGS_b[i] );    }

     // Harvest coverage from the segments.

     if (get_calls)
     {    for ( size_t w = 0; w < SEGS.size( ); w++ )
          {    int rpos = SEGS[w].rpos, u = SEGS[w].u; 
               int upos = SEGS[w].upos; 
               int64_t rid = SEGS[w].rid;
               for ( int z = rpos; z < reads[rid].isize( ); z++ )
               {    if ( upos >= 0 )
                    {    int b = reads[rid][z];
                         calls[b][u][upos]++;
                         qcalls[b][u][upos] += quals[rid][z];
                         qtotal += 2 * quals[rid][z];
                         if ( unibases[u][upos] == b ) qgood += 2 * quals[rid][z];
                         total += 2;    }
                    upos++;
                    if ( upos == unibases[u].isize( ) ) break;    }    }    }    }

// Compute unibase copy numbers.  Derived with as few changes as possible from
// UnibaseCopyNumber2.

void UnibaseCopyNumbersFromSegAligns( const int K, const vecbasevector& unibases, 
     const vec<int>& to_rc, const int PLOIDY, const String run_dir, 
     const String& FRAG_READS, const vec<segalign>& SEGS, 
     const vec<size_t>& S_START, vec<int>& CN )
{
     size_t nuni = unibases.size( );
     CN.resize(nuni);
     vec<int> readlen;
     {    cout << Date( ) << ": loading reads and computing read lengths" << endl;
          vecbasevector reads( run_dir + "/" + FRAG_READS + ".fastb" );
          readlen.resize( reads.size( ) );
          for ( size_t i = 0; i < reads.size( ); i++ )
               readlen[i] = reads[i].size( );    }
     Bool CN_VERBOSE = False;
     Ofstream( cn_out, "/dev/null" );
     int min_length, nlongest;
     GetMinLength( unibases, min_length, nlongest );
     double occCnPloidy = -1.0;
     vec<int64_t> n_kmer_hits, biasOccs( K+1, 1 ), biasInst( K+1, 0 );
     vec<double> biasCurveLoc( K+1, 1.0 );
     for ( int pass = 1; pass <= 2; pass++ )
     {    cout << Date( ) << ": start copy number loop 1, pass = " << pass << endl;
          n_kmer_hits.resize_and_set( nuni, 0 );
          #pragma omp parallel for
          for ( size_t u = 0; u < nuni; u++ )
          {    ForceAssertLt( u, unibases.size( ) );
               if ( pass == 1 && unibases[u].isize( ) < min_length ) continue;
               vec<int> count( unibases[u].isize( ) - K + 1, 0 );
               for ( size_t l = S_START[u]; l < S_START[u+1]; l++ )
               {    int rpos = SEGS[l].rpos, upos = SEGS[l].upos; 
                    size_t rid = SEGS[l].rid;
                    ForceAssertLt( rid, readlen.size( ) );
                    for ( int z = rpos; z <= readlen[rid] - K; z++ )
                    {    if ( upos > unibases[u].isize( ) - K ) break;
                         if ( upos >= 0 ) count[upos]++;
                         upos++;    }    }
               int gc = 0;
               for ( int j = 0; j < K; j++ )
                    if ( IsGC( unibases[u][j] ) ) gc++;
               #pragma omp critical
               {    for ( int z = 0; z <= unibases[u].isize( ) - K; z++ )
                    {    if ( z > 0 ) 
                         {    if ( IsGC( unibases[u][z-1] ) ) gc--;
                              if ( IsGC( unibases[u][z+K-1] ) ) gc++;    }
                         if ( pass == 1 )
                         {    n_kmer_hits[u] += count[z];
                              n_kmer_hits[ to_rc[u] ] += count[z];
                              biasOccs[gc] += 2 * count[z];
                              biasInst[gc]++;    }
                         else
                         {    double freq = double(count[z]) / biasCurveLoc[gc];
                              n_kmer_hits[u] += freq;
                              n_kmer_hits[ to_rc[u] ] += freq;    }    }    }    }
          if ( pass == 1 )
          {    cout << Date( ) << ": computing bias" << endl;
               vec<double> avg_kmer_freq;
               for ( size_t u = 0; u < nuni; u++ ) 
               {    double n_kmers = unibases[u].size() - K + 1;
                    avg_kmer_freq.push_back( n_kmer_hits[u] / n_kmers );    }
               occCnPloidy = double(Sum(avg_kmer_freq)) / (double) nlongest; 
               occCnPloidy /= (double) PLOIDY;
               ComputeBias( K, occCnPloidy, biasOccs, biasInst, biasCurveLoc, 
                    CN_VERBOSE, cn_out );    }    }
     VecPdfEntryVec CN_pdf;
     vec<longlong> cnHisto;
     cout << Date( ) << ": calling ComputeProbs" << endl;
     ComputeProbs( K, PLOIDY, occCnPloidy, THRESH_DEFAULT, ERR_RATE_DEFAULT, 
          unibases, n_kmer_hits, CN_pdf, cnHisto, CN_VERBOSE, cn_out );
     cout << Date( ) << ": calling GetMostLikelyValue" << endl;
     #pragma omp parallel for
     for ( size_t i = 0; i < nuni; i++ )
          GetMostLikelyValue( CN[i], CN_pdf[i] );    }

Bool cmp_rid( const segalign& a1, const segalign& a2 )
{    return a1.rid < a2.rid;    }

void MakeUlinks( const String& run_dir, const String& TMP, const String& FRAG_READS,
     const String& JUMP_READS, const int K, const int PLOIDY, const vec<int>& CN, 
     const vec<int>& to_rc, vec<segalign>& SEGS, vec<segalign>& JSEGS, 
     const vecbasevector& unibases, const vec<int>& innie_sep, 
     const vec<int>& innie_dev, const vec<double>& innie_percent,
     const double min_innie_percent, const Bool log_links_detailed, 
     const int min_kmers, const int max_devs, vec<ulink_with_uids>& ulinks, 
     const Bool CHECKPOINT )
{
     String ULINKS_file = run_dir + "/" + TMP + "/UnipathPatcher.ULINKS";
     vec<String> unknown_libs;
     if ( CHECKPOINT && IsRegularFile(ULINKS_file) && IsRegularFile(ULINKS_file) )
          BinaryRead3( ULINKS_file, ulinks );
     else
     {    // Only running for jump pairs for now.  For some reason using fragment
          // pairs results in a lot of noise.  Also note that now handling of 
          // "innies" would be nonsensical for fragment pairs - would need to fix.

          // for ( int pass = 1; pass <= 2; pass++ )
          for ( int pass = 2; pass <= 2; pass++ )
          {    cout << Date( ) << ": pass = " << pass << ", set up for link "
                    << "computation" << endl;
               String head = pass == 1 
                    ? run_dir + "/" + FRAG_READS : run_dir + "/" + JUMP_READS;
               vecbasevector reads( head + ".fastb" );
               uint64_t nreads = reads.size( );
               vec<int> readlen( reads.size( ) );
               for ( size_t i = 0; i < reads.size( ); i++ )
                    readlen[i] = reads[i].size( );
               if ( !log_links_detailed ) Destroy(reads);
               PairsManager pairs( head + ".pairs" );
               vec<segalign>& SUGS = ( pass == 1 ? SEGS : JSEGS );
               ParallelSort( SUGS, cmp_rid );
               vec<size_t> S_START_RID(nreads+1);
               {    size_t SEG_POS = 0;
                    for ( uint64_t rid = 0; rid <= nreads; rid++ )
                    {    while( SEG_POS < SUGS.size( ) && SUGS[SEG_POS].rid < rid ) 
                              ++SEG_POS;
                         S_START_RID[rid] = SEG_POS;    }    }
               cout << Date( ) << ": pass = " << pass << ", compute links from "
                    << ToStringAddCommas( pairs.nPairs( ) ) << " pairs" << endl;
               for ( size_t i = 0; i < pairs.nPairs( ); i++ )
               {    longlong id1 = pairs.ID1(i), id2 = pairs.ID2(i);
     
                    // Require that jump reads have only one placement.  This may be
                    // too stringent.

                    if ( pass == 2 ) 
                    {    if ( S_START_RID[id1+1] - S_START_RID[id1] > 1 ) continue;
                         if ( S_START_RID[id2+1] - S_START_RID[id2] > 1 ) 
                              continue;    }

                    for (size_t j1 = S_START_RID[id1]; j1 < S_START_RID[id1+1]; j1++)
                    {    for (size_t j2 = S_START_RID[id2]; 
                              j2 < S_START_RID[id2+1]; j2++)
                         {    int u1 = SUGS[j1].u, u2 = to_rc[ SUGS[j2].u ];
                              if ( u1 == u2 || u1 == to_rc[u2] ) continue;
                              if ( CN[u1] > PLOIDY ) continue;
                              if ( CN[u2] > PLOIDY ) continue;
                              if ( unibases[u1].isize( ) - K + 1 < min_kmers ) 
                                   continue;
                              if ( unibases[u2].isize( ) - K + 1 < min_kmers ) 
                                   continue;
                              int nu1 = unibases[u1].size( ); 
                              int nu2 = unibases[u2].size( );
                              int nr1 = readlen[id1], nr2 = readlen[id2];
                              int upos1 = SUGS[j1].upos, upos2 = SUGS[j2].upos;
                              if (log_links_detailed)
                              {    PRINT4( upos1, SUGS[j1].rpos, upos2, 
                                        SUGS[j2].rpos );    }
                              int lib = pairs.libraryID(i);
  		              // Library -> separation/deviation
                              int sep_in = innie_sep[lib] + nr1 + nr2
                                   - ( nu1 - upos1 + SUGS[j1].rpos + nu2 
                                   - upos2 + SUGS[j2].rpos );
                              int dev_in = innie_dev[lib];
                              if ( pass == 1 )
                              {    if ( sep_in + max_devs * dev_in >= -K+1 )
                                   {    ulinks.push( u1, u2, sep_in, dev_in,
                                             Min( nu1, upos1 + nr1 ),
                                             Max( 0, nu2 - upos2 - nr2 ), 1 );
                                        ulinks.push( to_rc[u2], to_rc[u1], sep_in, 
                                             dev_in, Min( nu2, upos2 + nr2 ), 
                                             Max( 0, nu1 - upos1 - nr1 ), 1 );
                                        if (log_links_detailed)
                                        {    PRINT4( u1, u2, sep_in, dev_in );
                                             PRINT4( to_rc[u2], to_rc[u1], sep_in, 
                                                  dev_in );
                                             reads[id1].Print( cout, id1 );
                                             reads[id2].Print( cout, id2 );
                                             cout << "\n";    }    }    }
                              else
                              {    int sep_out
                                        = pairs.getLibrarySep(pairs.libraryID(i))
                                        - ( upos1 - SUGS[j1].rpos
                                        + upos2 - SUGS[j2].rpos ) + nr1 + nr2;
                                   int dev_out
                                        = pairs.getLibrarySD(pairs.libraryID(i));
                                   if (log_links_detailed)
                                   {    PRINT4( u1, u2, sep_in, dev_in );
                                        PRINT4( u2, u1, sep_out, dev_out );
                                        PRINT4( to_rc[u2], to_rc[u1], 
                                             sep_in, dev_in );
                                        PRINT4( to_rc[u1], to_rc[u2], sep_out, 
                                             dev_out );
                                        reads[id1].Print( cout, id1 );
                                        reads[id2].Print( cout, id2 );
                                        cout << "\n";    }
			           
			           // Add this link to the ulinks list - in both fw
			           // and rc form, as both an innie and an outie.
     
                                   if ( innie_percent[lib] >= min_innie_percent )
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
                                        Max(0, upos2), 1 );    }    }    }    }    }
          cout << Date( ) << ": sorting " << ToStringAddCommas( ulinks.size( ) ) 
               << " unipath links" << endl;
          ParallelSort(ulinks);
          if (CHECKPOINT) BinaryWrite3( ULINKS_file, ulinks );    }    }

#define MEM ToStringAddCommas( MemUsageBytes( ) )

void BuildJoinData( const vec<int>& dead_fw, const vec<int>& dead_rc,
     const vec< pair<int,int> >& to_process_sepdev, const vecbasevector& unibases, 
     vec<int>& to_rc, const String& run_dir, const String& PATCHDIR,
     const String& FRAG_READS, const String& JUMP_READS, vec<segalign>& SEGS, 
     vec<segalign>& JSEGS, vec<size_t>& S_START,
     vec<size_t>& JS_START, const vec<int>& innie_sep,
     const vec<int>& innie_dev, const vec<double>& innie_percent, 
     const double min_innie_percent, vec<size_t>& join_data_offsets,
     const vec< pair<int,int> >& joiners, const String& JOINDATA_file )
{    
     cout << Date( ) << ": finding reads near dead ends, memory usage = "
          << MEM << endl;
     const int prox = 30;
     int nde = dead_fw.size( ) + dead_rc.size( );
     String patch = run_dir + "/" + PATCHDIR;
     String extenders_b_frag_file = patch + "/extenders_b_frag_file";
     String extenders_q_frag_file = patch + "/extenders_q_frag_file";
     String extenders_b_jump_file = patch + "/extenders_b_jump_file";
     String extenders_q_jump_file = patch + "/extenders_q_jump_file";
     vec< vec<opair> > extenders_frag(nde), extenders_jump(nde);
     int nd = dead_fw.size( );
     {    vecbasevector reads( run_dir + "/" + FRAG_READS + ".fastb" );
          {    PairsManager pairs( run_dir + "/" + FRAG_READS + ".pairs" );
               pairs.makeCache( );
               vec<int> INNIE_SEP, INNIE_DEV;
               for ( size_t i = 0; i < pairs.nLibraries( ); i++ )
               {    cout << Date( ) << ": see fragment library "
                         << pairs.getLibraryName(i) << ", "
                         << pairs.getLibrarySep(i) << " +/- "
                         << pairs.getLibrarySD(i) << endl;
                    INNIE_SEP.push_back( pairs.getLibrarySep(i) );
                    INNIE_DEV.push_back( pairs.getLibrarySD(i) );    }
               #pragma omp parallel for
               for ( int i = 0; i < dead_fw.isize( ); i++ )
               {    int u = dead_fw[i];
                    for ( size_t j = S_START[u]; j < S_START[u+1]; j++ )
                    {    const segalign& a = SEGS[j];
                         int stop = a.upos + reads[a.rid].size( );
                         if ( stop >= unibases[u].isize( ) - prox )
                         {    int p = pairs.getPairID(a.rid);
                              if ( p < 0 ) continue;
                              int lib = pairs.libraryID(p);
                              int64_t id1 = a.rid, id2 = pairs.getPartnerID(a.rid);
                              int sep = INNIE_SEP[lib], dev = INNIE_DEV[lib];
                              if ( !a.fw ) swap( id1, id2 );
                              extenders_frag[i].push( 
                                   id1, id2, sep, dev );    }    }
                         u = to_rc[ dead_fw[i] ];
                    for ( size_t j = S_START[u]; j < S_START[u+1]; j++ )
                    {    const segalign& a = SEGS[j];
                         int start = a.upos;
                         if ( start <= prox )
                         {    int p = pairs.getPairID(a.rid);
                              if ( p < 0 ) continue;
                              int64_t id1 = a.rid, id2 = pairs.getPartnerID(a.rid);
                              int lib = pairs.libraryID(p);
                              int sep = INNIE_SEP[lib], dev = INNIE_DEV[lib];
                              if ( a.fw ) swap( id1, id2 );
                                   extenders_frag[i].push( 
                                   id1, id2, sep, dev );    }    }    }
               #pragma omp parallel for
               for ( int i = 0; i < dead_rc.isize( ); i++ )
               {    int u = dead_rc[i];
                    for ( size_t j = S_START[u]; j < S_START[u+1]; j++ )
                    {    const segalign& a = SEGS[j];
                         int start = a.upos;
                         if ( start <= prox )
                         {    int p = pairs.getPairID(a.rid);
                              if ( p < 0 ) continue;
                              int64_t id1 = a.rid, id2 = pairs.getPartnerID(a.rid);
                              int lib = pairs.libraryID(p);
                              int sep = INNIE_SEP[lib], dev = INNIE_DEV[lib];
                              if ( !a.fw ) swap( id1, id2 );
                              extenders_frag[nd+i].push( 
                                   id1, id2, sep, dev );    }    }
                    u = to_rc[ dead_rc[i] ];
                    for ( size_t j = S_START[u]; j < S_START[u+1]; j++ )
                    {    const segalign& a = SEGS[j];
                         int stop = a.upos + reads[a.rid].size( );
                         if ( stop >= unibases[u].isize( ) - prox )
                         {    int p = pairs.getPairID(a.rid);
                              if ( p < 0 ) continue;
                              int64_t id1 = a.rid, id2 = pairs.getPartnerID(a.rid);
                              int lib = pairs.libraryID(p);
                              int sep = INNIE_SEP[lib], dev = INNIE_DEV[lib];
                              if ( a.fw ) swap( id1, id2 );
                              extenders_frag[nd+i].push( 
                                   id1, id2, sep, dev );    }    }    }
               #pragma omp parallel for
               for ( int i = 0; i < extenders_frag.isize( ); i++ )
                    UniqueSort( extenders_frag[i] );    }
          Destroy(SEGS), Destroy(S_START);
          cout << Date( ) << ": building extenders_b_frag, memory usage = "
               << MEM << endl;
          {    bvec3 extenders_b_frag(nde);
               // DO NOT parallelize this with omp! it'll be slower 
               for ( int i = 0; i < extenders_frag.isize( ); i++ )
               {    for ( int j = 0; j < extenders_frag[i].isize( ); j++ )
                    {    const opair& p = extenders_frag[i][j];
                         basevector b = reads[p.id2];
                         b.ReverseComplement( );
                         extenders_b_frag[i].push_back( reads[p.id1] );
                         extenders_b_frag[i].push_back( b );    }    }
               BinaryWriter::writeFile(extenders_b_frag_file.c_str(),
                                       extenders_b_frag);
               cout << Date( ) << ": now memory usage = " << MEM << endl;    }
          Destroy(reads);
          vecqualvector quals;
          vec<int64_t> reads_needed;
          for ( int i = 0; i < extenders_frag.isize( ); i++ )
          {    for ( int j = 0; j < extenders_frag[i].isize( ); j++ )
               {    const opair& p = extenders_frag[i][j];
                    reads_needed.push_back( p.id1, p.id2 );    }    }
          UniqueSort(reads_needed);
          quals.Read( run_dir + "/" + FRAG_READS + ".qualb", reads_needed );
          cout << Date( ) << ": building extenders_q_frag, memory usage = "
               << MEM << endl;
          {    qvec3 extenders_q_frag(nde);
               // DO NOT parallelize this with omp! it'll be slower 
               for ( int i = 0; i < extenders_frag.isize( ); i++ )
               {    for ( int j = 0; j < extenders_frag[i].isize( ); j++ )
                    {    const opair& p = extenders_frag[i][j];
                         size_t xid1 = BinPosition( reads_needed, p.id1 );
                         size_t xid2 = BinPosition( reads_needed, p.id2 );
                         extenders_q_frag[i].push_back( quals[xid1] );
                         extenders_q_frag[i].push_back( 
                              Reverse( quals[xid2] ) );    }    }    
               BinaryWriter::writeFile(extenders_q_frag_file.c_str(),
                                       extenders_q_frag);
               cout << Date( ) << ": now memory usage = " << MEM << endl;    }    }

     // Gather up jumps.

     {    vecbasevector jreads( run_dir + "/" + JUMP_READS + ".fastb" );
          for ( size_t i = 0; i < jreads.size( ); i++ )
               jreads[i].ReverseComplement( );
          {    PairsManager jpairs( run_dir + "/" + JUMP_READS + ".pairs" );
               jpairs.makeCache( );
               for ( size_t i = 0; i < jpairs.nLibraries( ); i++ )
               {    cout << Date( ) << ": see jump library "
                         << jpairs.getLibraryName(i) << ", "
                         << jpairs.getLibrarySep(i) << " +/- "
                         << jpairs.getLibrarySD(i) << endl;    }
               #pragma omp parallel for
               for ( int i = 0; i < dead_fw.isize( ); i++ )
               {    int u = dead_fw[i];
                    for ( size_t j = JS_START[u]; j < JS_START[u+1]; j++ )
                    {    const segalign& a = JSEGS[j];
                         int stop = a.upos + jreads[a.rid].size( );
                         if ( stop >= unibases[u].isize( ) - prox )
                         {    int p = jpairs.getPairID(a.rid);
                              if ( p < 0 ) continue;
                              int lib = jpairs.libraryID(p);
                              int64_t id1 = a.rid, id2 = jpairs.getPartnerID(a.rid);
                              int sep = innie_sep[lib], dev = innie_dev[lib];
                              if ( innie_percent[lib] < min_innie_percent )
                              {    sep = 1000000000, dev = 1000;    } // nonsense
                              if ( !a.fw ) swap( id1, id2 );
                              extenders_jump[i].push( 
                                   id1, id2, sep, dev );    }    }
                    u = to_rc[ dead_fw[i] ];
                    for ( size_t j = JS_START[u]; j < JS_START[u+1]; j++ )
                    {    const segalign& a = JSEGS[j];
                         int start = a.upos;
                         if ( start <= prox )
                         {    int p = jpairs.getPairID(a.rid);
                              if ( p < 0 ) continue;
                              int lib = jpairs.libraryID(p);
                              int64_t id1 = a.rid, id2 = jpairs.getPartnerID(a.rid);
                              int sep = innie_sep[lib], dev = innie_dev[lib];
                              if ( innie_percent[lib] < min_innie_percent )
                              {    sep = 1000000000, dev = 1000;    } // nonsense
                              if ( a.fw ) swap( id1, id2 );
                              extenders_jump[i].push( 
                                   id1, id2, sep, dev );    }    }    }
               #pragma omp parallel for
               for ( int i = 0; i < dead_rc.isize( ); i++ )
               {    int u = dead_rc[i];
                    for ( size_t j = JS_START[u]; j < JS_START[u+1]; j++ )
                    {    const segalign& a = JSEGS[j];
                         int start = a.upos;
                         if ( start <= prox )
                         {    int p = jpairs.getPairID(a.rid);
                              if ( p < 0 ) continue;
                              int lib = jpairs.libraryID(p);
                              int64_t id1 = a.rid, id2 = jpairs.getPartnerID(a.rid);
                              int sep = innie_sep[lib], dev = innie_dev[lib];
                              if ( innie_percent[lib] < min_innie_percent )
                              {    sep = 1000000000, dev = 1000;    } // nonsense
                              if ( !a.fw ) swap( id1, id2 );
                              extenders_jump[nd+i].push( 
                                   id1, id2, sep, dev );    }    }
                    u = to_rc[ dead_rc[i] ];
                    for ( size_t j = JS_START[u]; j < JS_START[u+1]; j++ )
                    {    const segalign a = JSEGS[j];
                         int stop = a.upos + jreads[a.rid].size( );
                         if ( stop >= unibases[u].isize( ) - prox )
                         {    int p = jpairs.getPairID(a.rid);
                              if ( p < 0 ) continue;
                              int lib = jpairs.libraryID(p);
                              int64_t id1 = a.rid, id2 = jpairs.getPartnerID(a.rid);
                              int sep = innie_sep[lib], dev = innie_dev[lib];
                              if ( innie_percent[lib] < min_innie_percent )
                              {    sep = 1000000000, dev = 1000;    } // nonsense
                              if ( a.fw ) swap( id1, id2 );
                              extenders_jump[nd+i].push( 
                                   id1, id2, sep, dev );    }    }    }
               #pragma omp parallel for
               for ( int i = 0; i < extenders_jump.isize( ); i++ )
                    UniqueSort( extenders_jump[i] );    }
          Destroy(JSEGS), Destroy(JS_START);
          cout << Date( ) << ": building extenders_b_jump, memory usage = "
               << MEM << endl;
          {    bvec3 extenders_b_jump(nde);
               // DO NOT parallelize this with omp! it'll be slower 
               for ( int i = 0; i < extenders_jump.isize( ); i++ )
               {    for ( int j = 0; j < extenders_jump[i].isize( ); j++ )
                    {    const opair& p = extenders_jump[i][j];
                         basevector b = jreads[p.id2];
                         b.ReverseComplement( );
                         extenders_b_jump[i].push_back( jreads[p.id1] );
                         extenders_b_jump[i].push_back( b );    }    }
               BinaryWriter::writeFile(extenders_b_jump_file.c_str(),
                                       extenders_b_jump);
               cout << Date( ) << ": now memory usage = " << MEM << endl;    }
          Destroy(jreads);
          cout << Date( ) << ": loading jquals" << endl;
          vecqualvector jquals;
          vec<int64_t> reads_needed;
          for ( int i = 0; i < extenders_jump.isize( ); i++ )
          {    for ( int j = 0; j < extenders_jump[i].isize( ); j++ )
               {    const opair& p = extenders_jump[i][j];
                    reads_needed.push_back( p.id1, p.id2 );    }    }
          UniqueSort(reads_needed);
          jquals.Read( run_dir + "/" + JUMP_READS + ".qualb", reads_needed );
          for ( size_t i = 0; i < jquals.size( ); i++ )
               jquals[i].ReverseMe( );
          cout << Date( ) << ": building extenders_q_jump, memory usage = "
               << MEM << endl;
          {    qvec3 extenders_q_jump(nde);
               // DO NOT parallelize this with omp! it'll be slower 
               for ( int i = 0; i < extenders_jump.isize( ); i++ )
               {    for ( int j = 0; j < extenders_jump[i].isize( ); j++ )
                    {    const opair& p = extenders_jump[i][j];
                         size_t xid1 = BinPosition( reads_needed, p.id1 );
                         size_t xid2 = BinPosition( reads_needed, p.id2 );
                         extenders_q_jump[i].push_back( jquals[xid1] );
                         extenders_q_jump[i].push_back( 
                              Reverse( jquals[xid2] ) );    }    }
               BinaryWriter::writeFile(extenders_q_jump_file.c_str(),
                                       extenders_q_jump);
               cout << Date( ) << ": now memory usage = " << MEM << endl;    }    }
     cout << Date( ) << ": found extenders, memory usage = " << MEM << endl;

     // Build input data for joins.

     bvec3 extenders_b_frag, extenders_b_jump;
     qvec3 extenders_q_frag, extenders_q_jump;
     BinaryReader::readFile( extenders_b_frag_file.c_str(),&extenders_b_frag );
     BinaryReader::readFile( extenders_b_jump_file.c_str(),&extenders_b_jump );
     BinaryReader::readFile( extenders_q_frag_file.c_str(),&extenders_q_frag );
     BinaryReader::readFile( extenders_q_jump_file.c_str(),&extenders_q_jump );
     cout << Date( ) << ": reloaded extending data, memory usage = " << MEM << endl;

     BinaryWriter jdWriter(JOINDATA_file.c_str(),false);
     size_t offset = jdWriter.write(PCottageJoinData::HEADER);
     PCottageJoinData joinData;
     size_t nnn = joiners.size();
     join_data_offsets.clear();
     join_data_offsets.reserve(nnn);

     for ( size_t i = 0; i < nnn; i++ )
     {
         size_t v1 = joiners[i].first;
         size_t v2 = joiners[i].second;

         joinData.sep = to_process_sepdev[v1].first;
         joinData.dev = to_process_sepdev[v1].second;
         joinData.L = unibases[dead_fw[v1]];
         joinData.R = unibases[dead_rc[v2-dead_fw.size()]];

          // Combine extenders for v1 and v2.

         vec<opair> const& extF1 = extenders_frag[v1];
         vec<opair> const& extF2 = extenders_frag[v2];
         vec<opair> const& extJ1 = extenders_jump[v1];
         vec<opair> const& extJ2 = extenders_jump[v2];
         size_t totSize = extF1.size()+extF2.size()+extJ1.size()+extJ2.size();
         joinData.reads.clear().reserve(2*totSize);
         joinData.quals.clear().reserve(2*totSize);
         joinData.pairs.clear();
         joinData.pairs.reserve(totSize);

         joinData.reads = extenders_b_frag[v1];
         joinData.quals = extenders_q_frag[v1];
         for ( size_t j = 0; j < extF1.size(); j++ )
         {    const opair& p = extF1[j];
              joinData.pairs.push( p.sep, p.dev );    }
         AssertEq(joinData.reads.size(),joinData.quals.size());
         AssertEq(joinData.reads.size(),2*joinData.pairs.size());


         bvec3::const_reference extBF2 = extenders_b_frag[v2];
         qvec3::const_reference extQF2 = extenders_q_frag[v2];
         for ( size_t j = 0; j < extF2.size( ); j++ )
         {    const opair& p = extF2[j];
              if ( BinMember(extF1,p) ) continue;
              size_t idx = 2*j;
              joinData.reads.push_back( extBF2[idx] );
              joinData.reads.push_back( extBF2[idx+1] );
              joinData.quals.push_back( extQF2[idx] );
              joinData.quals.push_back( extQF2[idx+1] );
              joinData.pairs.push( p.sep, p.dev );    }
         AssertEq(joinData.reads.size(),joinData.quals.size());
         AssertEq(joinData.reads.size(),2*joinData.pairs.size());

         bvec3::const_reference bj = extenders_b_jump[v1];
         joinData.reads.append( bj.begin(), bj.end() );
         qvec3::const_reference qj = extenders_q_jump[v1];
         joinData.quals.append( qj.begin(), qj.end() );
         for ( size_t j = 0; j < extJ1.size( ); j++ )
         {    const opair& p = extJ1[j];
              joinData.pairs.push( p.sep, p.dev );    }
         AssertEq(joinData.reads.size(),joinData.quals.size());
         AssertEq(joinData.reads.size(),2*joinData.pairs.size());

         bvec3::const_reference extBJ2 = extenders_b_jump[v2];
         qvec3::const_reference extQJ2 = extenders_q_jump[v2];
         for ( size_t j = 0; j < extJ2.size( ); j++ )
         {    const opair& p = extJ2[j];
              if ( BinMember(extJ1,p) ) continue;
              size_t idx = 2*j;
              joinData.reads.push_back( extBJ2[idx] );
              joinData.reads.push_back( extBJ2[idx+1] );
              joinData.quals.push_back( extQJ2[idx] );
              joinData.quals.push_back( extQJ2[idx+1] );
              joinData.pairs.push( p.sep, p.dev );    }
         AssertEq(joinData.reads.size(),joinData.quals.size());
         AssertEq(joinData.reads.size(),2*joinData.pairs.size());

         join_data_offsets.push_back(offset);
         offset += jdWriter.write(joinData); }
     jdWriter.close();
     Destroy(to_rc);    }

namespace
{
    class ResultsProc
    {
    public:
        ResultsProc( vec< vec<int> >& startStop, vec<String>& reports,
                     int& nPatches, vec<Bool>& joined, vec<Bool>& perfect,
                     bvec3& newStuff )
        : mStartStop(startStop), mReports(reports), mNPatches(nPatches),
          mJoined(joined), mPerfect(perfect), mNewStuff(newStuff) {}

        // compiler-supplied defaults for copying and the destructor are OK

        void operator()( PCottageResults const& results )
        {
            if ( results.report == "FAILED" )
                FatalErr("Inconsistent data discovered by item "
                            << results.itemNumber);

            mStartStop[results.itemNumber] = results.startStop;
            mReports[results.itemNumber] = results.report;

            ForceAssertEq(results.reads.size()%3,0ul);
            size_t nPatches = results.reads.size()/3;
            if ( nPatches )
            {
                mNPatches += nPatches;
                mJoined[results.itemNumber] = True;
                if ( results.report.Contains("JOIN 1") &&
                     results.report.RevAfter("JOIN 1").Contains("perfect") )
                    mPerfect[results.itemNumber] = True;
                bvec3::reference newStuff = mNewStuff[results.itemNumber];
                newStuff.reserve(newStuff.size()+6*nPatches);
                vecbvec::const_iterator iItr(results.reads.begin());
                for ( size_t iii = 0; iii < nPatches; ++iii )
                {
                    newStuff.append(iItr,iItr+3);
                    newStuff.append(iItr,iItr+3);
                    vecbvec::iterator end(newStuff.end());
                    for ( vecbvec::iterator itr(end-3); itr != end; ++itr )
                        itr->ReverseComplement();
                    iItr += 3;
                }
            }
        }

    private:
        vec< vec<int> >& mStartStop;
        vec<String>& mReports;
        int& mNPatches;
        vec<Bool>& mJoined;
        vec<Bool>& mPerfect;
        bvec3& mNewStuff;
    };
}

void MakeJoins( const vec< pair<int,int> >& Xjoiners, const vecbasevector& X,
     const String& progname, const int K, const String& data_dir, 
     const String& run_dir, const String& PATCHDIR, int& attempted_joins,
     vec<Bool>& joined, int& npatches, vec<Bool>& perfect,
     bvec3& new_stuff, vec< vec<int> >& start_stop,
     vec<String>& reports, unsigned num_procs, const vec<int>& joins,
     const String& LOG, const Bool log_attempted_joins, const Bool log_action, 
     const int MAX_PATHS, const int MAX_PATH_ITERATIONS, const int MAX_EDGES, 
     const int MAX_READS, const int MAX_JOINS, const int MIN_OVERLAP_END, 
     const String& JOINDATA, const vec<size_t>& join_data_offsets,
     const vec< pair<int,int> >& LR_to_process )
{
    std::cout << Date() << ": considering " << Xjoiners.size() << " joins"
              << std::endl;
    new_stuff.clear();
    new_stuff.resize(Xjoiners.size());
    start_stop.clear();
    start_stop.resize(Xjoiners.size());
    reports.clear();
    reports.resize(Xjoiners.size());
    attempted_joins = 0;
    perfect.resize_and_set(Xjoiners.size(), False);
    joined.resize_and_set(Xjoiners.size(), False);

    size_t nnn = Xjoiners.size();
    std::vector<PCottageWhichJoin> toDo;
    toDo.reserve(Xjoiners.size());
    for ( size_t i = 0; i < nnn; ++i )
    {
        if ( (joins.empty() || BinMember(joins, i)) && (LR_to_process.empty()
                || Member(LR_to_process, Xjoiners[i])) )
        {
            int u1 = Xjoiners[i].first;
            int u2 = Xjoiners[i].second;
            PCottageWhichJoin whichJoin;
            whichJoin.itemNumber = i;
            whichJoin.offset = join_data_offsets[i];
            whichJoin.u1 = u1;
            whichJoin.u2 = u2;
            toDo.push_back(whichJoin);
            if ( log_attempted_joins )
            {
                reports[i] += "\ntrying to join " + ToString(u1) + "["
                        + ToString(X[u1].size()) + " bases]" + " to "
                        + ToString(u2) + "[" + ToString(X[u2].size())
                        + " bases]" + "\n";
            }
        }
    }

    if ( toDo.empty() )
    {
        std::cout << "There are no joins to attempt." << std::endl;
        return;
    }

    attempted_joins = toDo.size();
    if ( toDo.size() < num_procs )
        num_procs = toDo.size();

    // Launch workers.
    std::cout << Date() << ": forking " << num_procs << " cottage workers to do "
              << toDo.size() << " joins." << std::endl;

    if ( true )
    {
        String ROOT = run_dir + "/" + PATCHDIR + "/" + progname + ".";
        vecString args;
        args.reserve(10);
        args.push_back(String("PatcherCottage"));
        args.push_back(String("K=")+ToString(K));
        args.push_back(String("JOINDATA=")+JOINDATA);
        args.push_back(String("ROOT=")+ROOT);
        args.push_back(String("data_dir=")+data_dir);
        args.push_back(String("MAX_READS=")+ToString(MAX_READS));
        args.push_back(String("MAX_JOINS=")+ToString(MAX_JOINS));
        args.push_back(String("MIN_OVERLAP_END=")+ToString(MIN_OVERLAP_END));
        args.push_back(String("LOG=")+LOG);
        args.push_back(String("NH=True"));
        typedef WorklistMP<PCottageWhichJoin, PCottageResults> WL;
        WL worklist(num_procs, args);
        npatches = 0;
        ResultsProc resultsProc(start_stop, reports, npatches, joined, perfect,
                new_stuff);

        std::cout << Date() << ": farming out work" << std::endl;
        if ( !worklist.doWork(toDo.begin(), toDo.end(), resultsProc) )
        {    cout << "\nSorry, the workers have rebelled.  Something has gone "
                  << "badly wrong.  Aborting....\n";
             TracebackThisProcess( );    }
    }

    std::cout << Date() << ": done" << std::endl;
}

void BuildAlignments( const vecbasevector& reads, const int seed, 
     const int min_align_length, const double max_align_errors, 
     vec<xalign>& xaligns )
{    xaligns.clear( );
     int nreads = reads.size( );
     static vec<int32_t> mmers, pos, ids;
     mmers.clear( ), pos.clear( ), ids.clear( );
     for ( int id = 0; id < nreads; id++ )
     {    for ( int j = 0; j <= reads[id].isize( ) - seed; j++ )
          {    int32_t x = reads[id].extractKmer( j, seed );
               mmers.push_back(x), ids.push_back(id), pos.push_back(j);    }    }
     SortSync( mmers, ids, pos );
     vec< set< triple<int32_t,int32_t,int32_t> > > yaligns(nreads);
     for ( int i = 0; i < mmers.isize( ); i++ )
     {    int j = mmers.NextDiff(i);
          for ( int l1 = i; l1 < j; l1++ )
          {    for ( int l2 = i+1; l2 < j; l2++ )
               {    int32_t id1 = ids[l1], id2 = ids[l2];
                    if ( id1 == id2 ) continue;
                    int pos1 = pos[l1], pos2 = pos[l2];
                    if ( pos2 > pos1 ) { swap(id1,id2); swap(pos1,pos2); }
                    int len1 = reads[id1].size(), len2 = reads[id2].size();
                    int start1 = pos1 - pos2, start2 = 0;
                    int stop1 = Min( len1, len2 + start1 );
                    int align_length = stop1 - start1;
                    if ( align_length < min_align_length ) continue;
                    triple<int32_t,int32_t,int32_t> t1( id2, start1, start2 );
                    if ( Member( yaligns[id1], t1 ) ) continue;
                    int stop2 = Min( len2, len1 - start1 ), errs = 0;
                    for ( int p1 = start1; p1 < stop1; p1++ )
                    {    int p2 = p1 - start1;
                         if ( reads[id1][p1] != reads[id2][p2] ) errs++;    } 
                    if ( double(errs)/double(align_length) > max_align_errors )
                         continue;
                    yaligns[id1].insert(t1);
                    triple<int32_t,int32_t,int32_t> t2( id1, start2, start1 );
                    yaligns[id2].insert(t2);    }    }
          i = j - 1;    }
     for ( int id1 = 0; id1 < nreads; id1++ )
     {    for ( set< triple<int32_t,int32_t,int32_t> >::iterator 
               j = yaligns[id1].begin( ); j != yaligns[id1].end( ); j++ )
          {    xaligns.push( id1, j->first, j->second, j->third );    }    }
     Sort(xaligns);    }

Bool Consistent( const align_chain& c, const xalign& a, const vec<xalign>& xaligns,
     const vec<int>& xaligns_index )
{
     int id1 = c.ids.back( ), id2 = a.id2;
     ForceAssertEq( id1, a.id1 );

     // Test for consistency of a with c.  As part of this, we initially set zpos1
     // and zpos2 to pos1 and pos2 for the alignment of id1 with id2.  Then as we 
     // work our way back through the reads id of c, we update zpos1 and zpos2 so 
     // that they are pos1 and pos2 for the implied alignment of id with id2.
     // Simultaneously, we let left_ext track the leftmost position on id2 that has
     // so far been covered.  As soon as left_ext = 0, we're done confirming.

     int zpos1 = a.pos1, zpos2 = a.pos2;
     int left_ext = a.pos2;
     if ( left_ext > 0 )
     {    for ( int j = c.ids.isize( ) - 2; j >= 0; j-- )
          {    int id = c.ids[j], pos1 = c.pos1[j], pos2 = c.pos2[j];

               // Update zpos1 and zpos2.  

               int zpos2_new_minus_zpos1_new = (zpos2 - zpos1) + (pos2 - pos1);
               if ( zpos2_new_minus_zpos1_new >= 0 )
                    zpos2 = zpos2_new_minus_zpos1_new, zpos1 = 0;
               else zpos1 = -zpos2_new_minus_zpos1_new, zpos2 = 0;

               // Update left_ext.  If read id doesn't go further to the left 
               // on id2 than the previous read, there's nothing to check.

               int left_ext_new = Max( 0, zpos2-zpos1 );
               if ( left_ext_new >= left_ext ) continue;
               left_ext = left_ext_new;

               // Test for existence of the implied alignment between id and
               // id2.  Note that we are guaranteed to meet the minimum 
               // overlap.  Note that this testing is also a bit inefficient
               // because during the entire process, we are not testing a
               // minimal number of times.

               Bool found = False;
               for ( int l = xaligns_index[id]; l < xaligns_index[id+1]; l++ )
               {    const xalign& a = xaligns[l];
                    if ( a.id2 == id2 && a.pos1 == zpos1 && a.pos2 == zpos2 )
                    {    found = True;
                         break;    }    }
               if ( !found ) return False;

               // Test for done.

               if ( left_ext == 0 ) break;    }    }

     // Consistent.

     return True;    }

// GetRights.  Given an align_chain, find all its right extensions, according to
// the following rules:
//
// 1. If x is the last read in the chain, the "overlap" between x and xi must be at 
// least half the overlap between the longest "overlap" x --> xj, amongst reads that
// extend x to the right by at least one base.  Here the "overlap" as defined
// is given by the number of asterisks:
//
//          *************************
//              ---------x---------->
//          -------------xj------------>
//
// 2. Only minimal extensions are returned.
// *** 1 changed, now MIN_OVERLAP_END is enough ***

void GetRights( const align_chain& c, const vec<xalign>& xaligns,
     const vec<int>& xaligns_index, const vec<Bool>& good_right,
     vec<align_chain>& c_next, const int nreads, const int MIN_OVERLAP_END )
{    
     // First find the consistent extensions.  Don't allow reuse of reads.

     int id1 = c.ids.back( );
     int start = xaligns_index[id1], stop = xaligns_index[id1+1];
     vec<Bool> usable( stop - start, False );
     vec<int> cids(c.ids);
     Sort(cids);
     for ( int i = start; i < stop; i++ )
     {    if ( !good_right[i] ) continue;
          if ( BinMember( cids, xaligns[i].id2 ) ) continue;
          Bool consistent = Consistent( c, xaligns[i], xaligns, xaligns_index );
          if (consistent) usable[i-start] = True;    }

     // Now test rule 1.

     int max_overlap = 0;
     for ( int i = 0; i < usable.isize( ); i++ )
     {    if ( !usable[i] ) continue;
          max_overlap = Max( max_overlap, xaligns[start+i].OverlapPlus( ) );    }
     int min_overlap = (max_overlap + 1) / 2;
     vec<int> ids, overlaps;
     for ( int i = 0; i < usable.isize( ); i++ )
     {    if ( !usable[i] ) continue;
          ids.push_back(i);
          overlaps.push_back( xaligns[start+i].OverlapPlus( ) );    }
     ReverseSortSync( overlaps, ids );
     int min_use = 4;
     int used = 0;
     for ( int z = 0; z < ids.isize( ); z++ )
     {    int i = ids[z];
          int min_overlap_this = min_overlap;
          if ( xaligns[start+i].id2 == nreads + 1 )
               min_overlap_this = Min( MIN_OVERLAP_END, min_overlap_this );
          if ( overlaps[z] < min_overlap_this && used >= min_use ) usable[i] = False;
          ++used;    }

     // Test rule 2 and record extensions.

     c_next.clear( );
     while(1)
     {    int min_right_ext = 1000000000, min_right_ext_id = -1;
          
          // We proceed in reverse order to insure that we'll use the last read
          // (the terminator) if we can.

          for ( int i = usable.isize( ) - 1; i >= 0; i-- )
          {    if ( !usable[i] ) continue;
               int right_ext = xaligns[start+i].RightExt2( );
               if ( right_ext < min_right_ext )
               {    min_right_ext = right_ext, min_right_ext_id = i;    }    }
          if ( min_right_ext_id < 0 ) break;
          const xalign& a = xaligns[ start + min_right_ext_id ];
          usable[min_right_ext_id] = False;
          align_chain c_new(c);
          c_new.AddAlign(a);
          c_next.push_back(c_new);

          // Find and exclude extensions of a.

          for ( int i = 0; i < usable.isize( ); i++ )
          {    if ( !usable[i] ) continue;
               const xalign& b = xaligns[start+i];
               for ( int j = xaligns_index[a.id2]; j < xaligns_index[a.id2+1]; j++ )
               {    const xalign& c = xaligns[j];
                    if ( c.id2 != b.id2 ) continue;
                    if ( b.Offset( ) != a.Offset( ) + c.Offset( ) ) continue;
                    usable[i] = False;
                    break;    }    }    }    }

vecbasevector* xalign::reads;
vecbasevector* align_chain::reads;

void CombinePairs( vecbasevector& reads, vecqualvector& quals, 
     const vec< pair<int,int> >& pairs, const double score_prox, 
     const int min_score, const int min_qdiff, const Bool verbosity, String& report )
{
     // Look for overlaps between the two reads in a pair, and where
     // appropriate, merge them.  The code here is very similar to that
     // used by CloseUnipathGaps, in a somewhat different context.

     for ( int ip = 0; ip < pairs.isize( ); ip++ )
     {    if ( reads[2*ip].size( ) == 0 || reads[2*ip+1].size( ) == 0 ) continue;
          basevector &b1 = reads[2*ip], &b2 = reads[2*ip+1];
          qualvector &q1 = quals[2*ip], &q2 = quals[2*ip+1];
          int n1 = b1.size( ), n2 = b2.size( ), maxscore = -1000000000;
          int ceil_overlap = -( pairs[ip].first - 3*pairs[ip].second );
          if ( ceil_overlap <= 0 ) continue;
          int ceil2_overlap = -( pairs[ip].first - 6*pairs[ip].second );
          int max_overlap = Min( n1, n2, ceil2_overlap );
          static vec<int> bests, scores, goods, goodscores;
          bests.clear( ), scores.clear( ), goods.clear( ), goodscores.clear( );
          for ( int o = 0; o <= max_overlap; o++ )
          {    int plus = 0, minus = 0;
               for ( int j = 0; j < o; j++ )
               {    if ( b1[n1-o+j] == b2[j] ) plus++;
                    else minus++;    }
               int score = plus - minus;
               if ( score > maxscore ) maxscore = score;
               if ( score >= double(maxscore) / score_prox ) 
               {    bests.push_back(o);
                    scores.push_back(score);    }    }
          for ( int j = 0; j < bests.isize( ); j++ )
          {    if ( scores[j] >= double(maxscore) / score_prox )
               {    goods.push_back( bests[j] );
                    goodscores.push_back( scores[j] );    }    }
          if ( goods.solo( ) && goodscores[0] >= min_score 
               && goods[0] <= ceil_overlap )
          {    int o = goods[0];
               Bool too_close_to_call = False;
               for ( int i = 0; i < o; i++ )
               {    if ( b1[n1-o+i] != b2[i] )
                    {    int qdiff = Abs( q1[n1-o+i] - q2[i] );
                         if ( qdiff < min_qdiff ) too_close_to_call = True;    }    }
               if (too_close_to_call) continue;
               if ( verbosity >= 1 )
               {    report += "joining reads " + ToString(2*ip) + " and "
                         + ToString(2*ip+1) + " of a pair along "
                         + "overlap of size " + ToString(o) + "\n";    }
               basevector join(n1+n2-o);
               qualvector join_q(n1+n2-o);
               for ( int i = 0; i < n1; i++ )
               {    join.Set( i, b1[i] );
                    join_q[i] = q1[i];    }
               for ( int i = 0; i < n2; i++ )
               {    join.Set( i+n1-o, b2[i] );
                    join_q[i+n1-o] = q2[i];    }
               for ( int i = 0; i < o; i++ )
               {    if ( b1[n1-o+i] == b2[i] )
                         join_q[n1-o+i] = Max( q1[n1-o+i], q2[i] );
                    else
                    {    if ( verbosity >= 2 )
                         {    report += "base " + ToString(n1-o+i) + ", choosing "
                                   + "between Q scores " + ToString(int(q1[n1-o+i]))
                                   + " and " + ToString(int(q2[i])) + "\n";    }
                         if ( q1[n1-o+i] >= q2[i] )
                         {    join.Set( n1-o+i, b1[n1-o+i] );
                              join_q[n1-o+i] = q1[n1-o+i] - q2[i];    }
                         else join_q[n1-o+i] = q2[i] - q1[n1-o+i];    }    }
               b1 = join; b2.clear( );
               q1 = join_q; q2.clear( );    }    }    }

void CorrectErrors( const vecbasevector& reads, const vecqualvector& quals,
     const vec<xalign>& xaligns, const int min_align_length,
     vecbasevector& reads_new, vecqualvector& quals_new,
     const Bool log_correct, String& report )
{
     reads_new = reads;
     quals_new = quals;
     ostringstream* pout = 0;
     if (log_correct) pout = new ostringstream;
     for ( int i = 0; i < xaligns.isize( ); i++ )
     {    int j; 
          for ( j = i + 1; j < xaligns.isize( ); j++ )
               if ( xaligns[j].id1 != xaligns[i].id1 ) break;
          int64_t id1 = xaligns[i].id1;
          const basevector& r1 = reads[id1];
          int indicators = 0;
          const int max_ext_print = 10;
          if (log_correct)
          {    *pout << "\nerror correction of read " << id1 << "\n";
               int max_left = 0;
               for ( int l = i; l < j; l++ )
                    max_left = Max( max_left, xaligns[l].LeftExt2( ) );
               max_left = Min( max_left, max_ext_print );
               *pout << String( max_left + 5, ' ' );
               for ( int l = 0; l < r1.isize( ); l++ )
               {    if ( l % 10 == 0 ) *pout << ( l/10 % 10 );
                    else *pout << " ";    }
               *pout << "\n" << String( max_left + 5, ' ' );
               for ( int l = 0; l < r1.isize( ); l++ )
                    *pout << ( l % 10 );
               *pout << "\n" << String( 4 - ToString(id1).isize( ), ' ' ) 
                    << id1 << " ";
               *pout << String( max_left, '.' ) << r1.ToString( ) << "\n";
               for ( int l = i; l < j; l++ )
               {    const xalign& a = xaligns[l];
                    if ( a.Overlap( ) < min_align_length ) continue;
                    int64_t id2 = a.id2;
                    *pout << String( 4 - ToString(id2).isize( ), ' ' ) << id2 << " ";
                    *pout << String( Max( 0, max_left + a.Offset( ) ), '.' );
                    for ( int u = 0; u < reads[id2].isize( ); u++ )
                    {    if ( max_left + a.Offset( ) + u < 0 ) continue;
                         if ( u + a.Offset( ) >= r1.isize( ) + max_ext_print ) break;
                         *pout << as_base( reads[id2][u] );    }
                    *pout << "\n";    }
               *pout << "\n";    }
          for ( int p1 = 0; p1 < r1.isize( ); p1++ )
          {    vec<int> score( 4, 0 ), ids( 4, vec<int>::IDENTITY );
               score[ r1[p1] ] += quals[id1][p1];
               for ( int l = i; l < j; l++ )
               {    if ( xaligns[l].Overlap( ) < min_align_length ) continue;
                    int64_t id2 = xaligns[l].id2;
                    int p2 = p1 + xaligns[l].pos2 - xaligns[l].pos1;
                    if ( p2 < 0 || p2 >= reads[id2].isize( ) ) continue;
                    score[ reads[id2][p2] ] += quals[id2][p2];    }
               ReverseSortSync( score, ids );
               const int max_score_to_change = 200;
               if ( score[1] > max_score_to_change ) indicators++;
               if ( ids[0] != r1[p1] && score[1] <= max_score_to_change )
               {    reads_new[id1].Set( p1, ids[0] );
                    if (log_correct)
                    {    int t;
                         for ( t = 1; t < 4; t++ )
                              if ( ids[t] == r1[p1] ) break;
                         *pout << "read " << id1 << ", position " << p1 << ", " 
                              << as_base( r1[p1] ) << "[" << score[t] << "=" 
                              << int(quals[id1][p1]);
                         for ( int l = i; l < j; l++ )
                         {    int64_t id2 = xaligns[l].id2;
                              int p2 = p1 + xaligns[l].pos2 - xaligns[l].pos1;
                              if ( p2 < 0 || p2 >= reads[id2].isize( ) ) continue;
                              if ( reads[id2][p2] != r1[p1] ) continue;
                              *pout << "+" << int(quals[id2][p2]);    }
                         *pout << "] --> " << as_base( ids[0] )
                              << "[" << score[0] << "]\n";    }
                    quals_new[id1][p1] = 0;    }    }
          if (log_correct) *pout << "indicators = " << indicators << "\n";
          i = j - 1;    }    
     if (log_correct) 
     {    report += pout->str( );
          delete pout;    }     }
