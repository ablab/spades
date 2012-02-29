///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// AlignReads.  This is experimental code to align reads to an assembly.
// Generates as output: sub_dir + "/" + ASSEMBLY + ".readlocs".
//
// TO DO:
// 1. Use kmer class.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "Basevector.h"
#include "Equiv.h"
#include "MainTools.h"
#include "PairsManager.h"
#include "ParseSet.h"
#include "PrintAlignment.h"
#include "Superb.h"
#include "VecUtilities.h"
#include "pairwise_aligners/SmithWatBandedA.h"
#include "paths/ReadLoc.h"
#include "random/Shuffle.h"
#include "system/SortInPlace.h"

// Horrors.  A char20 is really 20 bases.

class char20 {
     public:
     unsigned char x[20/4];
     friend Bool operator==( const char20& I1, const char20& I2 )
     {    return memcmp( I1.x, I2.x, 20/4 ) == 0;    }
     friend Bool operator!=( const char20& I1, const char20& I2 )
     {    return memcmp( I1.x, I2.x, 20/4 ) != 0;    }
     friend Bool operator<( const char20& I1, const char20& I2 )
     {    return memcmp( I1.x, I2.x, 20/4 ) < 0;    }
     friend Bool operator>( const char20& I1, const char20& I2 )
     {    return memcmp( I1.x, I2.x, 20/4 ) > 0;    }
     void GetFromChars( const vec<char>& b )
     {    for ( int i = 0; i < 20/4; i++ )
               x[i] = 0;
          for ( int i = 0; i < 20; i++ )
          {    if ( b[i] == 'A' );
               else if ( b[i] == 'C' ) x[i/4] ^= ( 1 << 2*(i%4) );
               else if ( b[i] == 'G' ) x[i/4] ^= ( 2 << 2*(i%4) );
               else if ( b[i] == 'T' ) x[i/4] ^= ( 3 << 2*(i%4) );
               else
               {    cout << "GetFromChars: illegal character " << b[i] << endl;
                    cout << "Abort." << endl;
                    exit(1);    }    }    }
               
};

typedef char20 FASTAVECTOR;

class kmer_pos_fw_id {

     public:

     kmer_pos_fw_id( ) : pos(0), fw(False), id(0) { }

     kmer_pos_fw_id( const FASTAVECTOR& kmer, const int pos, const Bool fw,
          const size_t id ) : kmer(kmer), pos(pos), fw(fw), id(id) { }

     FASTAVECTOR kmer;
     int pos;
     Bool fw;
     size_t id;

     friend Bool operator<( const kmer_pos_fw_id& x1, const kmer_pos_fw_id& x2 )
     {    if ( x1.kmer < x2.kmer ) return True;
          if ( x2.kmer < x1.kmer ) return False;
          if ( x1.id < x2.id ) return True;
          if ( x1.id > x2.id ) return False;
          if ( x1.fw < x2.fw ) return True;
          if ( x1.fw > x2.fw ) return False;
          if ( x1.pos < x2.pos ) return True;
          return False;    }

};

class pos_fw_id {

     public:

     pos_fw_id( ) : pos(0), fw(False), id(0) { }

     pos_fw_id( const int pos, const Bool fw,
          const size_t id ) : pos(pos), fw(fw), id(id) { }

     int pos;
     Bool fw;
     size_t id;

     friend Bool operator==( const pos_fw_id& x1, const pos_fw_id& x2 )
     {    return x1.pos == x2.pos && x1.fw == x2.fw && x1.id == x2.id;    }

     friend Bool operator<( const pos_fw_id& x1, const pos_fw_id& x2 )
     {    if ( x1.id < x2.id ) return True;
          if ( x1.id > x2.id ) return False;
          if ( x1.fw < x2.fw ) return True;
          if ( x1.fw > x2.fw ) return False;
          if ( x1.pos < x2.pos ) return True;
          return False;    }

};

String SingleBar( )
{    return "----------------------------------------------"
          "--------------------------------------";    }
String DoubleBar( )
{    return "=============================================="
          "======================================";    }

// Test for inconsistency.  These tests are not exactly right and should be refined.

Bool Inconsistent( const int tig1, const int tig2, const Bool FW1, const Bool FW2,
     const int offset1, const int offset2, const vecbasevector& tigs,
     const vec<superb>& scaffolds,
     const vec<int>& to_super, const vec<int>& to_super_pos )
{    int s1 = to_super[tig1], s2 = to_super[tig2];
     const superb &S1 = scaffolds[s1], &S2 = scaffolds[s2];
     int p1 = to_super_pos[tig1], p2 = to_super_pos[tig2];
     int L1 = tigs[tig1].size( ), L2 = tigs[tig2].size( );
     if ( s1 == s2 && FW1 == FW2 ) return True;
     if ( tig1 == tig2 && Abs( offset1 - offset2 ) > 500 ) return True;
     if ( tig1 != tig2 )
     {    if ( s1 != s2 ) 
          {    
               // In this case we require that we're at the end of one scaffold and
               // at the end of a contig on the other scaffold.

               Bool end_S1 = ( p1 == 0 && offset1 < 100 )
                    || ( p1 == S1.Ntigs( ) - 1 && offset1 >= L1 - 200 );
               Bool end_S2 = ( p2 == 0 && offset2 < 100 )
                    || ( p2 == S2.Ntigs( ) - 1 && offset2 >= L2 - 200 );
               Bool end_tig1 = ( offset1 < 100 || offset1 >= L1 - 200 );
               Bool end_tig2 = ( offset2 < 100 || offset2 >= L2 - 200 );
               if ( !( (end_S1 && end_tig2) || (end_S2 && end_tig1) ) )
                    return True;    }
          else if ( Abs( p1 - p2 ) != 1 ) return True;    }
     return False;     }

class align_plus {

     public:

     align_plus( ) : tig(-1) { }
     align_plus( const align& a, const int tig, const Bool fw )
          : a(a), tig(tig), fw(fw) { }

     align a;
     int tig;
     Bool fw;

};

vec<int> *TO_SUPER, *TO_SUPER_POS;

void PrintAlignmentNicely( const align_plus& ap, const basevector& TIG, 
     const vecbasevector& bases, const int64_t id )
{    ostringstream sout;
     const align& a = ap.a;
     cout << "id = " << id << ", tig = " << ap.tig << " (" << (*TO_SUPER)[ap.tig] 
          << "." << (*TO_SUPER_POS)[ap.tig] << "), " 
          << ( ap.fw ? "fw" : "rc" ) << "\n";
     sout << "read: " << a.pos1( ) << " - " << a.Pos1( ) << ", contig: " 
          << a.pos2( ) << " - " << a.Pos2( ) << " of " << TIG.size( ) << "\n";
     basevector b = bases[id];
     if ( !ap.fw ) b.ReverseComplement( );
     PrintVisualAlignment( False, sout, b, TIG, a );
     String vis = sout.str( );
     for ( int i = 0; i < 80; i++ )
          vis.GlobalReplaceBy( " \n", "\n" );
     vis.GlobalReplaceBy( "\n\n\n", "\n\n" );
     vis.GlobalReplaceBy( "\n\n\n", "\n\n" );
     if ( vis.Contains( "\n\n", -1 ) ) vis.resize( vis.isize( ) - 1 );
     cout << vis;    
     flush(cout);    }

void HashReadPair( const int K, 
                   const int64_t id1, 
                   const int64_t id2, 
                   const vecbasevector& bases, 
                   const vec<int64_t>& hash, 
                   const vec<char20>& kmers,
                   const vec<int32_t>& len, 
                   vec< vec<size_t> >& starts, 
                   vec< vec<size_t> >& stops,
                   vec< vec<kmer_pos_fw_id> >& kpfir )
{
     int64_t nhash = hash.size( );
     starts.resize(2), stops.resize(2), kpfir.resize(2);
     for ( int side = 0; side < 2; side++ )
     {    int64_t id = ( side == 0 ? id1 : id2 );
          int nkmers = bases[id].isize( ) - K + 1;
          if (nkmers <= 0) continue;
          starts[side].resize(nkmers), stops[side].resize(nkmers);
          kpfir[side].resize(nkmers);
          basevector b, brc;
          vec<char> akmer(K);
          FASTAVECTOR kmer;
          for ( int j = 0; j < nkmers; j++ )
          {    
               // Compute kmer.

               b.SetToSubOf( bases[id], j, K );
               brc = b;
               brc.ReverseComplement( );
               basevector& ab = ( b < brc ? b : brc );
               for ( int l = 0; l < K; l++ )
                    akmer[l] = as_base( ab[l] );
               kmer.GetFromChars(akmer);
               kpfir[side][j] = kmer_pos_fw_id(kmer, j, b < brc ? True : False, id);

               // Locate kmer.

               int64_t x = 0;
               for ( int l = 0; l < 5; l++ )
               {    x += int64_t( kmer.x[l] );
                    x = x << 8;    }
               int64_t key = x % nhash;
               int64_t h = hash[key];

               while(1)
               {    if ( h < 0 )
                    {    starts[side][j] = 0;
                         stops[side][j] = 0;
                         break;    }
                    if ( kmers[h] == kmer )
                    {    starts[side][j] = h;
                         stops[side][j] = h + len[h];
                         break;    }
                    key++;
                    if ( key == nhash ) key = 0;
                    h = hash[key];    }    }    }    }

void PrintCounts( const int side, const vec< vec<size_t> >& starts,
     const vec< vec<size_t> >& stops )
{    for ( int j = 0; j < starts[side].isize( ); j++ )
     {    if ( j > 0 ) cout << ",";
          int k;
          for ( k = j+1; k < starts[side].isize( ); k++ )
          {    if ( stops[side][k] - starts[side][k] !=
                    stops[side][j] - starts[side][j] )
               {    break;    }    }
          cout << stops[side][j] - starts[side][j];
          if ( k - j > 1 ) cout << "[" << k-j << "]";
          j = k - 1;    }
     cout << "\n";    }

void DumpAlignment( const vecbasevector& bases, const int64_t id,
     const align_plus& ap, const vecbasevector& tigs, const int side, 
     const vec< vec<size_t> >& starts, const vec< vec<size_t> >& stops,
     const vec<pos_fw_id>& places, const String& title )
{
     #pragma omp critical
     {    cout << "\n" << DoubleBar( ) << "\n" << title << ":\n";
          cout << bases[id].ToString( ) << "\n";
          PrintCounts( side, starts, stops );
          PRINT( places.size( ) );
          PrintAlignmentNicely( ap, tigs[ap.tig], bases, id );
          cout << DoubleBar( ) << "\n\n";    }    }

double Align( const basevector& b, const basevector& TIG, const int o1, 
     const int o2, const int tig, const Bool fw, align_plus& ap )
{    int errors;
     const int bandwidth_add = 5;
     int o = (o1 + o2)/2;
     SmithWatBandedA( b, TIG, -o, (o2 - o1)/2 + bandwidth_add, ap.a, errors );
     ap.tig = tig;
     ap.fw = fw;
     return double( ActualErrors( b, TIG, ap.a ) ) / double( ap.a.extent1( ) );    }

void IncrementCounter( int64_t& count )
{   
     #pragma omp critical
     count++;    }

void GetPlaces( const size_t max_freq, const int n, const int K, const int side,
     const vec< vec<size_t> >& starts, const vec< vec<size_t> >& stops,
     const vec< vec<kmer_pos_fw_id> >& kpfir, const vec<pos_fw_id>& pfi,
     vec< vec<pos_fw_id> >& places )
{
     int nkmers = n - K + 1;
     for ( int j = 0; j < nkmers; j++ )
     {    if ( stops[side][j] - starts[side][j] <= max_freq )
          {    for ( size_t u = starts[side][j]; u < stops[side][j]; u++ )
               {    int tig = pfi[u].id;
                    Bool fw = !( kpfir[side][j].fw ^ pfi[u].fw );
                    int offset;
                    if (fw) offset = pfi[u].pos - kpfir[side][j].pos;
                    else offset = kpfir[side][j].pos - ( n - pfi[u].pos - K );
                    places[side].push( offset, fw, tig );    }    }    }
     UniqueSort(places[side]);    }

void FindAcceptableAlignment( const int rpass, const double max_OK, 
     const vecbasevector& bases, const uint64_t id, const vecbasevector& tigs, 
     const int K, const int side, const vec< vec<size_t> >& starts, 
     const vec< vec<size_t> >& stops, const vec< vec<kmer_pos_fw_id> >& kpfir, 
     const vec<pos_fw_id>& pfi, vec< vec<pos_fw_id> >& places, 
     align_plus& ap, double& BADNESS )
{
     // Find the placements of the read that are defined by unique placements of a
     // read kmer on the assembly.  These placements have the form (offset, fw, tig)
     // where the read is placed on contig tig starting at position offset, with 
     // orientation fw.

     const size_t max_freq = 1;
     GetPlaces( max_freq, bases[id].size( ), K, side, starts, stops,
          kpfir, pfi, places );

     // Divide these places into groups.

     const int group_sep_bound = 100;
     vec<align_plus> aligns;
     vec<double> badnesses;
     for ( int i = 0; i < places[side].isize( ); i++ )
     {    int j;
          for ( j = i + 1; j < places[side].isize( ); j++ )
          {    if ( places[side][j].id != places[side][i].id ) break;
               if ( places[side][j].fw != places[side][i].fw ) break;
               if ( places[side][j].pos > places[side][i].pos + group_sep_bound ) 
                    break;    }

          // Carry out one alignment for each group of places.

          int o1 = places[side][i].pos, o2 = places[side][j-1].pos;
          int o = o1 + (o2-o1)/2;
          int tig = places[side][i].id;
          Bool fw = places[side][i].fw;
          basevector b = bases[id];
          if ( !fw ) b.ReverseComplement( );
          align_plus ap;
          double badness = Align( b, tigs[tig], o1, o2, tig, fw, ap );
          aligns.push_back(ap);
          badnesses.push_back(badness);
          i = j - 1;    }

     // Sort alignments.

     if ( aligns.empty( ) ) return; /* fail */
     SortSync( badnesses, aligns );

     // Special processing for long jumps.  We require that the first 10 bases align
     // with no errors, and that the first 20 bases align with at most four errors.
     // It would probably be better to compare to fastavector contigs as at present
     // this will lose reads landing on SNPs.  Note that the algorithm relies on the
     // fact that long jump reads are not reversed.  If the molecular biology 
     // changed this would need to be adjusted.

     if ( rpass == 2 )
     {    Bool long_jump_OK = False;
          const align_plus& ap = aligns[0];
          align a = ap.a;
          const basevector& TIG = tigs[ap.tig];
          basevector b = bases[id];
          if ( !ap.fw ) b.ReverseComplement( );
          if ( !ap.fw ) a.ReverseThis( b.size( ), TIG.size( ) );
          const int min_bases_perf = 10;
          const int min_bases = 20;
          const int max_errors = 4;
          if ( a.pos1( ) == 0 && a.Pos1( ) >= min_bases )
          {    long_jump_OK = True;
               int p1 = a.pos1( ), p2 = a.pos2( ), errs = 0;
               for ( int j = 0; j < a.Nblocks( ); j++ )
               {    if ( p1 >= min_bases ) break;
                    if ( a.Gaps(j) > 0 )
                    {    p2 += a.Gaps(j);
                         errs += a.Gaps(j);    }
                    if ( a.Gaps(j) < 0 )
                    {    p1 -= a.Gaps(j);
                         errs -= a.Gaps(j);    }
                    if ( errs > max_errors || ( p1 < min_bases_perf && errs > 0 ) )
                    {    long_jump_OK = False;
                         break;    }
                    for ( int x = 0; x < a.Lengths(j); x++ )
                    {    if ( p1 >= min_bases ) break;
                         char c1 = ( ap.fw ? b[p1] : 3 - b[ b.isize( ) - p1 - 1 ] );
                         char c2 = ( ap.fw ? TIG[p2]
                              : 3 - TIG[ TIG.isize( ) - p2 - 1 ] );
                         if ( c1 != c2 ) errs++;
                         if ( errs > max_errors
                              || ( p1 < min_bases_perf && errs > 0 ) )
                         {    long_jump_OK = False;
                              break;    }
                         ++p1; ++p2;    }
                    if ( !long_jump_OK ) break;    }    }
          if ( !long_jump_OK ) return; /* fail */    }

     // Decide if best alignment is good enough.

     if ( badnesses[0] > max_OK ) return; /* fail */
     if ( badnesses.size( ) >= 2 && badnesses[0] == badnesses[1] ) return; /* fail */

     // More special processing for jumps and long jumps....
     // Test to see if the alignment is suspicious.  This uses some half-ass
     // heuristics that seem to be good enough:
     // (1) We do nothing unless the number of read kmers having multiplicity > 1 in
     //     the genomes if >= the number having  multiplicity 1, minus 3.
     // (2) Otherwise, we look for read kmers that have multiplicity > 1 in the
     //     genome but are not aligned perfectly.  If that happens we call the
     //     alignment suspicious.
     // Note that this does lose some good alignments.  Probably the way to do this 
     // correctly is to actually check the alignments for read kmers having
     // multiplicity below some reasonable threshold.

     if ( rpass == 1 || rpass == 2 )
     {    const align_plus& ap = aligns[0];
          Bool suspicious = False;
          int count1 = 0, count2p = 0;
          int nkmers = bases[id].isize( ) - K + 1;
          for ( int j = 0; j < nkmers; j++ )
          {    int freq = stops[side][j] - starts[side][j];
               if ( freq == 1 ) count1++;
               else if ( freq > 1 ) count2p++;    }
          const int count_delta = 3;
          if ( count2p >= count1 - count_delta )
          {    basevector b = bases[id];
               if ( !ap.fw ) b.ReverseComplement( );
               const basevector& TIG = tigs[ap.tig];
               int p1 = ap.a.pos1( ), p2 = ap.a.pos2( );
               vec<Bool> matches( nkmers, False );
               for ( int j = 0; j < ap.a.Nblocks( ); j++ )
               {    if ( ap.a.Gaps(j) > 0 ) p2 += ap.a.Gaps(j);
                    if ( ap.a.Gaps(j) < 0 ) p1 -= ap.a.Gaps(j);
                    int start = 0, p1_begin = p1;
                    for ( int x = 0; x < ap.a.Lengths(j); x++ )
                    {    if ( b[p1] != TIG[p2] )
                         {    for ( int l = start; l <= x - K; l++ )
                                   matches[ p1_begin + l ] = True;
                              start = x + 1;    }
                         ++p1; ++p2;    }    
                    for ( int l = start; l <= ap.a.Lengths(j) - K; l++ )
                         matches[ p1_begin + l ] = True;    }
               for ( int j = 0; j < nkmers; j++ )
               {    if ( stops[side][j] - starts[side][j] > max_freq )
                    {    
                         // Does kmer starting at j match perfectly in the
                         // alignment?
     
                         if ( !matches[j] )
                         {    suspicious = True;
                              break;    }    }    }    }
          if (suspicious)
          {    // #pragma omp critical
               // cout << "alignment of read " << id << " is suspicious" << endl;
               return; /* fail */   }    }

     // Accept alignment.

     ap = aligns[0];
     BADNESS = badnesses[0];    }

void ProcessPairs(

     // pair and read data

     const int rpass,
     const PairsManager& pairs,
     const vecbasevector& bases,

     // assembly data

     const vecbasevector& tigs,
     const vec<superb>& scaffolds,
     const vec<int>& to_super,
     const vec<int>& to_super_pos,

     // hash info

     const int K,
     const vec<int64_t>& hash,
     const vec<char20>& kmers,
     const vec<pos_fw_id>& pfi,
     const vec<int32_t>& len,

     // logging control

     const Bool PRINT_STATUS,
     const Bool PRINT_STATS,
     const vec<int>& trace_ids,
     const Bool PRINT_BADS,
     const Bool PRINT_ALIGNMENTS,
     const int PRINT_ALIGNMENTS_TIG,
     int& PRINT_MISSING,
     const Bool PRINT_INCONSISTENT,

     // results

     vec<read_loc>& locs )
{
     // Track stats.

     int64_t neither_end = 0, one_end = 0;
     int64_t both_ends_consistently = 0, both_ends_inconsistently = 0;
     int64_t both_ends_after = 0;
     int64_t inconsistent_fixed = 0;

     // Go through pairs.

     if (PRINT_STATUS) 
     {    cout << Date( ) << ": going through pairs, mem usage = "
               << ToStringAddCommas( MemUsageBytes( ) ) << endl;    }
     double clock = WallClockTime( );
     const size_t batch_size = 10000;
     const size_t nbatches = ( pairs.nPairs( ) + batch_size - 1 ) / batch_size;
     size_t batches_processed = 0, dots_printed = 0;
     if (PRINT_STATUS)
     {    cout << Date( ) << ": ";
          flush(cout);    }
     #pragma omp parallel for
     for ( size_t batch = 0; batch < nbatches; batch++ )
     {    vec< vec<size_t> > starts, stops;
          vec< vec<kmer_pos_fw_id> > kpfir;
          vec< vec<pos_fw_id> > places(2);
          align_plus ap1, ap2;

          // <-----

     for ( size_t pi = batch * batch_size; 
          pi < Min( (batch+1) * batch_size, pairs.nPairs( ) ); pi++ )
     {    int64_t id1 = pairs.ID1(pi), id2 = pairs.ID2(pi);
          if ( trace_ids.nonempty( ) && !BinMember( trace_ids, id1 )
               && !BinMember( trace_ids, id2 ) )
          {    continue;    }
          int offset1 = 0, offset2 = 0;
          ap1.tig = -1, ap2.tig = -1;

          // Process pair.

          starts.clear( ), stops.clear( ), kpfir.clear( );
          places[0].clear( ), places[1].clear( );
          HashReadPair( K, id1, id2, bases, hash, kmers, len, starts, stops, kpfir );
          double max_OK = 0.5;
          if ( rpass == 2 ) max_OK = 0.25;
          for ( int side = 0; side < 2; side++ )
          {    align_plus& ap = ( side == 0 ? ap1 : ap2 );
               int64_t id = ( side == 0 ? id1 : id2 );
               int nkmers = bases[id].isize( ) - K + 1;
               if (nkmers <= 0) continue;
               double BADNESS;
               FindAcceptableAlignment( rpass, max_OK, bases, id, tigs, K, side, 
                    starts, stops, kpfir, pfi, places, ap, BADNESS );
               if ( ap.tig >= 0 )
               {    ( side == 0 ? offset1 : offset2 ) = -ap.a.Offset( );
                    if (PRINT_BADS)
                    {    const double max_good = 0.3;
                         if ( BADNESS > max_good )
                         {    DumpAlignment( bases, id, ap, tigs, side, starts,
                                   stops, places[side], 
                                   "BAD ALIGNMENT DETECTED" );    }    }
                    if ( PRINT_ALIGNMENTS || PRINT_ALIGNMENTS_TIG == ap.tig )
                    {    DumpAlignment( bases, id, ap, tigs, side, starts, stops,
                              places[side], "ALIGNMENT" );    }    }
               else if ( PRINT_MISSING > 0 )
               {
                    #pragma omp critical
                    {    cout << "\n" << DoubleBar( ) << "\n";
                         cout << "MISSING ALIGNMENT DETECTED:\n";
                         cout << id << ", " << bases[id].ToString( ) << "\n";
                         PrintCounts( side, starts, stops );
                         cout << "\n" << DoubleBar( ) << "\n\n";    
                         PRINT_MISSING--;    }    }    }

          // Print inconsistent results.

          Bool inconsistent = False;
          int tig1 = ap1.tig, tig2 = ap2.tig;
          Bool FW1 = ap1.fw, FW2 = ap2.fw;
          const align &a1 = ap1.a, &a2 = ap2.a;
          if ( tig1 >= 0 && tig2 >= 0 )
          {    inconsistent = Inconsistent( tig1, tig2, FW1, FW2, offset1, offset2, 
                    tigs, scaffolds, to_super, to_super_pos );    }
          if ( inconsistent && PRINT_INCONSISTENT )
          {    
               #pragma omp critical
               {    cout << "\n" << DoubleBar( ) << "\n";
                    cout << "INCONSISTENCY DETECTED (pair = " << pi << "):\n";
                    for ( int side = 0; side < 2; side++ )
                    {    if ( side == 1 ) cout << SingleBar( ) << "\n";
                         int64_t id = ( side == 0 ? id1 : id2 );
                         const align_plus& ap = ( side == 0 ? ap1 : ap2 );
                         cout << bases[id].ToString( ) << "\n";
                         PrintCounts( side, starts, stops );
                         PrintAlignmentNicely( ap, tigs[ap.tig], bases, id );    }
                    cout << DoubleBar( ) << "\n\n";    }    }

          // Record results.

          Bool left_end_only = ( tig1 >= 0 && tig2 < 0 );
          Bool right_end_only = ( tig1 < 0 && tig2 >= 0 );
          #pragma omp critical
          {    if ( tig1 < 0 && tig2 < 0 ) neither_end++;
               else if ( left_end_only || right_end_only ) one_end++;
               else if (inconsistent) both_ends_inconsistently++;
               else both_ends_consistently++;
               if ( tig1 >= 0 && tig2 >= 0 ) both_ends_after++;    }

          // If only one end aligned, look for home for partner.  Also if pair is
          // placed inconsistently, try to rescue.

          align_plus ap1_new, ap2_new;
          for ( int pass = 0; pass < 2; pass++ )
          {    
               // On pass = 0, look for left read placement; on pass = 1, right.

               if ( pass == 0 && !( ( tig1 < 0 && tig2 >= 0 ) || inconsistent ) )
                    continue;
               if ( pass == 1 && !( ( tig1 >= 0 && tig2 < 0 ) || inconsistent ) )
                    continue;
               const align_plus& ap = ( pass == 0 ? ap2 : ap1 );
               int tig = ap.tig;
               const align& a = ap.a;
               int side = pass;
               int64_t id = ( pass == 0 ? id1 : id2 );
               align_plus& ap_new = ( pass == 0 ? ap1_new : ap2_new );

               // Define where the other read should land.  Low and high are bounds
               // for the predicted starting point of the other read on the contig.

               Bool fw = !ap.fw;
               const int dev_mult = 3;
               int low, high, sep = pairs.sep(pi), dev = pairs.sd(pi);
               if (ap.fw)
               {    low = a.Pos2( ) + sep - dev_mult * pairs.sd(pi);
                    high = a.Pos2( ) + sep + dev_mult * pairs.sd(pi);     }
               else
               {
                    low = a.pos2( ) - sep - bases[id].isize( ) - dev_mult * dev;
                    high = a.pos2( ) - sep - bases[id].isize( ) + dev_mult * dev;   }
               vec<int> offsets;
               for ( int j = 0; j < starts[side].isize( ); j++ )
               {    int64_t start = starts[side][j], stop = stops[side][j];
                    Bool fwx = ( fw == kpfir[side][j].fw );
                    int add = ( fw ? j : bases[id].isize( ) - j - K );
                    pos_fw_id LOW( low+add, fwx, tig ), HIGH( high+add, fwx, tig );
                    int64_t startx = lower_bound( pfi.begin( ) + start, 
                         pfi.begin( ) + stop, LOW ) - pfi.begin( );
                    int64_t stopx = upper_bound( pfi.begin( ) + start, 
                         pfi.begin( ) + stop, HIGH ) - pfi.begin( );
                    for ( int64_t u = startx; u < stopx; u++ )
                    {    if (fw) 
                              offsets.push_back( pfi[u].pos - kpfir[side][j].pos );
                         else 
                         {    offsets.push_back( kpfir[side][j].pos 
                                   - ( bases[id].size( ) 
                                   - pfi[u].pos - K ) );    }    }    }
               UniqueSort(offsets);

               // Align it.

               if ( offsets.nonempty( ) )
               {    basevector b = bases[id];
                    if ( !fw ) b.ReverseComplement( );
                    double badness = Align( b, tigs[tig], offsets.front( ), 
                         offsets.back( ), tig, fw, ap_new );
                    if ( badness > max_OK ) ap_new.tig = -1;    }    }

          // Save alignments.

          const basevector &b1 = bases[id1], &b2 = bases[id2];
          basevector b1rc(b1), b2rc(b2);
          b1rc.ReverseComplement( ), b2rc.ReverseComplement( );
          if ( ap1_new.tig >= 0 && ap2_new.tig < 0 )
          {    if ( tig1 < 0 ) IncrementCounter(both_ends_after);
               else IncrementCounter(inconsistent_fixed);
               ap1 = ap1_new;
               if ( PRINT_ALIGNMENTS || PRINT_ALIGNMENTS_TIG == ap1.tig )
               {    DumpAlignment( bases, id1, ap1, tigs, 0, starts, stops,
                         places[0], "RECOVERED ALIGNMENT" );    }    }
          if ( ap1_new.tig < 0 && ap2_new.tig >= 0 )
          {    if ( tig2 < 0 ) IncrementCounter(both_ends_after);
               else IncrementCounter(inconsistent_fixed);
               ap2 = ap2_new;
               if ( PRINT_ALIGNMENTS || PRINT_ALIGNMENTS_TIG == ap2.tig )
               {    DumpAlignment( bases, id2, ap2, tigs, 1, starts, stops,
                         places[1], "RECOVERED ALIGNMENT" );    }    }
          if ( ap1_new.tig >= 0 && ap2_new.tig >= 0 )
          {    int e1 = ActualErrors( ( ap1.fw ? b1 : b1rc ), tigs[ap1.tig], ap1.a );
               int e1_new = ActualErrors( ( ap1_new.fw ? b1 : b1rc ), 
                    tigs[ap1_new.tig], ap1_new.a );
               int e2 = ActualErrors( ( ap2.fw ? b2 : b2rc ), tigs[ap2.tig], ap2.a );
               int e2_new = ActualErrors( ( ap2_new.fw ? b2 : b2rc ), 
                    tigs[ap2_new.tig], ap2_new.a );
               double badness_a = double( e1 + e2_new ) 
                    / double( ap1.a.extent1( ) + ap2_new.a.extent1( ) );
               double badness_b = double( e1_new + e2 ) / 
                    double( ap1_new.a.extent1( ) + ap2.a.extent1( ) );
               if ( badness_b < badness_a )
               {    ap1 = ap1_new;
                    if ( PRINT_ALIGNMENTS || PRINT_ALIGNMENTS_TIG == ap1.tig )
                    {    DumpAlignment( bases, id1, ap1, tigs, 0, starts, stops,
                              places[0], "RESCUED ALIGNMENT" );    }    }
               if ( badness_a < badness_b )
               {    ap2 = ap2_new;
                    if ( PRINT_ALIGNMENTS || PRINT_ALIGNMENTS_TIG == ap2.tig )
                    {    DumpAlignment( bases, id2, ap2, tigs, 1, starts, stops,
                              places[1], "RESCUED ALIGNMENT" );    }    }
               if ( badness_a != badness_b )
                    IncrementCounter(inconsistent_fixed);    }

          // Record alignments.

          #pragma omp critical
          {    int library_id = pairs.libraryID(pi);
               if ( ap1.tig >= 0 ) 
               {    locs.push( ap1.a, id1, ap1.tig, ap1.fw,
                         rpass, library_id, bases[id1].size( ) );
                    locs.back( ).SetPartnerReadId(id2);
                    locs.back( ).SetPartnerReadLength( bases[id2].size( ) );    }
               if ( ap2.tig >= 0 ) 
               {    locs.push( ap2.a, id2, ap2.tig, ap2.fw,
                         rpass, library_id, bases[id2].size( ) );
                    locs.back( ).SetPartnerReadId(id1);
                    locs.back( ).SetPartnerReadLength( bases[id1].size( ) );    }
               if ( ap1.tig >= 0 && ap2.tig >= 0 )
               {    read_loc& rl1 = locs[locs.size( ) - 2]; 
                    read_loc& rl2 = locs[locs.size( ) - 1];
                    rl1.SetPartnerPlaced( ); rl2.SetPartnerPlaced( );
                    rl1.SetPartnerContigId(ap2.tig); rl2.SetPartnerContigId(ap1.tig);
                    rl1.SetPartnerFw(ap2.fw); rl2.SetPartnerFw(ap1.fw);
                    rl1.SetPartnerStart( rl2.Start( ) );
                    rl2.SetPartnerStart( rl1.Start( ) );
                    rl1.SetPartnerBandwidth( rl2.Bandwidth( ) );
                    rl2.SetPartnerBandwidth( rl1.Bandwidth( ) );    }    }    }

     if (PRINT_STATUS)
     {
          #pragma omp critical
          {    batches_processed++;
               double done_percent 
                    = 100.0 * double(batches_processed) / double(nbatches);
               while ( done_percent >= dots_printed+1)
               {    if ( dots_printed % 10 == 0 && dots_printed > 0 
                         && dots_printed != 50 )
                    {    cout << " ";    }
                    cout << ".";
                    dots_printed++;
                    if ( dots_printed % 50 == 0 ) cout << endl;    
                    if ( dots_printed == 50 ) cout << Date( ) << ": ";
                    flush(cout);    }    }    }
     }

     // Report results.

     if (PRINT_STATS)
     {    int64_t np = pairs.nPairs( );
          cout << "\nPAIR STATISTICS\n";
          cout << "neither end aligned = " << neither_end << " (" 
               << PERCENT_RATIO( 3, neither_end, np ) << ")" << endl;
          cout << "one end aligned = " << one_end << " (" 
               << PERCENT_RATIO( 3, one_end, np ) << ")" << endl;
          cout << "both ends aligned consistently = " << both_ends_consistently 
               << " (" << PERCENT_RATIO( 3, both_ends_consistently, np ) 
               << ")" << endl;
          cout << "both ends aligned inconsistently = " 
               << both_ends_inconsistently << " (" 
               << PERCENT_RATIO( 3, both_ends_inconsistently, np ) << ")" << endl;
          cout << "both ends aligned consistently after fixing = " 
               << both_ends_after << " (" 
               << PERCENT_RATIO( 3, both_ends_after, np ) << ")" << endl;
          cout << "inconsistent pairs fixed = " << inconsistent_fixed << " (" 
               << PERCENT_RATIO( 3, inconsistent_fixed, np ) << ")\n" << endl;    }

     if (PRINT_STATUS) 
     {    cout << Date( ) << ": " << TimeSince(clock) 
               << " used going through reads" << endl;    }    }

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_String_OrDefault(SUBDIR, "test");
     CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
	   "Number of threads to use (use all available processors if set to 0)");
     CommandArgument_String_OrDefault(ASSEMBLY, "linear_scaffolds0.clean.patched");
     CommandArgument_Bool_OrDefault(PRINT_ALIGNMENTS, False);
     CommandArgument_Int_OrDefault(PRINT_ALIGNMENTS_TIG, -1);
     CommandArgument_Bool_OrDefault(PRINT_INCONSISTENT, False);
     CommandArgument_Int_OrDefault(PRINT_MISSING, 0);
     CommandArgument_Bool_OrDefault(PRINT_BADS, False);
     CommandArgument_Bool_OrDefault(PRINT_STATUS, True);
     CommandArgument_Bool_OrDefault(PRINT_STATS, True);
     CommandArgument_String_OrDefault_Doc(TRACE_IDS, "",
          "just trace what happens to these reads and their partners; "
          "automatically turns off writing");
     CommandArgument_String_OrDefault_Doc(READLOCS_PREFIX, "",
          "if specified, file extension is .READLOCS_PREFIX.readlocs "
          "instead of .readlocs");
     CommandArgument_String_OrDefault(LIB_TYPES, "{frag,jump,long_jump}");
     // Allow the program to accept read files input other
     // than the default hard coded files:
     // frag_reads_filt_cpd.fastb 
     // jump_reads_filt_cpd.fastb 
     // long_jump_reads_filt.fastb 
     CommandArgument_String_OrDefault(FRAG_READS, "frag_reads_filt_cpd");
     CommandArgument_String_OrDefault(JUMP_READS, "jump_reads_filt_cpd");
     CommandArgument_String_OrDefault(LONG_JUMP_READS, "long_jump_reads_filt");
     // to generate alignments of the scaffold_reads
     CommandArgument_Int_OrDefault(MAX_ASSEMBLY, -1);
     EndCommandArguments;

     // Thread control

     NUM_THREADS = configNumThreads(NUM_THREADS);
     omp_set_num_threads( NUM_THREADS );

     // Parse arguments.

     double begin_clock = WallClockTime( );
     vec<int> trace_ids;
     ParseIntSet( TRACE_IDS, trace_ids );
     vec<Bool> lib_types( 3, False );
     vec<String> lib_types_s;
     ParseStringSet( LIB_TYPES, lib_types_s );
     for ( int i = 0; i < lib_types_s.isize( ); i++ )
     {    if ( lib_types_s[i] == "frag" ) lib_types[0] = True;
          else if ( lib_types_s[i] == "jump" ) lib_types[1] = True;
          else if ( lib_types_s[i] == "long_jump" ) lib_types[2] = True;
          else
          {    cout << "Illegal LIB_TYPES." << endl;
               return 1;    }    }

     // Define directories.
  
     String data_dir = PRE + "/" + DATA, run_dir = PRE + "/" + DATA + "/" + RUN;
     String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR; 

     // Load assembly.

     if (PRINT_STATUS) cout << Date( ) << ": loading contigs" << endl;
     vecbasevector tigs( sub_dir + "/" + ASSEMBLY + ".contigs.fastb" );
     uint32_t ntigs = tigs.size( );
     vec<superb> scaffolds;
     ReadSuperbs( sub_dir + "/" + ASSEMBLY + ".superb", scaffolds );
     vec<int> to_super( ntigs, -1 ), to_super_pos( ntigs, -1 );
     for ( int i = 0; i < scaffolds.isize( ); i++ ) 
     {    for ( int j = 0; j < scaffolds[i].Ntigs( ); j++ ) 
          {    to_super[ scaffolds[i].Tig(j) ] = i;
               to_super_pos[ scaffolds[i].Tig(j) ] = j;    }    }
     TO_SUPER = &to_super, TO_SUPER_POS = &to_super_pos;

     // Check for oversized assembly.

     String head = sub_dir + "/" + ASSEMBLY;
     if ( READLOCS_PREFIX != "" ) head += "." + READLOCS_PREFIX;
     if ( MAX_ASSEMBLY >= 0 )
     {    int64_t total_size = 0;
          for ( size_t i = 0; i < tigs.size( ); i++ )
               total_size += tigs[i].size( );
          if ( total_size > MAX_ASSEMBLY )
          {    cout << "Assembly size exceeds MAX_ASSEMBLY, writing empty file." 
                    << endl;
               vec<read_loc> locs;
               WriteReadLocs( head, tigs.size( ), locs );
               return 0;    }    }

     // Create kmers in contigs.
     // (this has been ribeiro-fied)
     if (PRINT_STATUS) cout << Date( ) << ": start building kpfi" << endl;
     int K = 20;
     vec<kmer_pos_fw_id> kpfi;
     {
       const size_t n_tigs = tigs.size();
       vec<size_t> koffset(n_tigs + 1, 0);
       for (size_t i_tig = 0; i_tig != n_tigs; i_tig++) {
         const size_t n_kmers = (tigs[i_tig].size() >= unsigned(K)) ? 
           tigs[i_tig].size() - K + 1 : 0;
         koffset[i_tig + 1] = koffset[i_tig] + n_kmers;
       }
       kpfi.resize(koffset.back());

       // do we need this?
       vec<int> shuffle;
       Shuffle(n_tigs, shuffle);
       
       #pragma omp parallel for
       for (size_t ii = 0; ii < n_tigs; ii++) {
         const size_t i_tig = shuffle[ii];
         const size_t n_bases = tigs[i_tig].size();
         if (n_bases >= unsigned(K)) {
           size_t n_kmers = n_bases - K + 1;
           vec<char> akmer(K);
           FASTAVECTOR kmer;
           basevector b, brc;
           for (unsigned int ik = 0; ik < n_kmers; ik++) {
             b.SetToSubOf(tigs[i_tig], ik, K);
             brc = b;
             brc.ReverseComplement();
             basevector& ab = (b < brc ? b : brc);
             for (int l = 0; l < K; l++)
               akmer[l] = as_base(ab[l]);
             kmer.GetFromChars(akmer);
             kpfi[koffset[i_tig] + ik] = 
               kmer_pos_fw_id(kmer, ik, (b < brc ? True : False), i_tig);
           }
         }
       }
     }

     if (PRINT_STATUS) cout << Date( ) << ": starting parallel sort, mem usage = "
          << ToStringAddCommas( MemUsageBytes( ) ) << endl;
     sortInPlaceParallel( kpfi.begin( ), kpfi.end( ), NUM_THREADS );
     if (PRINT_STATUS) cout << Date( ) << ": sort complete, mem usage = "
          << ToStringAddCommas( MemUsageBytes( ) ) << endl;
     size_t N = kpfi.size( );

     // Split kpfi into kmers and pfi.

     if (PRINT_STATUS) cout << Date( ) << ": splitting kpfi" << endl;
     vec<char20> kmers(N); 
     vec<pos_fw_id> pfi(N);
     for ( size_t i = 0; i < N; i++ )
     {    kmers[i] = kpfi[i].kmer;
          pfi[i].pos = kpfi[i].pos;
          pfi[i].fw = kpfi[i].fw;
          pfi[i].id = kpfi[i].id;    }
     Destroy(kpfi);

     // Record lengths.

     if (PRINT_STATUS) cout << Date( ) << ": recording lengths" << endl;
     vec<int32_t> len(N);
     for ( size_t i = 0; i < N; i++ )
     {    size_t j;
          for ( j = i + 1; j < N; j++ )
               if ( kmers[j] != kmers[i] ) break;
          len[i] = j - i;
          i = j - 1;    }

     // Build hash table.  If an entry in the hash table is used, we use the
     // next entry, circling back to beginning if we reach the end of the table.  

     if (PRINT_STATUS) 
     {    cout << Date( ) << ": building hash table, N = " 
               << ToStringAddCommas(N) << endl;    }
     int64_t nhash = 5 * N;
     vec<int64_t> hash( nhash, -1 );
     int nproc = omp_get_max_threads( );
     vec<size_t> tstarts(nproc+1);
     for ( int i = 0; i < nproc; i++ )
     {    tstarts[i] = (i*N)/nproc;
          if ( i > 0 )
          {    while( tstarts[i] > tstarts[i-1]
                    && kmers[ tstarts[i] ] == kmers[ tstarts[i]-1 ] )
               {    tstarts[i]--;    }    }    }
     tstarts[nproc] = N;
     {    vec< vec< pair<int64_t,size_t> > > key_id(nproc);
          #pragma omp parallel for
          for ( int pi = 0; pi < nproc; pi++ )
          {    for ( size_t i = tstarts[pi]; i < tstarts[pi+1]; i++ )
               {    size_t j;
                    for ( j = i + 1; j < N; j++ )
                         if ( kmers[j] != kmers[i] ) break;
                    int64_t x = 0;
                    for ( int l = 0; l < 5; l++ )
                    {    x += int64_t( kmers[i].x[l] );
                         x = x << 8;    }
                    int64_t key = x % nhash;
                    key_id[pi].push( key, i );    
                    i = j - 1;    }
               Sort( key_id[pi] );    }
          for ( int pi = 0; pi < nproc; pi++ )
          {    for ( size_t l = 0; l < key_id[pi].size( ); l++ )
               {    int64_t key = key_id[pi][l].first;
                    while( hash[key] >= 0 ) 
                    {    key++;
                         if ( key == nhash ) key = 0;    }
                    hash[key] = key_id[pi][l].second;    }    }    }
     if (PRINT_STATUS)
     {    cout << Date( ) << ": memory in use after initial phase = "
               << ToStringAddCommas( MemUsageBytes( ) ) << endl;    }

     // Set up to track stuff.

     vec<read_loc> locs;

     // Go through three passes, for the three read types.

     for ( int rpass = 0; rpass < 3; rpass++ )
     {    if ( !lib_types[rpass] ) continue;
          String type_name 
               = ( rpass == 0 ? "frag" : ( rpass == 1 ? "jump" : "long_jump" ) );
          String lib_name 
               = (rpass == 0 ? FRAG_READS: (rpass == 1 ? JUMP_READS:LONG_JUMP_READS) );
          if ( !IsRegularFile( run_dir + "/" + lib_name + ".fastb" ) )
          {    continue;    }

          // Load reads.

          if (PRINT_STATUS) 
               cout << Date( ) << ": loading " << type_name << " data" << endl;
          vecbasevector bases( run_dir + "/" + lib_name + ".fastb" );
          PairsManager pairs( run_dir + "/" + lib_name + ".pairs" );

          // Go through the read pairs.

          size_t last_locs = locs.size( );
          ProcessPairs( rpass, pairs, bases, tigs, scaffolds, to_super, 
               to_super_pos, K, hash, kmers, pfi, len, PRINT_STATUS, PRINT_STATS,
               trace_ids, PRINT_BADS, PRINT_ALIGNMENTS, PRINT_ALIGNMENTS_TIG, 
               PRINT_MISSING, PRINT_INCONSISTENT, locs );

          // For long jumps, we sort and clean the read locations.

          if ( rpass == 2 )
          {    if (PRINT_STATUS) 
                    cout << Date( ) << ": sorting long jump read locs" << endl;
               sortInPlaceParallel( locs.begin( ) + last_locs, locs.end( ),
                                       NUM_THREADS );
               if (PRINT_STATUS)
                    cout << Date( ) << ": cleaning long jump read locs" << endl;
               vec<uint64_t> keep( bases.size( ), False );
               vec<size_t> locs_start( ntigs + 1 );
               {    size_t POS = last_locs;
                    for ( uint32_t u = 0; u <= ntigs; u++ )
                    {    while( POS < locs.size( ) && locs[POS].ContigId( ) < u ) 
                              ++POS;
                         locs_start[u] = POS;    }    }
               #pragma omp parallel for
               for ( uint it = 0; it < ntigs; it++ )
               {    size_t low = locs_start[it], high = locs_start[it+1];
                    for ( size_t i = low; i < high; i++ )
                    {    size_t j;
                         for ( j = i + 1; j < high; j++ )
                         {    if ( !SameStart( locs[j], locs[i] ) ) break;
                              if ( !SameFw( locs[j], locs[i] ) ) break;    }
                         const size_t max_pile = 1000;
                         if ( j - i > max_pile )
                         {    i = j - 1;
                              continue;    }
     
                         // Now all the reads between i and j are on the same contig
                         // and have the same start point and orientation.  Attempt
                         // to find the consensus of their partners, then pick the
                         // read that is closest to the consensus.
                         // 1. We are munging together all the long jump libraries 
                         //    here.  This isn't correct but compensates for 
                         //    incorrect breaking up of true libraries into 
                         //    sublibraries in input data.
                         // 2. This gets executed from both sides so two different
                         //    pairs may be saved.
                         // 3. This is totally brute force and could be done much
                         //    more subtly.

                         vec<basevector> partners;
                         int maxp = 0;
                         for ( size_t u = i; u < j; u++ )
                         {    const basevector& b = bases[ locs[u].PartnerReadId() ];
                              partners.push_back(b);
                              maxp = Max( maxp, b.isize( ) );    }
                         basevector consensus(maxp);
                         for ( int l = 0; l < maxp; l++ )
                         {    vec<int> count( 4, 0 ), id( 4, vec<int>::IDENTITY );
                              for ( int m = 0; m < partners.isize( ); m++ )
                              {     if ( partners[m].isize( ) > l )
                                        count[ partners[m][l] ]++;    }
                              ReverseSortSync( count, id );
                              consensus.Set( l, id[0] );    }
                         int max_agrees = -1, best = -1;
                         for ( int l = 0; l < partners.isize( ); l++ )
                         {    int agrees = 0;
                              for ( int m = 0; m < partners[l].isize( ); m++ )
                                   if ( partners[l][m] == consensus[m] ) agrees++;
                              if ( agrees > max_agrees )
                              {    max_agrees = agrees;
                                   best = l;    }    }
                         // I think we can get away without making the following
                         // two lines single-threaded, but am not 100% sure.
                         keep[ locs[i+best].ReadId( ) ] = True;
                         keep[ locs[i+best].PartnerReadId( ) ] = True;    
                         i = j - 1;    }    }
               size_t count = last_locs;
               for ( size_t x = last_locs; x < locs.size( ); x++ )
               {    if ( keep[ locs[x].ReadId( ) ] )
                    {    if ( count != x ) locs[count] = locs[x];
                         count++;    }    }
               locs.resize(count);

                // Check for duplicates that have nearby start points.

               {    size_t POS = last_locs;
                    for ( uint32_t u = 0; u <= ntigs; u++ )
                    {    while( POS < locs.size( ) && locs[POS].ContigId( ) < u ) 
                              ++POS;
                         locs_start[u] = POS;    }    }
               equiv_rel_64 e( bases.size( ) );
               const int max_start_diff = 5;
               #pragma omp parallel for
               for ( uint it = 0; it < ntigs; it++ )
               {    size_t low = locs_start[it], high = locs_start[it+1];
                    for ( size_t i = low; i < high; i++ )
                    {    size_t j, j_next = high;
                         for ( j = i + 1; j < high; j++ )
                         {    if ( j_next == high && !SameStart( locs[j], locs[i] ) )
                                   j_next = j;
                              if ( locs[j].Start( ) - locs[i].Start( ) 
                                   > max_start_diff )
                              {    break;    }    }
                         const size_t max_pile = 1000;
                         if ( j - i > max_pile )
                         {    i = j_next - 1;
                              continue;    }
                         for ( size_t u1 = i; u1 < j; u1++ )
                         {    const read_loc& rl1 = locs[u1];
                              for ( size_t u2 = u1 + 1; u2 < j; u2++ )
                              {    const read_loc& rl2 = locs[u2];
                                   if ( !SameFw( rl1, rl2 ) ) continue;
                                   int64_t pid1 = rl1.PseudoPairId( );
                                   int64_t pid2 = rl2.PseudoPairId( );
                                   if (!rl1.PartnerPlaced() && !rl2.PartnerPlaced())
                                        e.Join( pid1, pid2 ); 
                                   if ( !rl1.PartnerPlaced( ) ) continue;
                                   if ( !rl2.PartnerPlaced( ) ) continue;
                                   if ( !SamePartnerContigId( rl1, rl2 ) ) continue;
                                   if ( Abs( rl1.PartnerStart( )
                                        - rl2.PartnerStart( ) ) > max_start_diff ) 
                                   {    continue;    }
                                   if ( !SamePartnerFw( rl1, rl2 ) ) continue;
                                   e.Join( pid1, pid2 );    }    }
                         i = j_next - 1;    }    }
               vec<Bool> to_remove( locs.size( ) - last_locs, False );
               for ( size_t x = last_locs; x < locs.size( ); x++ )
               {    if ( !e.Representative( locs[x].PseudoPairId( ) ) )
                         to_remove[ x - last_locs ] = True;    }
               count = last_locs;
               for ( size_t x = last_locs; x < locs.size( ); x++ )
               {    if ( !to_remove[ x - last_locs ] )
                    {    if ( count != x ) locs[count] = locs[x];
                         count++;    }    }
               locs.resize(count);    }

          // Announce status.

          if (PRINT_STATUS)
          {    cout << Date( ) << ": locs so far = " 
                    << ToStringAddCommas( locs.size( ) ) << "; mem usage = "
                    << ToStringAddCommas( MemUsageBytes( ) ) << endl;    }    }

     // Sort, clean and write read locations.  The cleaning step removes some
     // duplicates, but it should remove more, and should select which to keep,
     // rather than just keep the first one.  Also there are issues with how
     // start points are defined.
     //
     // For long jumps, we remove cross-library duplicates.  This isn't really
     // right but compensates for incorrect library definitions that present one
     // true library as several libraries.

     if ( TRACE_IDS == "" )
     {    if (PRINT_STATUS) cout << Date( ) << ": sorting read locations" << endl;
          sortInPlaceParallel( locs.begin( ), locs.end( ), NUM_THREADS );
          if (PRINT_STATUS) cout << Date( ) << ": removing duplicates" << endl;
          vec<Bool> to_remove( locs.size( ), False );
          vec<size_t> locs_start( ntigs + 1 );
          {    size_t POS = 0;
               for ( uint32_t u = 0; u <= ntigs; u++ )
               {    while( POS < locs.size( ) && locs[POS].ContigId( ) < u ) ++POS;
                    locs_start[u] = POS;    }    }
          #pragma omp parallel for
          for ( uint it = 0; it < ntigs; it++ )
          {    size_t low = locs_start[it], high = locs_start[it+1];
               for ( size_t i = low; i < high; i++ )
               {    size_t j;
                    for ( j = i + 1; j < high; j++ )
                    {    if ( !SameStart( locs[j], locs[i] ) ) break;
                         if ( !SameFw( locs[j], locs[i] ) ) break;    }
                    const size_t max_pile = 1000;
                    if ( j - i > max_pile )
                    {    i = j - 1;
                         continue;    }

                    // Now all the reads between i and j are on the same contig
                    // and have the same start point and orientation.

                    for ( size_t u1 = i; u1 < j; u1++ )
                    {    const read_loc& rl1 = locs[u1];
                         if ( !rl1.PartnerPlaced( ) ) continue;
                         for ( size_t u2 = u1 + 1; u2 < j; u2++ )
                         {    const read_loc& rl2 = locs[u2];
                              if ( !rl2.PartnerPlaced( ) ) continue;
                              if ( rl1.ReadClass( ) != rl2.ReadClass( ) ) continue;
                              if ( !rl1.LongJump( ) && !SameLibraryId( rl1, rl2 ) )
                                   continue;
                              if ( !SamePartnerContigId( rl1, rl2 ) ) continue;
                              if ( !SamePartnerStart( rl1, rl2 ) ) continue;
                              if ( !SamePartnerFw( rl1, rl2 ) ) continue;
                              uint64_t low_id = Min( rl1.ReadId( ), 
                                   rl1.PartnerReadId( ), rl1.ReadId( ), 
                                   rl2.PartnerReadId( ) );
                              Bool keep1 = ( low_id == rl1.ReadId( ) 
                                   || low_id == rl1.PartnerReadId( ) );
                              if (keep1) to_remove[u2] = True;
                              else to_remove[u1] = True;    }    }

                    // For long jumps, if a partner of one of the reads is placed
                    // on a different contigs, kill all pairs having an unplaced
                    // partner.

                    Bool long_partner_on_diff = false;
                    for ( size_t u1 = i; u1 < j; u1++ )
                    {    const read_loc& rl1 = locs[u1];
                         if ( !rl1.LongJump( ) ) continue;
                         if ( !rl1.PartnerPlaced( ) ) continue;
                         if ( rl1.PartnerContigId( ) != rl1.ContigId( ) )
                         {    long_partner_on_diff = True;
                              break;    }    }
                    if (long_partner_on_diff)
                    {    for ( size_t u1 = i; u1 < j; u1++ )
                         {    const read_loc& rl1 = locs[u1];
                              if ( !rl1.LongJump( ) ) continue;
                              if ( rl1.PartnerPlaced( ) ) continue;
                              to_remove[u1] = True;    }    }

                    // Kill long jump pairs that are very short.

                    const int min_sep = 1000;
                    for ( size_t u = i; u < j; u++ )
                    {    const read_loc& rl = locs[u];
                         if ( !rl.LongJump( ) ) continue;
                         if ( !rl.PartnerPlaced( ) ) continue;
                         if ( rl.PartnerContigId( ) != rl.ContigId( ) ) continue;
                         if ( rl.Fw( ) == rl.PartnerFw( ) ) continue;
                         if ( Abs( rl.ActualSep( ) ) < min_sep )
                              to_remove[u] = True;    }

                    i = j - 1;    }    }

          EraseIf( locs, to_remove );
          if (PRINT_STATUS) cout << Date( ) << ": writing read locations" << endl;
          WriteReadLocs( head, tigs.size( ), locs );
          if (PRINT_STATUS) cout << Date( ) << ": done" << endl;    }

     // Display pileups.

     int64_t uncovered = 0, agree_lt_50 = 0, total = 0;

     // Done.

     if (PRINT_STATUS) 
     {    cout << "\n" << Date( ) << ": done, time used = " 
               << TimeSince(begin_clock) << endl;    }    }
