///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// AssemblyAccuracy.  Given an assembly consisting of scaffolds, and a reference
// sequence for the genome, assess the accuracy of the contigs that comprise
// the scaffolds as follows.  Take all contigs.  Break contigs >= 20 kb into equal 
// chunks of size ~10 kb.  Find the best alignment of each chunk to the reference:
// the one with the smallest number of errors (substitution or indel bases).
//
// Then divide chunks into six classes:
// 1. perfect - no errors
// 2. good    - error rate greater than 0, less than 0.1%
// 3. fair    - error rate at least 0.1%, less than 1.0%
// 4. bad     - error rate at least 1.0%, less than 10.0%
// 5. evil    - error rate at least 10.0%, but has 100-mer match
// 6. novel   - has no 100-mer match (indicative of DNA not present in the reference)
//
// Report the distribution of chunks by class.
//
// Chromosomes are allowed to be circular.
//
// Gaps in the assembly and reference should be denoted by a sequence of one
// or more ns.  (Note lower case.)
//
// Sequences of ns or Ns are treated as 'pseudogaps' and their size is optimized
// to improve the alignment score, increased or decreased arbitrarily.
//
// Heuristic note: now we only look at offsets whose maximum perfect stretch is 
// at least 1/4 of the best maximum perfect stretch.  This speeds up the code but
// could change the answer in rare cases.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <strstream>

#include <omp.h>

#include "Alignment.h"
#include "Basevector.h"
#include "Bitvector.h"
#include "Fastavector.h"
#include "FastIfstream.h"
#include "FetchReads.h"
#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "ParseSet.h"
#include "PrintAlignment.h"
#include "STLExtensions.h"
#include "Superb.h"
#include "VecUtilities.h"
#include "efasta/EfastaTools.h"
#include "math/Functions.h"
#include "math/HoInterval.h"
#include "pairwise_aligners/SmithWatFree.h"
#include "random/Random.h"
#include "random/Shuffle.h"
#include "reporting/PerfStat.h"
#include "system/SortInPlace.h"

// Horrors.  A char100 is really 100 bases.

class char100 {
     public:
     char x[25];
     friend Bool operator==( const char100& I1, const char100& I2 )
     {    return memcmp( I1.x, I2.x, 25 ) == 0;    }
     friend Bool operator!=( const char100& I1, const char100& I2 )
     {    return memcmp( I1.x, I2.x, 25 ) != 0;    }
     friend Bool operator<( const char100& I1, const char100& I2 )
     {    return memcmp( I1.x, I2.x, 25 ) < 0;    }
     void GetFromChars( const vec<char>& b )
     {    for ( int i = 0; i < 25; i++ )
               x[i] = 0;
          for ( int i = 0; i < 100; i++ )
          {    if ( b[i] == 'A' );
               else if ( b[i] == 'C' ) x[i/4] ^= ( 1 << 2*(i%4) );
               else if ( b[i] == 'G' ) x[i/4] ^= ( 2 << 2*(i%4) );
               else if ( b[i] == 'T' ) x[i/4] ^= ( 3 << 2*(i%4) );
               else
               {    cout << "GetFromChars: illegal character " << b[i] << endl;
                    cout << "Abort." << endl;
                    exit(1);    }    }    }
               
};

Bool TryAdd( const fastavector& c, fastavector& r, ho_interval& g, 
     const int add, int& shift, int& errs, int& best_loc, alignment& a,
     const Bool go_to_empty = False )
{    unsigned int sh = 0, start = g.Start( );
     while( ( g.Length( ) > -add || ( g.Length( ) >= -add && go_to_empty ) )
          && c.size( ) <= r.size( ) + add )
     {    fastavector rnew;
          if ( add > 0 )
          {    for ( unsigned int z = 0; z < r.size( ); z++ )
               {    if ( z == start )
                    {    for ( int j = 0; j < add; j++ )
                              rnew.push_back( 'n' );    }
                    rnew.push_back( r[z] );    }    }
          else
          {    for ( unsigned int z = 0; z < r.size( ); z++ )
               {    if ( z >= start && z < start - add ) continue;
                    rnew.push_back( r[z] );    }    }
          int best_loc_tmp;
          alignment a_tmp;
          int errs2 = SmithWatFree(c, rnew, best_loc_tmp, a_tmp, false, False, 1, 1);
          if ( errs2 >= errs ) break;
          r = rnew, errs = errs2, best_loc = best_loc_tmp, a = a_tmp;
          sh += add;
          g.AddToStop(add);    }
     shift += sh;
     return sh != 0;    }

Bool TrySub( const fastavector& c, fastavector& r, ho_interval& g, 
     const int sub, int& shift, int& errs, int& best_loc, alignment& a,
     const Bool go_to_empty = False )
{    return TryAdd( c, r, g, -sub, shift, errs, best_loc, a, go_to_empty );    }

typedef char100 FASTAVECTOR;

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
          if ( x1.pos < x2.pos ) return True;
          if ( x1.pos > x2.pos ) return False;
          if ( x1.fw < x2.fw ) return True;
          if ( x1.fw > x2.fw ) return False;
          if ( x1.id < x2.id ) return True;
          return False;    }

          // return x1.kmer < x2.kmer;    }

};

void LoadAssembly( const String& ASSEMBLY, const String& ASSEMBLY_TRUSTED,
     vec<fastavector>& assembly, const vecbitvector& assembly_trusted, 
     vec<bitvector>& atrusted, vec< triple<int,int,int> >& contig_origin )
{   
     fast_ifstream in(ASSEMBLY);
     int scaffold_count = 0;
     char c;
     fastavector fv;
     String line;
     while( !in.fail( ) ) 
     {    getline( in, line );
          ForceAssertEq( line[0], '>' );
          fv.clear( );
          while ( !in.fail( ) ) 
          {    getline( in, line );
               for ( int i = 0; i < line.isize( ); i++ )
                    fv.push_back( line[i] );
               in.peek(c);
               if ( c == '>' ) break;    }
          for ( size_t j = 0; j < fv.size( ); j++ )
               if ( fv[j] != 'n' ) fv[j] = toupper( fv[j] );
          vec<fastavector> fvs = fv.SplitOnGaps( );
          int contig_count = 0;
          for ( unsigned int j = 0; j < fv.size( ); j++ )
          {    while( fv[j] == 'n' ) 
               {    ++j;
                    if ( j == fv.size( ) ) break;    }
               if ( j == fv.size( ) ) break;
               unsigned int k;
               for ( k = j + 1; k < fv.size( ); k++ )
                    if ( fv[k] == 'n' ) break;
               contig_origin.push( scaffold_count, contig_count, j );
               assembly.push_back( fvs[contig_count] );
               if ( ASSEMBLY_TRUSTED != "" )
               {    bitvector T(k-j);
                    for ( unsigned int u = j; u < k; u++ )
                         T.Set( u-j, assembly_trusted[scaffold_count][u] );
                    atrusted.push_back(T);    }
               contig_count++;
               j = k;    }
          scaffold_count++;    }    }

void LoadAssemblyFromEfasta( const String& ASSEMBLY, vec<efasta>& contigs,
     vec< triple<int,int,int> >& contig_origin )
{    vec<efasta> scaffolds;
     LoadEfastaIntoStrings( ASSEMBLY, scaffolds );
     vec<superb> superbs;
     SplitEfastaIntoContigs( scaffolds, contigs, superbs );
     for ( size_t i = 0; i < superbs.size( ); i++ )
     {    const superb& s = superbs[i];
          int pos = 0;
          for ( int j = 0; j < s.Ntigs( ); j++ )
          {    contig_origin.push( i, j, pos );    
               pos += s.Len(j);
               if ( j < s.Ntigs( ) - 1 ) pos += s.Gap(j);    }    }    }

void GetMaxPerf( const fastavector& c0, const vec<fastavector>& ref,
     const vec< triple<int,int,Bool> >& offset, vec<int>& max_perf )
{    max_perf.resize( offset.size( ), 0 );
     for ( int j = 0; j < offset.isize( ); j++ )
     {    int o = offset[j].second;
          const fastavector& R = ref[ offset[j].first ];
          fastavector c(c0);
          if ( !offset[j].third ) c.ReverseComplement( );
          for ( int l = 0; l < (int) c.size( ); l++ )
          {    if ( o+l < 0 ) continue;
               if ( o+l >= (int) R.size( ) ) break;
               if ( !GeneralizedBase::fromChar(c[l])
                         .matches( GeneralizedBase::fromChar(R[o+l]) ) )
               {    continue;    }
               int l2 = l+1;
               while(1)
               {    if ( l2 >= (int) c.size( ) || o+l2 >= (int) R.size( ) ) 
                         break;
                    if ( !GeneralizedBase::fromChar(c[l2])
                         .matches( GeneralizedBase::fromChar(R[o+l2]) ) )
                    {    break;    }
                    l2++;    }
               max_perf[j] = Max( max_perf[j], l2 - l );
               l = l2 - 1;    }    }    }

int MatchingBases( const alignment& a, const fastavector& c, const fastavector& r )
{    int matching = 0;
     align A(a);
     int p1 = A.pos1( ), p2 = A.pos2( );
     for ( int z = 0; z < A.Nblocks( ); z++ )
     {    if ( A.Gaps(z) > 0 ) p2 += A.Gaps(z);
          if ( A.Gaps(z) < 0 ) p1 -= A.Gaps(z);
          for ( int x = 0; x < A.Lengths(z); x++ )
          {    if ( c[p1] == r[p2] ) ++matching;
               ++p1; ++p2;    }    }
     return matching;    }

void TwiddleGapSizes( const fastavector& c, fastavector& r,
     vec<ho_interval>& gaps, const vec<Bool>& faraway, 
     alignment& a, int& errs, int& best_loc )
{    for ( int x = 0; x < gaps.isize( ); x++ )
     {    ho_interval g = gaps[x];
          if (faraway[x]) continue;
          int shift = 0;
          Bool sub100 = False, sub10 = False, sub1 = False;
          Bool add100 = False, add10 = False, add1 = False;
          sub100 = TrySub( c, r, g, 100, shift, errs, best_loc, a );
          sub10 = TrySub( c, r, g, 10, shift, errs, best_loc, a );
          sub1 = TrySub( c, r, g, 1, shift, errs, best_loc, a );
          if ( !sub100 )
               add100 = TryAdd( c, r, g, 100, shift, errs, best_loc, a );
          if ( !sub10 )
               add10 = TryAdd( c, r, g, 10, shift, errs, best_loc, a );
          if ( !sub1 )
               add1 = TryAdd( c, r, g, 1, shift, errs, best_loc, a );
          if ( add100 && !add10 && !add1 )
               sub10 = TrySub( c, r, g, 10, shift, errs, best_loc, a );
          if ( add10 && !add1 )
               sub1 = TrySub( c, r, g, 1, shift, errs, best_loc, a );
          for ( int y = x + 1; y < gaps.isize( ); y++ )
               gaps[y].Shift(shift);    }    }

void ImproveGapSizes( alignment& a, const fastavector& c, fastavector& r,
     int& errs, int& best_loc )
{
     // Find number of matching bases.

     int matching = MatchingBases( a, c, r ), total = c.size( );
	       
     // Find gaps in r, then try to improve the gap sizes.  Don't do
     // this for truly terrible alignments, because it's expensive.

     if ( double(matching)/double(total) >= 0.75 ) 
     {    vec<ho_interval> gaps;
          for ( unsigned int x = 0; x < r.size( ); x++ )
          {    if ( toupper(r[x]) == 'N' )
               {    unsigned int y;
                    for ( y = x + 1; y < r.size( ); y++ )
    		         if ( toupper(r[y]) != 'N' ) break;
                    gaps.push( x, y );
                    x = y;    }    }
	       
          // Determine which gaps are proximate to indels.
          // We only mess with these gaps.

          vec<Bool> faraway( gaps.size( ), False );
          const int far = 100;
          for ( int x = 0; x < gaps.isize( ); x++ ) 
          {    align A(a);
               int p2 = A.pos2( );
               for ( int z = 0; z < A.Nblocks( ); z++ ) 
               {    if ( A.Gaps(z) > 0 ) p2 += A.Gaps(z);
                    if ( p2 <= gaps[x].Start( ) - far 
 	                 && p2 + A.Lengths(z) >= gaps[x].Stop( ) + far ) 
                    {    faraway[x] = True;
                         break;    }
                    p2 += A.Lengths(z);    }    }
	          
          TwiddleGapSizes( c, r, gaps, faraway, a, errs, best_loc );
          for ( int x = 0; x < gaps.isize( ); x++ )
          {    ho_interval g = gaps[x];
               if (faraway[x]) continue;
               int shift = 0;
               if ( g.Length( ) == 1 )
               {    TrySub( c, r, g, 1, shift, errs, best_loc, 
                         a, True );    }
               for ( int y = x + 1; y < gaps.isize( ); y++ )
                    gaps[y].Shift(shift);    }     }     }

void CorrectErrorCountUsingTrusted( const alignment& a_best, 
     const bitvector& T_best, const fastavector& c_best,
     const fastavector& r_best, int& min_errs )
{    align A(a_best);
     min_errs = 0;
     unsigned int p1 = A.pos1( ), p2 = A.pos2( );
     for ( int j = 0; j < A.Nblocks( ); j++ )
     {    if ( A.Gaps(j) > 0 )
          {    Bool trusted = T_best[p1]
                    && ( p1+1 == T_best.size( ) || T_best[p1+1] );
               if ( !trusted ) min_errs += A.Gaps(j);
               p2 += A.Gaps(j);    }
          if ( A.Gaps(j) < 0 )
          {    Bool trusted = True;
               for ( unsigned int k = p1; k < p1 - A.Gaps(j); k++ )
                    if ( !T_best[k] ) trusted = False;
               if ( !trusted ) min_errs -= A.Gaps(j);
               p1 -= A.Gaps(j);    }
          for ( int x = 0; x < A.Lengths(j); x++ )
          {    if ( !GeneralizedBase::fromChar(c_best[p1])
                    .matches( GeneralizedBase::fromChar(r_best[p2]) ) )
               {    if ( !T_best[p1] ) min_errs++;    }
               ++p1; ++p2;    }    }    }

void GetSuperAmbiguities( const efasta& e, vec< pair<int,int> >& samb, 
     const int flank )
{
     // First find the coordinates of the ambiguities, then the
     // super-ambiguities.  Extend up to end if distance to end is less than flank.

     vec< pair<int,int> > amb;
     for ( int j = 0; j < e.isize( ); j++ )
     {    if ( e[j] != '{' ) continue;
          int k;
          for ( k = j + 1; k < e.isize( ); k++ )
               if ( e[k] == '}' ) break;
          amb.push( j, k+1 );
          j = k;    }
     samb.clear( );
     for ( int j = 0; j < amb.isize( ); j++ )
     {    int k;
          for ( k = j + 1; k < amb.isize( ); k++ )
               if ( amb[k].first - amb[k-1].second >= flank ) break;
          if ( amb[j].first < flank ) amb[j].first = 0;
          if ( e.isize( ) - amb[k-1].second < flank ) amb[k-1].second = e.size( );
          samb.push( amb[j].first, amb[k-1].second );
          j = k - 1;    }    }

class efasta_match {

     public:

     efasta_match( ) { }
     efasta_match( const int matches, const int mismatches, const int path_id,
          const int e1, const int e2 ) : matches(matches), mismatches(mismatches), 
          path_id(path_id), e1(e1), e2(e2) { }

     friend Bool operator<( const efasta_match& m1, const efasta_match& m2 )
     {    if ( m1.mismatches < m2.mismatches ) return True;
          if ( m1.mismatches > m2.mismatches ) return False;
          if ( m1.e1 < m2.e1 ) return True;
          if ( m1.e1 > m2.e1 ) return False;
          if ( m1.e2 < m2.e2 ) return True;
          if ( m1.e2 > m2.e2 ) return False;
          return m1.matches > m2.matches;    }

     int matches;
     int mismatches;
     int path_id;
     int e1, e2;

};


int main(int argc, char **argv)
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_Doc(ASSEMBLY,
          "fasta file for assembly; ns treated as gaps");
     CommandArgument_String_Doc(REF, "fasta file for reference" );
     CommandArgument_UnsignedInt_OrDefault(CHUNK_SIZE, 3000);
     CommandArgument_UnsignedInt_OrDefault(MIN_CONTIG, 1000);
     CommandArgument_Int_OrDefault(VERBOSITY, 0);
     CommandArgument_Bool_OrDefault_Doc(IMPERFECT_ONLY, False,
          "if VERBOSITY=2, only print imperfectly aligning chunks");
     CommandArgument_Bool_OrDefault_Doc(LOG_PROGRESS, False,
          "list chunks as they are aligned");
     CommandArgument_Int_OrDefault_Doc(SAMPLE, 10000,
          "if nonzero, pick this many chunks randomly instead of using them all");
     CommandArgument_UnsignedInt_OrDefault_Doc(MIN_TRUSTED_STRETCH, 15, 
          "if a sequence of "
          "less than this many is bracketed by Ns or ns, change to Ns");
     CommandArgument_String_OrDefault_Doc(QUERIES_TO_USE, "", "use this list of "
          "assembly records if specified");
     CommandArgument_UnsignedInt_OrDefault_Doc(MIN_CHUNK, 0,
          "ignore chunks smaller than this");
     CommandArgument_String_OrDefault(ASSEMBLY_TRUSTED, "");
     CommandArgument_Bool_OrDefault(EFASTA, False);
     CommandArgument_Bool_OrDefault(EFASTA_VERBOSE, False);
     CommandArgument_Int_OrDefault(MAX_INSTANTIATIONS, 4);
     CommandArgument_Int_OrDefault(MAX_COUNT, 250);
     CommandArgument_Int_OrDefault(CHUNK_ID, -1);
     CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
	   "Number of threads to use (use all available processors if set to 0)");

     EndCommandArguments;

     // Check arguments, etc.

     double clock = WallClockTime( );
     if ( SAMPLE > 0 && QUERIES_TO_USE != "" )
     {    cout << "SAMPLE and QUERIES_TO_USE options cannot be used "
               << "simultaneously." << endl;
          exit(1);    }

     // Thread control

     NUM_THREADS = configNumThreads(NUM_THREADS);
     omp_set_num_threads( NUM_THREADS );

     // Load the reference and build version that accommodates circularity.

     vec<fastavector> ref;
     LoadFromFastaFile( REF, ref );
     for ( size_t i = 0; i < ref.size( ); i++ )
     {    for ( size_t j = 0; j < ref[i].size( ); j++ )
               if ( ref[i][j] != 'n' ) ref[i][j] = toupper( ref[i][j] );    }
     for ( int i = 0; i < ref.isize( ); i++ )
     {    int add = Min( ref[i].size( ), 2*CHUNK_SIZE + CHUNK_SIZE/5 );
          fastavector fv2;
          fv2.SetToSubOf( ref[i], 0, add );
          ref[i] = Cat( ref[i], fv2 );    }
     size_t nref = ref.size( );
     vec<String> refnames;
     fast_ifstream rin(REF);
     String line;
     while(1)
     {    getline( rin, line );
          if ( rin.fail( ) ) break;
          if ( !line.Contains( ">", 0 ) ) continue;
          if ( line.Contains( " " ) ) line = line.Before( " " );
          refnames.push_back( line.After( ">" ) );    }

     // Load trusted bases on assembly.

     vecbitvector assembly_trusted;
     if ( ASSEMBLY_TRUSTED != "" ) assembly_trusted.ReadAll(ASSEMBLY_TRUSTED);

     // Clear untrusted stretches.

     for ( size_t t = 0; t < nref; t++ )
     {    for (unsigned int i = 0; i < ref[t].size( ); i++ )
          {    if ( toupper(ref[t][i]) == 'N' ) continue;
               unsigned int j;
               for ( j = i + 1; j < ref[t].size( ); j++ )
                    if ( toupper(ref[t][j]) == 'N' ) break;
               if ( j - i < MIN_TRUSTED_STRETCH )
               {    for ( unsigned int k = i; k < j; k++ )
                         ref[t][k] = 'N';    }
               i = j;    }    }

     // Load assembly and break it into contigs.  Note that code is in parallel
     // for the fasta and efasta cases.

     vec<fastavector> assembly;
     vec<efasta> contigs;
     vec<bitvector> atrusted;
     vec< triple<int,int,int> > contig_origin;
     if ( ASSEMBLY.Contains( ".efasta", -1 ) ) EFASTA = True;
     if ( !EFASTA )
     {    LoadAssembly( ASSEMBLY, ASSEMBLY_TRUSTED, assembly, assembly_trusted, 
               atrusted, contig_origin );    }
     else
     {    ForceAssertEq( ASSEMBLY_TRUSTED, "" );
          LoadAssemblyFromEfasta( ASSEMBLY, contigs, contig_origin );    }

     // Divide assembly into chunks.

     cout << Date() << ": breaking assembly into chunks" << endl;
     vec<fastavector> chunks;
     vec<efasta> chunkse;
     vec<bitvector> tchunks;
     vec< pair<int,int> > chunk_origin;
     vec<int> contig_sizes;
     if ( !EFASTA )
     {    for ( int i = 0; i < assembly.isize( ); i++ )
          {    const fastavector& A = assembly[i];
               contig_sizes.push_back( A.size( ) );
               if ( A.size( ) < MIN_CONTIG ) continue;
               if ( A.size( ) < CHUNK_SIZE )
               {    chunks.push_back(A);
                    if ( ASSEMBLY_TRUSTED != "" ) tchunks.push_back( atrusted[i] );
                    chunk_origin.push( i, 0 );    }
               else
               {    unsigned int n = A.size( ) / CHUNK_SIZE;
                    unsigned int ch = 1 + A.size( ) / n;
                    for ( unsigned int j = 0; j < A.size( ); j += ch )
                    {    fastavector b;
                         b.SetToSubOf( A, j, Min( ch, A.size( ) - j ) );
                         chunks.push_back(b);
                         if ( ASSEMBLY_TRUSTED != "" )
                         {    bitvector t;
                              t.SetToSubOf( 
                                   atrusted[i], j, Min( ch, A.size( ) - j ) );
                              tchunks.push_back(t);    }
                         chunk_origin.push( i, j/ch );    }    }    }    }
     else
     {    for ( int i = 0; i < contigs.isize( ); i++ )
          {    const efasta& A = contigs[i];
               contig_sizes.push_back( A.Length1( ) );
               if ( A.Length1( ) < (int) MIN_CONTIG ) continue;
               if ( A.Length1( ) < (int) CHUNK_SIZE )
               {    chunkse.push_back(A);
                    chunk_origin.push( i, 0 );    }
               else
               {    int LEN = A.Length1( );
                    int n = LEN/CHUNK_SIZE;
                    int ch = 1 + LEN/n, chunk_count = 0;
                    for ( int j = 0; j < LEN; j++ )
                    {    efasta b;
                         b.reserve(ch);
                         int len = 0;
                         while( j < A.isize( ) && len < ch )
                         {    if ( A[j] == '{' ) 
                              {    b.push_back( A[j++] );
                                   while( j < A.isize( ) && A[j] != ',' )
                                   {    b.push_back( A[j++] );
                                        len++;    }
                                   while( j < A.isize( ) && A[j] != '}' )
                                        b.push_back( A[j++] );
                                   b.push_back( A[j++] );    }
                              else
                              {    b.push_back( A[j++] );
                                   len++;     }    }
                         j--;
                         chunkse.push_back(b);
                         chunk_origin.push( i, chunk_count++ );    }    }    }    }
     int nchunks = ( !EFASTA ? chunks.size( ) : chunkse.size( ) );
     cout << Date( ) << ": " << contig_sizes.size( ) << " contigs, total size " 
          << ToStringAddCommas( BigSum(contig_sizes) ) << ", N50 = " 
          << N50(contig_sizes) << endl;
     cout << Date() << ": found " << nchunks << " chunks" << endl;

     // Select chunks.

     vec<Bool> use( nchunks, True );
     if ( SAMPLE == 0 ) SAMPLE = nchunks;
     int USING = nchunks;
     for ( int i = 0; i < nchunks; i++ )
     {    if ( !EFASTA && chunks[i].size( ) < MIN_CHUNK )
          {    use[i] = False;
               USING--;    }
          if ( EFASTA && chunkse[i].Length1( ) < (int) MIN_CHUNK )
          {    use[i] = False;
               USING--;    }    }
     while( SAMPLE < USING )
     {    int x = randomx( ) % nchunks;
          if ( use[x] )
          {    use[x] = False;
               USING--;    }    }
     if ( QUERIES_TO_USE != "" )
     {    vec<int> quse;
          ParseIntSet( QUERIES_TO_USE, quse );
          for ( int i = 0; i < nchunks; i++ )
          {    int m = chunk_origin[i].first;
               if ( !BinMember( quse, contig_origin[m].first ) )
                    use[i] = False;    }    }
     vec<Bool> nuse(nchunks);
     for ( int i = 0; i < nchunks; i++ )
          nuse[i] = !use[i];
     if ( !EFASTA ) EraseIf( chunks, nuse );
     else EraseIf( chunkse, nuse );
     if ( ASSEMBLY_TRUSTED != "" ) EraseIf( tchunks, nuse );
     EraseIf( chunk_origin, nuse );
     if (EFASTA)
     {    chunks.resize( chunkse.size( ) );
          for ( size_t i = 0; i < chunkse.size( ); i++ )
          {    fastavector v;
               chunkse[i].FlattenTo(v);
               chunks[i] = v;    }    }
     nchunks = chunks.size( );
     cout << Date() << ": using " << nchunks << " chunks" << endl;
     if ( VERBOSITY >= 1 ) cout << endl;

     // For EFASTA, for each ambiguity, attempt to determine which choice is right,
     // then modify chunkse and chunks accordingly.  Note that this computation
     // is somewhat redundant with the computation that follows.
     // 
     // The way this works is that for each ambiguity, we find the 100-mers that
     // flank it.  Ambiguities closer than 100 bases are lumped together as single
     // super-ambiguities.  Then we search the reference sequence for the flanking 
     // 100-mers.  The sequences between them are then used to determine which 
     // choices are kept.  There are some heuristic choices.

     typedef vec<kmer_pos_fw_id>::iterator Itr;
     typedef Comparator<kmer_pos_fw_id> Comp;
     InPlaceParallelSorter<Itr,Comp> ipps(NUM_THREADS);
     if (EFASTA)
     {    const int flank = 100;
          vec< triple< int, int, pair<basevector,basevector> > > flanks;
          for ( size_t i = 0; i < chunkse.size( ); i++ )
          {    const efasta& e = chunkse[i];

               // Define the super-ambiguities and then the flanks.

               vec< pair<int,int> > samb;
               GetSuperAmbiguities( e, samb, flank );
               for ( int j = 0; j < samb.isize( ); j++ )
               {    Bool at_left_end = ( samb[j].first == 0 );
                    Bool at_right_end = ( samb[j].second == e.isize( ) );
                    basevector left(flank), right(flank);
                    if (at_left_end) left.resize(0);
                    if (at_right_end) right.resize(0);
                    for ( int l = 0; l < flank; l++ )
                    {    if ( !at_left_end )
                              left.Set( l, as_char( e[samb[j].first - flank + l] ) );
                         if ( !at_right_end )
                              right.Set( l, as_char( e[samb[j].second + l] ) );    }
                    flanks.push( i, j, make_pair( left, right ) );    }    }
          size_t nref = ref.size( );
          vec<basevector> all;
          for ( size_t i = 0; i < nref; i++ )
               all.push_back( ref[i].ToBasevector( ) );
          vec<basevector> chunklets;
          for ( size_t i = 0; i < flanks.size( ); i++ )
               chunklets.push_back( flanks[i].third.first, flanks[i].third.second );
          size_t nchunklets = chunklets.size( );
          all.append(chunklets);
          vec<kmer_pos_fw_id> kpfi;
          vec<size_t> kcount( all.size( ), 0 ), koffset( all.size( ) + 1, 0 );
          vec<int> shuffle;
          Shuffle( all.size( ), shuffle );
          int K = 100;
          for ( int pass = 1; pass <= 2; pass++ )
          {    if ( pass == 2 )
               {    for ( int i = 1; i <= all.isize( ); i++ )
                         koffset[i] = koffset[i-1] + kcount[i-1];
                    kpfi.resize( koffset.back( ) );
                    for ( int i = 0; i < all.isize( ); i++ )
                         kcount[i] = 0;    }
               #pragma omp parallel for
               for ( size_t ip = 0; ip < all.size( ); ip++ )
               {    int i = shuffle[ip];
                    if ( all[i].size( ) == 0 ) continue;
                    vec<char> akmer(100);
                    FASTAVECTOR kmer;
                    for ( unsigned int j = 0; j <= all[i].size( ) - K; j++ )
                    {    if ( (size_t) i < nref )
                         {    Bool includes_N = False;
                              for ( int r = 0; r < K; r++ )
                                   if ( ref[i][j+r] == 'N' ) includes_N = True;
                              if (includes_N) continue;    }
                         basevector b, brc;
                         if ( pass == 2 )
                         {    b.SetToSubOf( all[i], j, K );
                              brc = b;
                              brc.ReverseComplement( );
                              basevector& ab = ( b < brc ? b : brc );
                              #pragma omp critical
                              for ( int l = 0; l < K; l++ )
                                   akmer[l] = as_base( ab[l] );
                              kmer.GetFromChars(akmer);
                              kpfi[ koffset[i] + kcount[i] ] = kmer_pos_fw_id(
                                   kmer, j, b < brc ? True : False, i );    }
                         kcount[i]++;    }    }    }
          cout << Date( ) << ": starting parallel sort 1, mem usage = "
               << ToStringAddCommas( MemUsageBytes( ) ) << endl;
          ipps.sort( kpfi.begin( ), kpfi.end( ) );
          cout << Date( ) << ": sort complete, mem usage = "
               << ToStringAddCommas( MemUsageBytes( ) ) << endl;
          vec< vec< triple<int,int,Bool> > > offset(nchunklets);
          for ( size_t i = 0; i < kpfi.size( ); i++ )
          {    size_t j;
               for ( j = i + 1; j < kpfi.size( ); j++ )
                    if ( kpfi[j].kmer != kpfi[i].kmer ) break;
               for ( size_t u1 = i; u1 < j; u1++ ) // index of reference contig
               {    if ( kpfi[u1].id >= nref ) continue;
                    for ( size_t u2 = i; u2 < j; u2++ ) // index of chunklet
                    {    if ( kpfi[u2].id < nref ) continue;
                         const basevector& chunklet 
                              = chunklets[ kpfi[u2].id - nref ];
                         int o;
                         if ( ( kpfi[u1].fw && kpfi[u2].fw ) 
                              || ( !kpfi[u1].fw && !kpfi[u2].fw ) )
                         {    o = kpfi[u1].pos - kpfi[u2].pos;    }
                         else 
                         {    o = kpfi[u1].pos 
                                   - (chunklet.size( ) - kpfi[u2].pos - K);    }
                         offset[ kpfi[u2].id - nref ].push( kpfi[u1].id, o, 
                              !( kpfi[u1].fw ^ kpfi[u2].fw ) );    }    }
               i = j - 1;    }
          vec< vec<int> > founds( flanks.size( ) );
          vec<efasta> chunkse_orig(chunkse);
          cout << Date( ) << ": searching flanks" << endl;
          for ( size_t i = 0; i < flanks.size( ); i++ )
          {    int fi = flanks[i].first, fj = flanks[i].second;
               if ( CHUNK_ID >= 0 && CHUNK_ID != fi ) continue;
               const efasta& e = chunkse_orig[fi];
               vec< pair<int,int> > samb;
               GetSuperAmbiguities( e, samb, flank );
               efasta sa;
               for ( int j = samb[fj].first; j < samb[fj].second; j++ )
                    sa.push_back( e[j] );
               if (EFASTA_VERBOSE)
               {    cout << "\n=================================================="
                         << "==================================\n\n";
                    cout << "chunk " << fi << ", ambiguity " << sa 
                         << ", left flank placed "
                         << offset[2*i].size( ) << " times, right flank placed "
                         << offset[2*i+1].size( ) << " times\n";    }
               vec< vec<basevector> > G;
               sa.MakeGraph(G);
               int maxlen = 0;
               for ( int j = 0; j < G.isize( ); j++ )
               {    int m = 0;
                    for ( int k = 0; k < G[j].isize( ); k++ )
                         m = Max( m, G[j][k].isize( ) );
                    maxlen += m;    }

               // Limit the number of locations for the flanks.

               const int max_locs = 100;
               if ( offset[2*i].isize( ) * offset[2*i+1].isize( ) > max_locs ) 
                    continue;

               // Limit the number of paths through the super-ambiguity.

               const int max_count = 1000;
               int count = 1;
               for ( int j = 0; j < G.isize( ); j++ )
               {    count *= G[j].size( );
                    if ( count > max_count ) break;    }
               if ( count > max_count ) continue;

               // Find all paths through the super-ambiguity.

               vec<basevector> paths;
               sa.ExpandTo(paths);

               // Go through all locations for the flanks and determine which paths
               // are compatible with them.  The empty flank cases (at contig end)
               // requires special handling.

               vec<efasta_match> MATCHES;
               Bool at_left_end = ( samb[fj].first == 0 );
               Bool at_right_end = ( samb[fj].second == e.isize( ) );
               if ( !at_left_end && !at_right_end )
               {    for ( int e1 = 0; e1 < offset[2*i].isize( ); e1++ )
                    for ( int e2 = 0; e2 < offset[2*i+1].isize( ); e2++ )
                    {    if (EFASTA_VERBOSE)
                         {    cout << "\nexamining placement pair (" << e1
                                   << "," << e2 << ")\n";    }
                         int t1 = offset[2*i][e1].first; 
                         int t2 = offset[2*i+1][e2].first;
                         int p1 = offset[2*i][e1].second; 
                         int p2 = offset[2*i+1][e2].second;
                         Bool fw1 = offset[2*i][e1].third; 
                         Bool fw2 = offset[2*i+1][e2].third;
                         if ( t1 != t2 || fw1 != fw2 ) continue;
                         int len = ( fw1 ? p2 - ( p1 + K ) : p1 - ( p2 + K ) );
                         if ( len < 0 ) continue;
                         if ( len > maxlen + 100 ) continue;
                         if (EFASTA_VERBOSE) 
                              PRINT6( t1, p1, int(fw1), t2, p2, int(fw2) );
                         int start = ( fw1 ? p1 + K : p2 + K );
                         basevector m;
                         m.SetToSubOf( all[t1], start, len );
                         if ( !fw1 ) m.ReverseComplement( );
                         if (EFASTA_VERBOSE) PRINT( m.ToString( ) );
                         for ( int l = 0; l < paths.isize( ); l++ )
                         {    if ( m == paths[l] ) 
                              {    if (EFASTA_VERBOSE) 
                                   {    cout << "MATCHES " << l << ", matches = " 
                                             << len << ", errs = 0" << endl;    }
                                   MATCHES.push( len, 0, l, 
                                        e1, e2 );    }    }    }    }
               const int smith_wat_extend = 10;
               if ( at_left_end && !at_right_end )
               {    for ( int e2 = 0; e2 < offset[2*i+1].isize( ); e2++ )
                    {    if (EFASTA_VERBOSE)
                         {    cout << "\n";
                              PRINT(e2);    }
                         int t2 = offset[2*i+1][e2].first;
                         int p2 = offset[2*i+1][e2].second;
                         Bool fw2 = offset[2*i+1][e2].third;
                         for ( int l = 0; l < paths.isize( ); l++ )
                         {    int matches = 0, mismatches = 0;
                              basevector p = paths[l];
                              if ( !fw2 ) p.ReverseComplement( );
                              int np = p.size( ), N2 = all[t2].size( );
                              int best_loc;
                              alignment al;
                              basevector T2;
                              int start, stop;
                              if (fw2)
                              {
                                   // Does p match all[t2], ending at p2?

                                   start = Max( 0, p2 - np - smith_wat_extend );
                                   stop = p2;    }
                              else
                              {
                                   // Does p match all[t2], starting at p2+flank?

                                   start = p2+flank;
                                   stop = start + np + smith_wat_extend;    }
                              if ( p.size( ) > 0 && p.isize( ) <= stop - start )
                              {    T2.SetToSubOf( all[t2], start, stop - start );
                                   SmithWatFree( p, T2, best_loc, al, !fw2, fw2 );
                                   align a(al);
                                   mismatches = a.Errors( p, T2 );
                                   matches = a.MatchingBases( p, T2 );    }
                              else mismatches = matches = 0;
                              if (EFASTA_VERBOSE) 
                              {    cout << "MATCHES " << l << ", matches = "
                                        << matches << ", errs = "
                                        << mismatches << endl;    }
                              MATCHES.push( matches, mismatches, 
                                   l, -1, e2 );    }    }    }
               if ( !at_left_end && at_right_end )
               {    for ( int e1 = 0; e1 < offset[2*i].isize( ); e1++ )
                    {    if (EFASTA_VERBOSE)
                         {    cout << "\n";
                              PRINT(e1);    }
                         int t1 = offset[2*i][e1].first; 
                         int p1 = offset[2*i][e1].second; 
                         Bool fw1 = offset[2*i][e1].third; 
                         for ( int l = 0; l < paths.isize( ); l++ )
                         {    int matches = 0, mismatches = 0;
                              basevector p = paths[l];
                              if ( !fw1 ) p.ReverseComplement( );
                              int np = p.size( ), N1 = all[t1].size( );
                              int best_loc;
                              alignment al;
                              basevector T1;
                              int start, stop;
                              if ( !fw1 )
                              {
                                   // Does p match all[t1], ending at p1?

                                   start = Max( 0, p1 - np - smith_wat_extend );
                                   stop = p1;    }
                              else
                              {
                                   // Does p match all[t1], starting at p1+flank?

                                   start = p1+flank;
                                   stop = start + np + smith_wat_extend;    }
                              T1.SetToSubOf( all[t1], start, stop - start );
                              if ( p.size( ) > 0 && p.isize( ) <= stop - start )
                              {    SmithWatFree( p, T1, best_loc, al, fw1, !fw1 );
                                   align a(al);
                                   mismatches = a.Errors( p, T1 );
                                   matches = a.MatchingBases( p, T1 );    }
                              else mismatches = matches = 0;
                              if (EFASTA_VERBOSE) 
                              {    cout << "MATCHES " << l << ", matches = "
                                        << matches << ", errs = "
                                        << mismatches << endl;    }
                              MATCHES.push( matches, mismatches, 
                                   l, e1, 0 );    }    }    }
               if ( at_left_end && at_right_end ) continue; // must be rare
               Sort(MATCHES);
               int el;
               for ( el = 0; el < MATCHES.isize( ); el++ )
                    if ( MATCHES[el].mismatches > MATCHES[0].mismatches ) break;
               MATCHES.resize(el);
               vec<int> equals;

               for ( int u = 0; u < MATCHES.isize( ); u++ )
               {    int v;
                    for ( v = u + 1; v < el; v++ )
                    {    if ( MATCHES[v].e1 != MATCHES[u].e1 ) break;
                         if ( MATCHES[v].e2 != MATCHES[u].e2 ) break;    }
                    for ( int w = u; w < v; w++ )
                    {    if ( MATCHES[w].matches < MATCHES[u].matches ) break;
                         equals.push_back( MATCHES[w].path_id );    }
                    u = v - 1;    }

               /*
               int el;
               for ( el = 0; el < MATCHES.isize( ); el++ )
               {    if ( MATCHES[el].matches != MATCHES[0].matches
                         || MATCHES[el].mismatches != MATCHES[0].mismatches )
                    {    break;    }    }
               Bool too_close = False;
               for ( int u = el; u < MATCHES.isize( ); u++ )
               {    if ( !( MATCHES[u].matches <= MATCHES[0].matches ) ) 
                         too_close = True;
                    if ( !( MATCHES[u].mismatches <= MATCHES[0].mismatches ) )
                         too_close = True;    }
               if (too_close) continue;
               for ( int j = 0; j < el; j++ )
                    equals.push_back( MATCHES[j].path_id );
               */

               UniqueSort(equals);
               if (EFASTA_VERBOSE)
               {    cout << "\nFOUND:";
                    for ( int j = 0; j < equals.isize( ); j++ )
                         cout << " " << equals[j];
                    cout << endl;    }
               if ( equals.empty( ) ) continue;

               // Expand sa into the list v of all sequences it represents and the
               // list vx of the origins of those sequences.

               vec<basevector> v;
               vec< vec<int> > vx;
               v.push_back( basevector( ) ), vx.push_back( vec<int>( ) );
               for ( int l = 0; l < G.isize( ); l++ )
               {    vec<basevector> w;
                    vec< vec<int> > wx;
                    for ( size_t r = 0; r < v.size( ); r++ )
                    {    for ( int j = 0; j < G[l].isize( ); j++ )
                         {    w.push_back( Cat( v[r], G[l][j] ) );    
                              vec<int> z( vx[r] );
                              z.push_back(j); 
                              wx.push_back(z);    }    }
                    v = w, vx = wx;    }

               // Now find the origins foundsx of the sequences that were found.

               vec< vec<int> > foundsx;
               for ( int l = 0; l < equals.isize( ); l++ )
                    foundsx.push_back( vx[ equals[l] ] );

               // Transpose.

               vec< vec<int> > foundsy( G.size( ) );
               for ( int l = 0; l < foundsx.isize( ); l++ )
               {    for ( int r = 0; r < foundsx[l].isize( ); r++ )
                         foundsy[r].push_back( foundsx[l][r] );    }
               for ( int l = 0; l < foundsy.isize( ); l++ )
                    UniqueSort( foundsy[l] );

               // Find the sequence sanew that should replace sa.

               efasta sanew;
               for ( int l = 0; l < G.isize( ); l++ )
               {    vec<Bool> to_delete( G[l].size( ), True );
                    for ( int m = 0; m < foundsy[l].isize( ); m++ )
                         to_delete[ foundsy[l][m] ] = False;
                    EraseIf( G[l], to_delete );    }
               sanew.MakeFromGraph(G);

               // Substitute sanew into chunkse, padding with blanks.

               for ( int l = 0; l < sanew.isize( ); l++ )
                    chunkse[fi][ l + samb[fj].first ] = sanew[l];
               for (int l = sanew.isize( ); l < samb[fj].second-samb[fj].first; l++)
                    chunkse[fi][ l + samb[fj].first ] = ' ';    }

          // Clean up blanks.

          for ( int r = 0; r < chunkse.isize( ); r++ )
          {    chunkse[r].GlobalReplaceBy( " ", "" );
               if ( VERBOSITY >= 4 )
                    chunkse[r].Print( cout, "chunk_" + ToString(r) );
               chunkse[r].FlattenTo( chunks[r] );    }    }

     // Find all 100-mer matches between the chunks and the reference.  From this
     // compute all possible offsets: the inferred start positions of the chunk
     // or its reverse complement on the reference.

     cout << Date() << ": finding matches, chunks vs reference, mem = " 
          << ToStringAddCommas( MemUsageBytes( ) ) << endl;
     int K = 100;
     vec<fastavector> all(ref);
     all.append(chunks);
     vec<kmer_pos_fw_id> kpfi;
     vec<size_t> kcount( all.size( ), 0 ), koffset( all.size( ) + 1, 0 );
     vec<int> shuffle;
     Shuffle( all.size( ), shuffle );
     for ( int pass = 1; pass <= 2; pass++ )
     {    if ( pass == 2 )
          {    for ( int i = 1; i <= all.isize( ); i++ )
                    koffset[i] = koffset[i-1] + kcount[i-1];
               kpfi.resize( koffset.back( ) );
               for ( int i = 0; i < all.isize( ); i++ )
                    kcount[i] = 0;    }
          #pragma omp parallel for
          for ( size_t ip = 0; ip < all.size( ); ip++ )
          {    int i = shuffle[ip];
               if ( all[i].size( ) == 0 ) continue;
               vec<char> akmer(100);
               FASTAVECTOR kmer;
               fastavector f;
               String choices;
               vec<fastavector> X(MAX_INSTANTIATIONS), X2(MAX_INSTANTIATIONS);
               unsigned int nextamb;
               for ( nextamb = 0; nextamb < all[i].size( ); nextamb++ )
               {    char c = toupper( all[i][nextamb] );
                    if ( c != 'A' && c != 'C' && c != 'G' && c != 'T' ) break;    }
               for ( unsigned int j = 0; j <= all[i].size( ) - K; j++ )
               {    Bool thisamb = ( j == nextamb );
                    if (thisamb)
                    {    for ( nextamb = j+1; nextamb < all[i].size( ); nextamb++ )
                         {    char c = toupper( all[i][nextamb] );
                              if ( c != 'A' && c != 'C' && c != 'G' && c != 'T' ) 
                                   break;    }    }
                    fastavector b, brc;
                    if ( !thisamb && nextamb - j > static_cast<unsigned>(K) )
                    {    if ( pass == 2 )
                         {    b.SetToSubOf( all[i], j, K );
                              brc = b;
                              brc.ReverseComplement( );
                              fastavector& ab = ( b < brc ? b : brc );
                              for ( int l = 0; l < K; l++ )
                                   akmer[l] = ab[l];
                              kmer.GetFromChars(akmer);
                              kpfi[ koffset[i] + kcount[i] ]
                                   = kmer_pos_fw_id(
                                   kmer, j, b < brc ? True : False, i );    }
                         kcount[i]++;
                         continue;    }
                    b.SetToSubOf( all[i], j, K );
                    int n = 0;
                    for ( int l = 0; l < K; l++ )
                    {    if ( toupper(b[l]) == 'N' ) ++n;
                         else break;    }
                    if ( n == K ) continue;
                    brc = b;
                    brc.ReverseComplement( );
                    fastavector& ab = ( b < brc ? b : brc );
                    int instantiations = 1;
                    for ( int l = 0; l < K; l++ )
                    {    const char& c = ab[l];
                         if ( c == 'A' || c == 'C' || c == 'G' || c == 'T' ) 
                              continue;
                         if ( c == 'N' ) instantiations *= 4;
                         else
                         {    choices = GeneralizedBase::fromChar(ab[l]).bases( );
                              instantiations *= choices.isize( );    }
                         if ( instantiations > MAX_INSTANTIATIONS ) break;    }
                    if ( instantiations > MAX_INSTANTIATIONS ) continue;
                    if ( instantiations == 1 )
                    {    if ( pass == 2 )
                         {    for ( int l = 0; l < K; l++ )
                                   akmer[l] = ab[l];
                              kmer.GetFromChars(akmer);
                              kpfi[ koffset[i] + kcount[i] ]
                                   = kmer_pos_fw_id( 
                                   kmer, j, b < brc ? True : False, i );    }
                         kcount[i]++;    }
                    else
                    {    int Xcount = 1, X2count = 0;
                         X[0].clear( );
                         for ( int l = 0; l < K; l++ )
                         {    choices = GeneralizedBase::fromChar(ab[l]).bases( );
                              for ( int m = 0; m < Xcount; m++ )
                              {    for ( int n = 0; n < choices.isize( ); n++ )
                                   {    f = X[m];
                                        f.push_back( choices[n] );
                                        X2[X2count++] = f;    }    }
                              for ( int m = 0; m < X2count; m++ )
                                   X[m] = X2[m];
                              Xcount = X2count, X2count = 0;    }
                         for ( int m = 0; m < Xcount; m++ )
                         {    if ( pass == 2 )
                              {    for ( int l = 0; l < K; l++ )
                                        akmer[l] = X[m][l];
                                   kmer.GetFromChars(akmer);
                                   kpfi[ koffset[i] + kcount[i] ]
                                        = kmer_pos_fw_id(
                                        kmer, j, b < brc ? True : False, i );    }
                              kcount[i]++;    }    }    }    }    }
     cout << Date( ) << ": starting parallel sort, mem usage = "
          << ToStringAddCommas( MemUsageBytes( ) ) << endl;
     ipps.sort( kpfi.begin( ), kpfi.end( ) );
     cout << Date( ) << ": sort complete, mem usage = "
          << ToStringAddCommas( MemUsageBytes( ) ) << endl;
     vec< vec< triple<int,int,Bool> > > offset(nchunks);
     for ( size_t i = 0; i < kpfi.size( ); i++ )
     {    size_t j;
          for ( j = i + 1; j < kpfi.size( ); j++ )
               if ( kpfi[j].kmer != kpfi[i].kmer ) break;
          for ( size_t u1 = i; u1 < j; u1++ ) // index of reference contig
          {    if ( kpfi[u1].id >= nref ) continue;
               for ( size_t u2 = i; u2 < j; u2++ ) // index of chunk
               {    if ( kpfi[u2].id < nref ) continue;
                    const fastavector& chunk = chunks[ kpfi[u2].id - nref ];
                    int o;
                    if ( ( kpfi[u1].fw && kpfi[u2].fw ) 
                         || ( !kpfi[u1].fw && !kpfi[u2].fw ) )
                    {    o = kpfi[u1].pos - kpfi[u2].pos;    }
                    else o = kpfi[u1].pos - ( chunk.size( ) - kpfi[u2].pos - K );
                    offset[ kpfi[u2].id - nref ].push(
                         kpfi[u1].id, o, !( kpfi[u1].fw ^ kpfi[u2].fw ) );    }    }
          i = j - 1;    }

     // Compute minimum number of errors in each chunk.

     cout << Date( ) << ": finding errors in assembly chunks" << endl;
     longlong total = 0;
     longlong perfect = 0, good = 0, fair = 0, bad = 0, evil = 0, novel = 0;
     longlong base_err_I_II = 0, base_total_I_II = 0;
     longlong base_err_I_II_III = 0, base_total_I_II_III = 0;
     longlong base_err_all = 0, base_total_all = 0;
     vec<String> chunk_reports(nchunks);
     int parallel_chunk_size = Max( 1, nchunks / ( 10 * omp_get_max_threads( ) ) );
     int chunks_to_go = nchunks;
     #pragma omp parallel for schedule(dynamic, parallel_chunk_size)
     for ( int i = 0; i < nchunks; i++ )
     {    if ( CHUNK_ID >= 0 && i != CHUNK_ID ) continue;
          double iclock = WallClockTime( );
          int m = chunk_origin[i].first;
          if (LOG_PROGRESS)
          {
               #pragma omp critical
               {    cout << "start aligning chunk " << i+1 << " = chunk " 
                         << chunk_origin[i].second << " from scaff.tig " 
                         << contig_origin[m].first << "." 
                         << contig_origin[m].second << endl;    }    }

          UniqueSort( offset[i] );

          // Expand ambiguities in chunks, within reason.

          vec<fastavector> these_chunks;
          if ( !EFASTA ) these_chunks.push_back( chunks[i] );
          else
          {    Bool OK = chunkse[i].ExpandTo( these_chunks, MAX_COUNT );
               if ( !OK ) these_chunks.push_back( chunks[i] );    }

          // Set up to track best alignment.  There are several defects in how we
          // do this in the EFASTA case, including:
          // - We only compute n once.
          // - The offsets may be a bit off.

          const vec< triple<int,int,Bool> >& offseti = offset[i];
          int n;
          if ( !EFASTA ) n = chunks[i].size( );
          else n = chunkse[i].Length1( );
          int min_errs = n;
          fastavector r, r_best, c, c_best;
          bitvector T, T_best;
          Bool rc = False;
          alignment a, a_best;
          int t_best = -1, start_best = -1;

          // For each offset, determine the maximum perfect stretch.  We only do
          // this for the first chunk, which isn't quite right, however in order
          // to do it for all the chunks we would have to fix the offsets.

          vec<int> max_perf;
          GetMaxPerf( these_chunks[0], ref, offseti, max_perf );
          int MAX_PERF = 0;
          if ( max_perf.nonempty( ) ) MAX_PERF = Max(max_perf);
          if ( VERBOSITY >= 3 )
          {
               #pragma omp critical
               PRINT2( i, MAX_PERF );    }

          // Go through the chunks.

          for ( size_t tc = 0; tc < these_chunks.size( ); tc++ )
          {
               // If we already know there's a perfect alignment, we're done.
               // However, this only works for the case tc = 0 because of how 
               // MAX_PERF is calculated.
     
               if ( tc == 0 && MAX_PERF == (int) these_chunks[tc].size( ) ) 
               {    int j;
                    for ( j = 0; j < offseti.isize( ); j++ )
                         if ( max_perf[j] == MAX_PERF ) break;
                    min_errs = 0;
                    t_best = offseti[j].first;
                    c_best = these_chunks[tc];
                    if ( !offseti[j].third ) c_best.ReverseComplement( );
                    start_best = offseti[j].second;
                    const fastavector& R = ref[ offseti[j].first ];
                    r_best.SetToSubOf( R, start_best, c_best.size( ) );
                    avector<int> gaps(1), lengths(1);
                    gaps(0) = 0;
                    lengths(0) = c_best.size( );
                    align al( 0, 0, gaps, lengths );
                    a_best.Set( al, 0 );
                    break;    }
     
               // Otherwise, we have to Smith-Waterman.

               for ( int j = 0; j < offseti.isize( ); j++ )
               {    int j2;
                    for ( j2 = j+1; j2 < offseti.isize( ); j2++ )
                    {    if ( offseti[j2].first != offseti[j].first ) break;
                         if ( offseti[j2].third != offseti[j].third ) break;
                         if ( offseti[j2].second - offseti[j2-1].second > 500 ) 
                              break;    }

                    // Only consider offsets whose maximum perfect stretch is
                    // at least 1/4 of the best possible.

                    int mp = 0;
                    for ( int z = j; z < j2; z++ )
                         mp = Max( mp, max_perf[z] );
                    if ( mp < MAX_PERF/4 )
                    {    j = j2 - 1;
                         continue;    }

                    // Now align.

                    const fastavector& R = ref[ offseti[j].first ];
                    int start = Max( 0, offseti[j].second - n/10 );
                    int odiff = offseti[j2-1].second - offseti[j].second;
                    r.SetToSubOf( R, start, Min( n + n/5 + odiff, 
                         static_cast<int>(R.size()) - start ) );
                    c = these_chunks[tc];
                    if ( !offseti[j].third ) c.ReverseComplement( );
                    if ( ASSEMBLY_TRUSTED != "" )
                    {    T = tchunks[i];
                         if ( !offseti[j].third ) T.ReverseMe( );    }
                    int best_loc;
                    if ( c.size( ) > r.size( ) ) continue;
                    int errs = SmithWatFree( 
                         c, r, best_loc, a, false, False, 1, 1 );

                    // Try to improve gap sizes.

                    ImproveGapSizes( a, c, r, errs, best_loc );

                    // Record alignment if best so far.
     
                    if ( errs < min_errs )
                    {    min_errs = errs;
                         a_best = a, c_best = c, T_best = T, r_best = r;    
                         rc = !offseti[j].third;    
                         t_best = offseti[j].first;    
                         start_best = start;    }
                    j = j2 - 1;    }
     
               // Correct error count, taking into account trusted bases.

               if ( ASSEMBLY_TRUSTED != "" )
               {    CorrectErrorCountUsingTrusted( a_best, T_best, c_best,
                         r_best, min_errs );    }    }

          // Report results.

          double err_rate = double(min_errs) / double(n);
          if ( VERBOSITY >= 1 )
          {    ostringstream out;
               out << setprecision(3) << 100.0 * err_rate << "% error rate";
               out << ", " << min_errs << " errors of " << n;
               out << ", chunk " << i << " = chunk " << chunk_origin[i].second
                    << " from scaff.tig " << contig_origin[m].first << "."
                    << contig_origin[m].second;
               if ( VERBOSITY >= 2 && t_best >= 0 )
               {    out << "\n\naligns " << ( rc ? "rc" : "fw" ) << " to reference "
                         << "contig " << t_best << " (" << refnames[t_best] 
                         << ") starting at " << start_best + a_best.pos2( ) << "\n";
                    PrintVisualAlignment( True, out, c_best, r_best, a_best );    }
               chunk_reports[i] = out.str( );    }
          #pragma omp critical
          {    total += n;
               if ( min_errs == 0 ) perfect += n;
               else if ( err_rate < 0.001 ) good += n;
               else if ( err_rate < 0.01 ) fair += n;
               else if ( err_rate < 0.1 ) bad += n;
               else if ( min_errs < n ) evil += n;
               else novel += n;
               if ( err_rate < 0.001 )
               {    base_total_I_II += n;
                    base_err_I_II += min_errs;    }
               if ( err_rate < 0.01 )
               {    base_total_I_II_III += n;
                    base_err_I_II_III += min_errs;    }
               base_total_all += n;
               base_err_all += min_errs;    
               if (LOG_PROGRESS)
               {    cout << "done - chunk " << i+1 << " = chunk " 
                         << chunk_origin[i].second << " of contig " 
                         << contig_origin[m].first << "." 
                         << contig_origin[m].second << ", used "
                         << TimeSince(iclock) << ", " << --chunks_to_go
                         << " chunks left" << endl;    }    }    }
     int report_counter = 1;
     if ( VERBOSITY >= 1 )
     {    for ( int i = 0; i < nchunks; i++ )
          {    if ( CHUNK_ID >= 0 && CHUNK_ID != i ) continue;
               if ( IMPERFECT_ONLY && chunk_reports[i].Contains( " 0 errors" ) )
                    continue;
               cout << "====================================================="
                    << "===========================\n\n";
               chunk_reports[i].GlobalReplaceBy( "\n\n\n", "\n\n" );
               if ( chunk_reports[i].Contains( "\n\n", -1 ) )
                    chunk_reports[i].resize( chunk_reports[i].isize( ) - 1 );
               cout << "[" << report_counter++ << "] "
                    << chunk_reports[i] << endl;    }    }

     // Generate summary statistics.

     cout << Date( ) << ": generating summary" << endl;
     double perf_chunks = 100.0 * double(perfect) / double(total);
     cout << "\nperfect   " << PERCENT_RATIO( 3, perfect, total )
          << "   " << perfect << endl;
     cout << "good      " << PERCENT_RATIO( 3, good, total )
          << "   " << good << endl;
     cout << "fair      " << PERCENT_RATIO( 3, fair, total )
          << "   " << fair << endl;
     cout << "bad       " << PERCENT_RATIO( 3, bad, total )
          << "   " << bad << endl;
     cout << "evil      " << PERCENT_RATIO( 3, evil, total )
          << "   " << evil << endl;
     cout << "novel     " << PERCENT_RATIO( 3, novel, total )
          << "   " << novel << endl;
     cout << "\nbase accuracy from chunks I, II: Q";
     if ( base_err_I_II > 0 )
     {    cout << -10.0 * log10(
               double(base_err_I_II)/double(base_total_I_II) ) << "\n";     }
     else cout << "infinity\n";
     double base_acc = 
          -10.0 * log10( double(base_err_I_II_III)/double(base_total_I_II_III) );
     cout << "base accuracy from chunks I, II, III: Q";
     if ( base_err_I_II_III > 0 ) cout << base_acc << "\n";
     else cout << "infinity\n";
     cout << "total error rate (including novel) = "
          << 10000.0 * double(base_err_all) / double(base_total_all)
          << " x 10^-4" << endl;
     cout << "\n";
     if ( perfect == total )
     {    PerfStat::log( ) 
               << PerfStat( "perf_chunks", "% of contig chunks that are perfect",
               100 );    }
     else
     {    PerfStat::log( ) << std::fixed << std::setprecision(1)
               << PerfStat( "perf_chunks", "% of contig chunks that are perfect",
               perf_chunks );    }
     double bad_evil_chunks = 100.0 * double(bad+evil) / double(total);
     PerfStat::log( ) << std::fixed << std::setprecision(1)
          << PerfStat( "bad_evil_chunks", "% of contig chunks with 1%+ errors",
          bad_evil_chunks );
     PerfStat::log( ) << std::fixed << std::setprecision(1)
          << PerfStat( "base_accuracy", "estimated base accuracy (Q)", base_acc );
     cout << "\ntotal time used = " << TimeSince(clock) << endl;
     cout << Date( ) << ": done" << endl;    }
