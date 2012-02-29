///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

const char *DOC =
  "Cleans up an assembly by removing small contigs and scaffolds, removing "
  "duplicate scaffolds, reordering scaffolds and renumber contigs, etc...";


// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "CoreTools.h"
#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "VecUtilities.h"
#include "math/Functions.h"
#include "paths/AssemblyCleanupTools.h"

/**
 * CleanAssembly
 *
 * Modify an assembly (remove small contigs and scaffolds, reorder according to length and
 * position in scaffolds).
 *
 * HEAD_IN: it loads <HEAD_IN>.{efasta,superb}
 * HEAD_OUT: it saves <HEAD_OUT>.{efasta,superb,fastb}
 * MIN_CONTIG_SIZE_SOLO: minimum solo contig size to keep
 * MIN_CONTIG_SIZE_IN: minimum contig size within scaffolds to keep
 * MIN_SCAFFOLD_SIZE: minimum scaffold size (ungapped) to keep
 * REORDER: order scaffolds numbering according to size and contig numbering according to 
 *        appearance in scaffolds
 * DEDUP: remove duplicate scaffolds (size and content)
 */

class char40 {
     public:
     char x[10];
     friend Bool operator==( const char40& I1, const char40& I2 )
     {    return memcmp( I1.x, I2.x, 10 ) == 0;    }
     friend Bool operator!=( const char40& I1, const char40& I2 )
     {    return memcmp( I1.x, I2.x, 10 ) != 0;    }
     friend Bool operator<( const char40& I1, const char40& I2 )
     {    return memcmp( I1.x, I2.x, 10 ) < 0;    }
     void GetFromString( const String& b )
     {    for ( int i = 0; i < 10; i++ )
               x[i] = 0;
          for ( int i = 0; i < 40; i++ )
          {    if ( b[i] == 'A' );
               else if ( b[i] == 'C' ) x[i/4] ^= ( 1 << 2*(i%4) );
               else if ( b[i] == 'G' ) x[i/4] ^= ( 2 << 2*(i%4) );
               else if ( b[i] == 'T' ) x[i/4] ^= ( 3 << 2*(i%4) );
               else
               {    cout << "GetFromString: illegal character " << b[i] << endl;
                    cout << "Abort." << endl;
                    exit(1);    }    }    }
               
};

int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandDoc(DOC);
  CommandArgument_String_Doc( HEAD_IN,
    "Assembly files to clean <HEAD_IN>.{superb,contigs.{efasta|fasta}}" );
  CommandArgument_String_Doc( HEAD_OUT,
    "Cleaned assembly files");
  CommandArgument_Int_OrDefault_Doc( MIN_CONTIG_SIZE_SOLO, 1000, 
    "Remove solo contigs that are smaller than MIN_CONTIG_SIZE_SOLO");
  CommandArgument_Int_OrDefault_Doc( MIN_CONTIG_SIZE_IN, 200, 
    "Remove contigs within scaffolds that are smaller than MIN_CONTIG_SIZE_IN");
  CommandArgument_Int_OrDefault_Doc( MIN_SCAFFOLD_SIZE, 1000, 
    "Remove scaffolds with ungapped size smaller than MIN_SCAFFOLD_SIZE");
  CommandArgument_Int_OrDefault_Doc( MIN_INDIRECT_GAP, 0, 
    "If specified, and the gap between contigs n and n+2 is less than the value, "
    "delete the intervening contig.");
  CommandArgument_Bool_OrDefault_Doc( REORDER, True,
    "Sort scaffolds by size (largest first) and renumber consecutive contigs");
  CommandArgument_Bool_OrDefault_Doc( DEDUP, True,
    "Remove Duplicate Singleton Scaffolds");
  CommandArgument_Int_OrDefault_Doc( MIN_UNIQUE, 0,
    "Delete nearly identical scaffolds that differ by less than MIN_UNIQUE kmers" );
  CommandArgument_Int_OrDefault_Doc( MAX_SWALLOW, 0,
    "Delete contigs smaller than MAX_SWALLOW that are adjacent to predicted "
    " overlaps that are at least as big. " );
  CommandArgument_UnsignedInt_OrDefault_Doc( NUM_THREADS, 0, 
    "Number of threads to use (use all available processors if set to 0)");
  EndCommandArguments;

  // Thread control
  
  NUM_THREADS = configNumThreads(NUM_THREADS);
  omp_set_num_threads( NUM_THREADS );
  
  // File names.
  String in_efasta_file = HEAD_IN + ".contigs.efasta";
  String in_fasta_file  = HEAD_IN + ".contigs.fasta";
  String in_superb_file = HEAD_IN + ".superb";

  String out_efasta_file        = HEAD_OUT + ".contigs.efasta";
  String out_fasta_file         = HEAD_OUT + ".contigs.fasta";
  String out_fastb_file         = HEAD_OUT + ".contigs.fastb";
  String out_superb_file        = HEAD_OUT + ".superb";

  String out_contig_map_file    = HEAD_OUT + ".contigs.mapping";
  String out_scaffold_map_file  = HEAD_OUT + ".superb.mapping";

  // Loading scaffolds
  cout << Date( ) << ": loading superb file" << endl;
  vec<superb> scaffoldsIn;
  ReadSuperbs( in_superb_file, scaffoldsIn );

  // Loading contigs
  vec<efasta> efastasIn;
  if ( IsRegularFile( in_efasta_file ) ){
    cout << Date() << ": Loading efasta" << endl;
    LoadEfastaIntoStrings(in_efasta_file, efastasIn);
  }else{
    cout << Date() << ": Loading fasta" << endl;
    if ( ! IsRegularFile( in_fasta_file ) )
      FatalErr(" neither " + in_efasta_file + " nor " + in_fasta_file + " was found!" );
    vec<fastavector> fastas;
    LoadFromFastaFile( in_fasta_file, fastas );
    efastasIn.resize( fastas.size() );
    for ( size_t i = 0; i < fastas.size(); i++ )
      efastasIn[i] = efasta( fastas[i] );
  }
  Assembly assembly( scaffoldsIn, efastasIn );

  size_t origScaffoldsTotLen = assembly.scaffoldsTotLen();
  PRINT( origScaffoldsTotLen );
  size_t origScaffoldsRedLen = assembly.scaffoldsRedLen();
  PRINT( origScaffoldsRedLen );
  size_t origScaffoldsNtigs = assembly.scaffoldsNtigs();
  PRINT( origScaffoldsNtigs );

  assembly.remove_unused_contigs();
  assembly.check_integrity();

  assembly.remove_small_contigs( MIN_CONTIG_SIZE_SOLO, MIN_CONTIG_SIZE_IN );
  assembly.remove_small_scaffolds( MIN_SCAFFOLD_SIZE );

  if (DEDUP) assembly.dedup( );

  if (REORDER) assembly.reorder( );

     // Delete tiny contigs that are adjacent to predicted overlaps that 
     // are at least as big.  Also remove tiny contigs in small indirect gaps.

     vec<Bool> to_remove( assembly.efastas.size( ), False );
     if ( MAX_SWALLOW > 0 )
     {    cout << Date( ) << ": deleting tiny contigs that are adjacent to "
               << "predicted overlaps that are at least as big" << endl;
          for ( int i = 0; i < assembly.scaffolds.isize( ); i++ )
          {    const superb& S = assembly.scaffolds[i];
               for ( int j = 0; j < S.Ntigs( ); j++ )
               {    if ( S.Len(j) > MAX_SWALLOW ) continue;
                    int min_gap = 1000000000;
                    if ( j > 0 ) min_gap = Min( min_gap, S.Gap(j-1) );
                    if ( j < S.Ngaps( ) - 1 ) min_gap = Min( min_gap, S.Gap(j) );
                    if ( S.Len(j) < -min_gap )
                         to_remove[ S.Tig(j) ] = True;    }    }    }
     if ( MIN_INDIRECT_GAP > 0 )
     {    for ( int i = 0; i < assembly.scaffolds.isize( ); i++ )
          {    const superb& S = assembly.scaffolds[i];
               for ( int j = 0; j < S.Ntigs( ) - 2; j++ )
               {    int indirect_gap 
                         = Max( 0, S.Gap(j) ) + Max( 0, S.Gap(j+1) ) + S.Len(j+1);
                    if ( indirect_gap < MIN_INDIRECT_GAP && !to_remove[j]
                         && !to_remove[j+2] )
                    {    to_remove[j+1] = True;    }    }    }    }
     assembly.remove_contigs(to_remove);

     // Delete scaffolds that do not have enough unique sequence in them.  In effect
     // we find pairs of scaffolds (s1, s2) for which the number of kmers in s1
     // that are not in s2 is less than a threshold.  Then provided that s2 is 
     // sufficiently bigger than s1, we delete s1.  We don't do this if s2 has 
     // already been deleted.

     if ( MIN_UNIQUE > 0 )
     {    cout << Date( ) << ": deleting scaffolds that lack enough unique "
               << "sequence" << endl;

          // Define heuristics.

          const int K = 40; // can't be changed without changing class char40
          const double min_size_mult = 2.0;
          const Bool unique_verbose = True;

          // Set up.

          int min_unique = MIN_UNIQUE - K + 1;
          vec<char40> kmers;
          vec<int> scaff_id;
          int N = assembly.efastas.size( ), ns = assembly.scaffolds.size( );
          vec<int> to_scaffold( N, -1 );
          for ( int s = 0; s < ns; s++ )
          {    const superb& S = assembly.scaffolds[s];
               for ( int j = 0; j < S.Ntigs( ); j++ )
                    to_scaffold[ S.Tig(j) ] = s;    }
          vec<int> len(N), tid( N, vec<int>::IDENTITY );
          for ( int i = 0; i < N; i++ )
               len[i] = assembly.efastas[i].size( );
          ReverseSortSync( len, tid );

          // Go through the contigs.

          #pragma omp parallel for schedule( dynamic, 1 )
          for ( int ti = 0; ti < N; ti++ )
          {    int t = tid[ti];
               vec<char40> tkmers;
               vec<int> tscaff_id;
               int s = to_scaffold[t];
               char40 b, brc;
               String x, xrc;
               const efasta& e = assembly.efastas[t];
               for ( int l = 0; l <= (int) e.size( ) - K; l++ )
               {    x = e.substr( l, K );
                    Bool legal = True;
                    for ( int i = 0; i < K; i++ )
                    {    if ( x[i] != 'A' && x[i] != 'C' && x[i] != 'G'
                              && x[i] != 'T' )
                         {    legal = False;    }    }
                    if ( !legal ) continue;
                    b.GetFromString(x);
                    StringReverseComplement( x, xrc );
                    brc.GetFromString(xrc);
                    tkmers.push_back( b, brc );
                    tscaff_id.push_back( s, s );    }
               #pragma omp critical
               {    kmers.append(tkmers);
                    scaff_id.append(tscaff_id);    }    }

          // Find candidates for deletion.

          ParallelSortSync( kmers, scaff_id );
          vec<int> total( ns, 0 ), unique( ns, 0 );
          for ( int64_t i1 = 0; i1 < (int64_t) kmers.size( ); i1++ )
          {    int s1 = scaff_id[i1];
               total[s1]++;
               Bool shared = False;
               for ( int64_t i2 = i1 - 1; i2 >= 0; i2-- )
               {    if ( kmers[i2] != kmers[i1] ) break;
                    if ( scaff_id[i2] != s1 )
                    {    shared = True;
                         break;    }    }
               for ( int64_t i2 = i1 + 1; i2 < (int64_t) kmers.size( ); i2++ )
               {    if ( kmers[i2] != kmers[i1] ) break;
                    if ( scaff_id[i2] != s1 )
                    {    shared = True;
                         break;    }    }
               if ( !shared ) unique[s1]++;    }
          vec<int> candidates;
          for ( int s = 0; s < ns; s++ )
          {    unique[s] /= 2, total[s] /= 2;
               if ( unique[s] < min_unique ) candidates.push_back(s);    }

          // Find scaffolds to be deleted.

          vec< vec<int> > owners( candidates.size( ) );
          for ( int64_t i1 = 0; i1 < (int64_t) kmers.size( ); i1++ )
          {    int s1 = scaff_id[i1];
               int p1 = BinPosition( candidates, s1 );
               if ( p1 < 0 ) continue;
               vec<int> o;
               for ( int64_t i2 = i1 - 1; i2 >= 0; i2-- )
               {    if ( kmers[i2] != kmers[i1] ) break;
                    int s2 = scaff_id[i2];
                    if ( total[s2] < min_size_mult * double( total[s1] ) ) continue;
                    if ( s2 != s1 ) o.push_back(s2);    }
               for ( int64_t i2 = i1 + 1; i2 < (int64_t) kmers.size( ); i2++ )
               {    if ( kmers[i2] != kmers[i1] ) break;
                    int s2 = scaff_id[i2];
                    if ( total[s2] < min_size_mult * double( total[s1] ) ) continue;
                    if ( s2 != s1 ) o.push_back(s2);    }
               UniqueSort(o);
               owners[p1].append(o);    }
          vec<Bool> to_delete( ns, False );
          for ( int c = 0; c < candidates.isize( ); c++ )
          {    int s1 = candidates[c];
               if ( to_delete[s1] ) continue;
               Sort( owners[c] );
               for ( int i = 0; i < owners[c].isize( ); i++ )
               {    int j = owners[c].NextDiff(i);
                    int n = ( j - i ) / 2;
                    int s2 = owners[c][i];
                    if ( !to_delete[s2] && total[s1] - n < min_unique )
                    {    to_delete[s1] = True;
                         if (unique_verbose)
                         {    cout << s1 << "[l=" << total[s1] << "] << " << s2 
                                   << "[l=" << total[s2] << "] (" << total[s1] - n 
                                   << " remaining)" << endl;    }
                         break;    }
                    i = j - 1;    }    }
          EraseIf( assembly.scaffolds, to_delete );
          EraseIf( assembly.scaffMap, to_delete );
          assembly.remove_unused_contigs( );    }

  assembly.check_integrity();

  size_t finalScaffoldsTotLen = assembly.scaffoldsTotLen();
  PRINT( finalScaffoldsTotLen );
  size_t finalScaffoldsRedLen = assembly.scaffoldsRedLen();
  PRINT( finalScaffoldsRedLen );

  // writing output
  cout << Date() << ": writing output files" << endl;
  assembly.WriteAll( HEAD_OUT );

  cout << Date( ) << ": input =  " << scaffoldsIn.size( ) << " scaffolds, "
       << efastasIn.size( ) << " contigs" << endl;
  cout << Date( ) << ": output = " << assembly.scaffolds.size( ) << " scaffolds, "
       << assembly.scaffoldsNtigs( ) << " contigs" << endl;
  cout << Date( ) << ": CleanAssembly done!\n" << endl;
}

