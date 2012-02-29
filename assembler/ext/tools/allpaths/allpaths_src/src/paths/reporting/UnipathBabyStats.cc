///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// UnipathBabyStats: compute a small number of unipath statistics appropriate for
// a regional assembly.  Note that genome_extended.fasta must exist.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

// MakeDepend: dependency HyperToReftigs
// MakeDepend: dependency AssemblyAccuracy

#include "Basevector.h"
#include "MainTools.h"
#include "math/Functions.h"
#include "paths/UnibaseUtils.h"
#include "paths/UnipathFixerTools.h"

int main( int argc, char** argv ) 
{
  
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String( PRE );
     CommandArgument_String( DATA );
     CommandArgument_String( RUN );
     CommandArgument_Int( K );
     CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
       "Number of threads to use (use all available processors if set to 0)");
     CommandArgument_String_OrDefault( READS, "all_reads" );
     CommandArgument_String_OrDefault( UNIBASES, "unibases" );
     EndCommandArguments;

     // Thread control (Uses OMP in UnipathFixerTools)
   
     NUM_THREADS = configNumThreads(NUM_THREADS);
     omp_set_num_threads( NUM_THREADS );

     // Define directories.
  
     String data_dir = PRE + "/" + DATA;
     String run_dir = data_dir + "/" + RUN;

     // Sanity check.

     ForceAssert( IsRegularFile( data_dir + "/../genome_extended.fasta" ) );

     // Set up output stream.

     Ofstream( sout, run_dir + "/UnipathBabyStats.out" );

     // Load unibases and compute involution.

     String unibases_head = run_dir + "/" + READS + "." + UNIBASES;
     String kK = ".k" + ToString(K);
     vecbasevector unibases( unibases_head + kK );
     vec<int> to_rc;
     UnibaseInvolution( unibases, to_rc, K );

     // Generate kmer summary stats.

     vecbasevector genome( data_dir + "/genome.fastb" );
     vec<basevector> missing_genomic_kmers;
     int64_t non_genomic_kmers = 0, N50_unibase = 0;
     if ( K == 96 )
     {    UnibaseSummaryStats<96>( unibases, genome, missing_genomic_kmers, 
               non_genomic_kmers, N50_unibase );    }
     else ForceAssert( 0 == 1 );
     longlong genome_size = 0;
     for ( size_t t = 0; t < genome.size( ); t++ )
          genome_size += genome[t].size( );
     sout << "missing = " << setiosflags(ios::fixed) << setprecision(1) << 10000.0
          * double( missing_genomic_kmers.size( ) ) / double( 2 * genome_size )
          << "%%, ";
     sout << "non-genomic = " << setiosflags(ios::fixed) << setprecision(1)
          << 10000.0 * double(non_genomic_kmers) / double( 2 * genome_size ) 
          << "%%, ";
     // sout << "N50 unibase = " << N50_unibase << endl;

     // Run AssemblyAccuracy.

     String unibases_fn = run_dir + "/" + READS + "." + UNIBASES + kK
          + ".onedir.tmp.fasta";
     {    Ofstream( out, unibases_fn );
          for ( int u = 0; u < (int) unibases.size( ); u++ )
               if ( u <= to_rc[u] ) unibases[u].Print( out, u );    }
     sout << "error rate ="
          << LineOfOutput( "AssemblyAccuracy" + ARG(ASSEMBLY, unibases_fn)
			   + ARG(REF, data_dir + "/../genome_extended.fasta")
			   + ARG(NUM_THREADS, NUM_THREADS)
			   + ARG(CHUNK_SIZE, "1G") + ARG(NH, True) 
			   + " | grep \"error rate\"" )
       .After( "=" ).Before( " x" ) << "%%, ";
     Remove(unibases_fn);

     // Run HyperToReftigs.

     sout << LineOfOutput( "HyperToReftigs" + ARGC(PRE) + ARGC(DATA) + ARGC(RUN)
          + ARGC(K) + ARGC(READS) + ARG(PRINT_REFTIGS, False) + ARG(NH, True)
          + " | grep N50" ) << "\n";

     // Mirror results in standard output.

     sout.close( );
     CpAppend( run_dir + "/UnipathBabyStats.out", cout );    }
