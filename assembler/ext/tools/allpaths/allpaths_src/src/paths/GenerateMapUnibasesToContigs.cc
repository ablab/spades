/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "CoreTools.h"
#include "MainTools.h"

#include "Basevector.h"
#include "Intvector.h"
#include "paths/MapUnibasesToContigs.h"

/**
 * GenerateMapUnibasesToContigs
 *
 * Call MapUnibasesToContigs to generate a unibases to contigs map,
 * which is saved as a feudal UInt64VecVec object.
 *
 * READS: the unibases, ie <READS>.unibases.k<K> in <RUN>
 * ASSEMBLY: the contigs, ie <ASSEMBLY>.contigs.fastb in <SUBDIR>
 * OUTPUT: UInt64VecVec output, ie <ASSEMBLY>.u2c
 */ 
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_Int( K );
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_String_OrDefault( READS, "all_reads" );
  CommandArgument_String_OrDefault( SUBDIR, "test" );
  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
    "Number of threads to use (use all available processors if set to 0)");
  CommandArgument_String_OrDefault(ASSEMBLY, "linear_scaffolds0.clean.patched");
  EndCommandArguments;
  
  // Thread control (MapUnibasesToContigs calls SearchFastb which uses OMP)
  
  NUM_THREADS = configNumThreads(NUM_THREADS);
  omp_set_num_threads( NUM_THREADS );


  String run_dir = PRE + "/" + DATA + "/" + RUN;
  String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;
  
  String str_kK = "k" + ToString( K );
  String uni_file = run_dir + "/" + READS + ".unibases." + str_kK;
  String head_file = sub_dir + "/" + ASSEMBLY;
  String cg_file =  head_file + ".contigs.fastb";
  
  bool CACHE = false;
  MapUnibasesToContigs( K, uni_file, cg_file, head_file, &cout );
  
}
