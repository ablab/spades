// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

// Program: GenomeToPaths
// 
// Given sequence for a genome as a fastb file, generate one
// <KmerPath> for each record in the fastb file (i.e. for each <genome part>).
// This KmerPath may contain adjacent and concatenatable intervals (i.e. not
// be the "minimal" / canonical / most-compressed representation of this sequence
// of kmers).
//
// Command-line arguments: HEAD, K.
//
// Input:
// HEAD.fastb = genome fastb file taken as input
//
// Output:
// HEAD.paths.k* = output paths for genome, where * is K.

#include "Basevector.h"
#include "MainTools.h"
#include "paths/KmerPath.h"
#include "paths/ReadsToPathsCoreX.h"

int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String_OrDefault(HEAD, "genome");
  CommandArgument_Int(K);
  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
    "Number of threads to use (use all available processors if set to 0)");
  CommandArgument_Bool_OrDefault_Doc(CANONICALIZE, False,
				     "canonicalize the genome paths" );
  EndCommandArguments;
  
  // Thread control
   
  NUM_THREADS = configNumThreads(NUM_THREADS);

  cout << Date( ) << ": GenomeToPaths is pathing the genome..." << endl;
  vecbasevector genome( HEAD + ".fastb" );
  vecKmerPath paths;
  ReadsToPathsCoreY( genome, K, paths, HEAD + ".GenomeToPaths", NUM_THREADS );
  if ( CANONICALIZE )
    for ( size_t i = 0; i < paths.size(); i++ )
      paths[i].Canonicalize();
  
  paths.WriteAll( HEAD + ".paths.k" + ToString(K) );
  cout << Date( ) << ": Done!" << endl;
}
