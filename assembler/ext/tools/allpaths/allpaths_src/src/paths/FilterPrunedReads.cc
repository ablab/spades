/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2010) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/* FilterPrunedReads
 *
 * Remove from a dataset of reads/quals/paths all of the reads containing kmers
 * not in a set of unipaths.
 *
 * This module is mostly a wrapper for the function MarkReadsWithNovelKmers.
 * With default parameters, it removes from the filled reads (the output of
 * FillFragments) all reads containing kmers that have been shaved off during
 * the extending/patching/etc. of unipaths.
 *
 *
 * Required input files:
 * <PRE>/<DATA>/<RUN>/<READS_IN>.paths.k<K>
 * <PRE>/<DATA>/<RUN>/<UNIPATHS>.unipathsdb.k<K>
 *
 * Optional input files: (if these exist, will create corresponding outputs)
 * <PRE>/<DATA>/<RUN>/<READS_IN>.{fastb,qualb,pairs}
 *
 * Output files:
 * <PRE>/<DATA>/<RUN>/<READS_OUT>.paths.k<K>
 * <PRE>/<DATA>/<RUN>/<READS_OUT>.{fastb,qualb,pairs} (may not be created)
 *
 *
 * Josh Burton
 * June 2010
 *
 *****************************************************************************/



#include "MainTools.h"
#include "String.h"
#include "Basevector.h"
#include "PairsManager.h"
#include "Qualvector.h"
#include "paths/KmerPath.h"
#include "util/ReadTracker.h"
#include "graph/Digraph.h"



int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_Int( K );
  CommandArgument_String_OrDefault( READS_IN, "filled_reads_ext" );
  CommandArgument_String_OrDefault_Doc( UNIPATHS, "extended", "Make sure READS and UNIPATHS use the same kmer numbering system!" );
  CommandArgument_Bool_OrDefault_Doc( USE_GRAPH, False, "use information from unipath adjacency graph" );
  CommandArgument_String_OrDefault( READS_OUT, "filled_reads_filt" );
  EndCommandArguments;
  
  
  // Filenames.
  String run_dir = PRE + "/" + DATA + "/" + RUN;
  String  reads_in = run_dir + "/" + READS_IN;
  String reads_out = run_dir + "/" + READS_OUT;
  String bases_in_file = reads_in + ".fastb";
  String quals_in_file = reads_in + ".qualb";
  String pairs_in_file = reads_in + ".pairs";
  String paths_in_file = reads_in + ".paths.k" + ToString(K);
  String unipathsdb_file = run_dir + "/" + UNIPATHS + ".unipathsdb.k" + ToString(K);
  String uadj_graph_file = run_dir + "/" + UNIPATHS + ".unipath_adjgraph.k" + ToString(K);
  
  if ( USE_GRAPH )
    ForceAssert( IsRegularFile( uadj_graph_file ) );

  // Load files (though don't load all the files yet - to save memory.)
  cout << Date() << ": Loading files..." << endl;
  vecKmerPath paths( paths_in_file );
  size_t n_reads = paths.size();
  
  vec<tagged_rpint> unipathsdb;
  BinaryRead3( unipathsdb_file, unipathsdb );
  
  digraph AG;
  if ( USE_GRAPH ){
    cout << Date() << ": reading adjacency graph" << endl;
    BinaryRead( uadj_graph_file, AG );
  }

  // Main algorithm!
  // This whole module is a wrapper for MarkReadsWithNovelKmers.
  vec<Bool> novel;
  if ( ! USE_GRAPH ){
    cout << Date() << ": Calling MarkReadsWithNovelKmers (main algorithm without adjacency graph)" << endl;
    MarkReadsWithNovelKmers( paths, unipathsdb, novel );
  }else{
    cout << Date() << ": Calling MarkReadsWithNovelKmers2 (main algorithm using adjacency graph)" << endl;
    MarkReadsWithNovelKmers2( paths, unipathsdb, AG, novel );
  }
  Destroy( unipathsdb );

  
  
  // Report results of main algorithm.
  // NOTE: If too many reads have been marked as novel, you've probably
  // mismatched the kmer numbering systems.
  size_t n_novel = Sum( novel );
  double pct_novel = 100.0 * n_novel / double( n_reads );
  cout << Date() << ": Number of reads containing kmers not in the unipaths: " << n_novel << " / " << n_reads << " = " << pct_novel << "%" << endl;
  
  // Remove the reads containing novel kmers from the dataset.
  cout << Date() << ": Removing these reads from the dataset" << endl;
  paths.EraseIf( novel );
  size_t n_reads_new = paths.size();
  
  // Write output paths file.
  cout << Date() << ": Writing output file: paths" << endl;
  paths.WriteAll( reads_out + ".paths.k" + ToString(K) );
  Destroy( paths );
  
  
  // Remove reads from the bases (if the file exists.)
  if ( IsRegularFile( bases_in_file ) ) {
    cout << Date() << ": Writing output file: fastb" << endl;
    vecbasevector bases( bases_in_file );
    ForceAssertEq( bases.size(), n_reads );
    bases.EraseIf( novel );
    ForceAssertEq( bases.size(), n_reads_new );
    bases.WriteAll( reads_out + ".fastb" );
  }
  else cout << "\tNot writing an output fastb file because the corresponding input does not exist: " << bases_in_file << endl;
  
  // Remove reads from the quals (if the file exists.)
  if ( IsRegularFile( quals_in_file ) ) {
    cout << Date() << ": Writing output file: qualb" << endl;
    vecqualvector quals( quals_in_file );
    ForceAssertEq( quals.size(), n_reads );
    quals.EraseIf( novel );
    ForceAssertEq( quals.size(), n_reads_new );
    quals.WriteAll( reads_out + ".qualb" );
  }
  else cout << "\tNot writing an output qualb file because the corresponding input does not exist: " << quals_in_file << endl;
  
  // Remove reads from the pairs (if the file exists.)
  if ( IsRegularFile( pairs_in_file ) ) {
    cout << Date() << ": Writing output file: pairs" << endl;
    PairsManager pairs( pairs_in_file );
    ForceAssertEq( pairs.nReads(), n_reads );
    pairs.removeReads( novel, true );
    ForceAssertEq( pairs.nReads(), n_reads_new );
    pairs.Write( reads_out + ".pairs" );
  }
  else cout << "\tNot writing an output pairs file because the corresponding input does not exist: " << pairs_in_file << endl;
  
  
  // Write a ReadTracker.
  cout << Date() << ": Writing output file: readtrack" << endl;
  ReadTracker rt;
  rt.AddReadSet( reads_in, novel );
  rt.Dump( reads_out ); // creates the file <reads_out>.readtrack
  
  
  
  // Done!
  cout << Date() << ": Done with FilterPrunedReads!" << endl;
  return 0;
}
