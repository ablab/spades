///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * AlignUnibasesToRef
 *
 * Align a set of unibases to a reference genome.  This is mostly a wrapper to
 * GetAlignsFast, but it also does some basic analysis at the end, after the
 * output file has been written (so if the analysis crashes, we're ok.)
 *
 * This module is designed to replace FindUnipathLocs in the RunAllPathsLG
 * pipeline.
 *
 *
 * INPUT: <RUN>/<READS>.unibases.k<K>, <DATA>/genome.lookup
 * OUTPUT: <RUN>/<READS>.unipaths.k<K>.locs
 *
 *
 * Josh Burton
 * April 2010
 *
 ******************************************************************************/

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "MainTools.h"
#include "Intvector.h" // VecIntVec
#include "lookup/LookAlign.h"
#include "math/Functions.h" // N50
#include "paths/reporting/ReftigUtils.h" // GetAlignsFast
#include "paths/simulation/Placement.h"



void
AnalyzeUnibaseAlignments( ostream & out,
			  const VecPlacementVec & placements,
			  const String & unibases_file,
			  const size_t n_unipaths,
			  const int K )
{
  // Maximum copy number that we track (affects output reporting only)
  // We leave this at 10, because KEEP_BEST=10 in QueryLookupTable
  static const size_t MAX_CN = 10;
  
  
  // Load in unibases to find their lengths.
  vec<int> ulen( n_unipaths );
  {
    vecbasevector unibases( unibases_file );
    ForceAssertEq( n_unipaths, unibases.size() );
    for ( size_t u = 0; u < n_unipaths; u++ )
      ulen[u] = unibases[u].size();
  }
  
  
  
  // Make a histogram of the copy numbers of unipaths, keyed by their lengths.
  // Note that we convert here from uniBASE to uniPATH lengths.
  VecIntVec CNs_histo( MAX_CN + 1 );
  size_t n_placements = 0;
  for ( size_t u = 0; u < n_unipaths; u++ ) {
    size_t CN = placements[u].size();
    n_placements += CN;
    if ( CN > MAX_CN ) CN = MAX_CN;
    CNs_histo[CN].push_back( ulen[u] - K + 1 );
  }
  
  
  // Make a chart header line.
  vec< vec<String> > rows;
  vec<String> row1, row2;
  row1.push_back( "CN", "N unibases", "Total unipath len", "N50 unipath len" );
  row2.push_back( "--", "----------", "-----------------", "---------------" );
  rows.push_back( row1, row2 );
  
  // For each element in the histogram of copy numbers, make a line of info.
  for ( size_t i = 0; i < CNs_histo.size(); i++ ) {
    
    // We have the convert the SerfVec into a vec, because N50 isn't defined
    // for SerfVecs.  Sigh.
    vec<int> lengths( CNs_histo[i].size() );
    for ( size_t j = 0; j < lengths.size(); j++ )
      lengths[j] = CNs_histo[i][j];
    
    String i_str = ToString(i);
    if ( i == MAX_CN ) i_str += "+";
    size_t N50_len = lengths.empty() ? 0 : N50( lengths );
    
    vec<String> row;
    row.push_back( i_str,
		   ToString( lengths.size() ),
		   ToString( Sum( lengths ) ),
		   ToString( N50_len ) );
    rows.push_back(row);
  }
  
  
  // Write the report!
  out << "\n\nREPORT ON ALIGNMENTS" << endl;
  out << n_unipaths << " unibases produced " << n_placements
      << " alignments to reference." << endl << endl;
  
  out << endl;
  out << "Distribution of true unibase copy numbers:" << endl;
  PrintTabular( out, rows, 4, "rrrr" );
  
  out << endl;
}
  




int main( int argc, char *argv[] )
{
  RunTime();
  
  BeginCommandArguments;
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_Int( K );
  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
    "Number of threads to use (use all available processors if set to 0)");
  CommandArgument_String_OrDefault( READS, "reads" );
  CommandArgument_Bool_OrDefault( REPORT, True );
  EndCommandArguments;

  // Thread control (OMP used by SearchFastb used by GetAlignsFast)
   
  NUM_THREADS = configNumThreads(NUM_THREADS);
  omp_set_num_threads( NUM_THREADS );

  // Find the names of the files to be used.
  String data_dir = PRE + "/" + DATA;
  String run_dir = PRE + "/" + DATA + "/" + RUN;
  String kK = ".k" + ToString(K);
  
  String genome_file = data_dir + "/genome.lookup";
  String unibases_file = run_dir + "/" + READS + ".unibases" + kK;
  String aligns_file = run_dir + "/" + READS + ".unipaths" + kK + ".locs";
  
  // Create temp dir
  Mkdir777(run_dir + "/tmp");
  
  // Call GetAlignsFast.  It gets aligns fast.
  // It uses aligns_file as a temporary storage location for alignments.
  cout << Date() << ": AlignUnibasesToRef is calling GetAlignsFast" << endl;
  int K_search = K <= 144 ? K : 144; 
  vec<look_align> aligns;
  GetAlignsFast( K_search, unibases_file, genome_file, aligns_file, aligns, false, run_dir + "/tmp" );
  
  
  // Convert the look_aligns to Placements.
  cout << Date() << ": Converting look_aligns to Placements" << endl;
  size_t n_unipaths = MastervecFileObjectCount( unibases_file );
  VecPlacementVec placements( n_unipaths );
  
  for ( size_t i = 0; i < aligns.size(); i++ ) {
    placement p( aligns[i].target_id, aligns[i].pos2(), aligns[i].Pos2(), aligns[i].Rc1() );
    placements[ aligns[i].query_id ].push_back( p );
  }
  Destroy( aligns );
  
  
  // Write placement file.  This overwrites the aligns_file which we created
  // earlier during GetAlignsFast.
  cout << Date() << ": Writing VecPlacementVec to disk" << endl;
  placements.WriteAll( aligns_file );
  cout << Date() << ": Output file is now in place - no need to re-run AlignUnibasesToRef if analysis crashes" << endl;
  
  
  // Write a basic report on these alignments.
  cout << Date() << ": Performing basic analysis" << endl;
  AnalyzeUnibaseAlignments( cout, placements, unibases_file, n_unipaths, K );
  
  cout << Date() << ": Done with AlignUnibasesToRef!" << endl;
}
