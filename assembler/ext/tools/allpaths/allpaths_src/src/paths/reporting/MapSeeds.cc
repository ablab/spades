///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MapSeeds.  Generate a file "seeds.map" showing the true locations of seeds
// on a reference sequence. 

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "Basevector.h"
#include "MainTools.h"
#include "VecUtilities.h" // SortSync
#include "lookup/LookAlign.h"
#include "math/Functions.h" // N50
#include "paths/KmerPath.h"
#include "paths/UnipathNhoodCommon.h" // digraphE<fsepdev>
#include "paths/simulation/Placement.h"




// MakeDepend: dependency QueryLookupTable

int main( int argc, char *argv[] )
{
  
  RunTime();
  
  BeginCommandArguments;
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_String_OrDefault( SUBDIR, "test" );
  CommandArgument_Int(K);
  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
    "Number of threads to use (use all available processors if set to 0)");
  CommandArgument_Bool_OrDefault_Doc(KEEP_ALIGNS, False, 
    "Optionally keep aligns file (saved in ..<SUBDIR>/seeds.qltout)");
  CommandArgument_String_OrDefault( READS, "reads" );
  EndCommandArguments;

  // Thread control
   
  NUM_THREADS = configNumThreads(NUM_THREADS);
  omp_set_num_threads( NUM_THREADS );
  
  cout << Date() << ": Loading files" << endl;
  String data_dir = PRE + "/" + DATA;
  String run_dir = PRE + "/" + DATA + "/" + RUN;
  String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;
  String temp_dir = sub_dir + "/tmp";
  Mkdir777(temp_dir);
  
  String KS = ToString(K);
  String unibases_file = run_dir + "/" + READS + ".unibases.k" + KS;
  size_t n_unipaths = MastervecFileObjectCount( unibases_file );
  
  // Load unipath link graph to get a mapping of unipath IDs to component IDs.
  vec<int> u_components( n_unipaths, -1 );
  {
    digraphE<fsepdev> LG;
    BinaryRead( sub_dir + "/unipath_link_graph.seeds.k" + KS, LG );
    ForceAssertEq( n_unipaths, (size_t)LG.N() );
    
    vec< vec<int> > components;
    LG.Components( components );
    for ( size_t i = 0; i < components.size(); i++ )
      for ( size_t j = 0; j < components[i].size(); j++ )
	u_components[ components[i][j] ] = i;
  }
  
  // Load seeds.ids.
  Ifstream( in, sub_dir + "/seeds.ids" );
  vec<size_t> ids;
  int n_seeds;
  in >> n_seeds;
  ids.reserve( n_seeds );
  ids.ReadFromTextStream( in );
  if ( ids.nonempty( ) ) ForceAssertLt( ids.back(), n_unipaths );
  
  // Open seeds.map for printing.
  Ofstream( out, sub_dir + "/seeds.map" );
  command.PrintTheCommandPretty( out );
  
  
  // Prepare for parallel calls to QueryLookupTable.
  int N = ids.size();
  for ( size_t j = 0; j < NUM_THREADS; j++ ) {
    Ofstream( idsout, temp_dir + "/seeds.ids.alt." + ToString(j) );
    for ( int i = 0; i < N; i++ )
      if ( i % NUM_THREADS == j ) idsout << ids[i] << "\n";
  }
  
  // Call QueryLookupTable in parallel to align the seeds to reference.
  cout << Date() << ": Aligning seeds to reference" << endl;
#pragma omp parallel for
  for ( size_t j = 0; j < NUM_THREADS; j++ ) {
    SystemSucceed( "QueryLookupTable K=12 MM=12 MC=0.15"
		   + ARG(SEQS, unibases_file)
		   + ARG(SEQS_IS_FASTB, True) + ARG(AI, True)
		   + ARG(L, data_dir + "/genome.lookup") + ARG(VISUAL, True)
		   + ARG(SEQS_TO_PROCESS, "@" + temp_dir + "/seeds.ids.alt." + ToString(j) )
		   + ARG(PARSEABLE, True) + ARG(TMP_DIR, temp_dir)
		   + " > " + temp_dir + "/seeds.aligns." + ToString(j) );
  }
  
  // Load (and in case save as a single file) the alignments.
  cout << Date() << ": Loading alignments" << endl;
  vec<look_align> aligns;
  for ( size_t j = 0; j < NUM_THREADS; j++ ) {
    vec<look_align> alignsj;
    LoadLookAligns( temp_dir + "/seeds.aligns." + ToString(j), alignsj );
    aligns.append(alignsj);
  }
  if ( KEEP_ALIGNS ) {
    String qlt_file = sub_dir + "/seeds.qltout";
    cout << Date( ) << ": Saving alignments" << endl;
    WriteLookAligns( qlt_file, aligns );
  }
  
  // Convert the alignments to a VecPlacementVec.
  cout << Date() << ": Converting alignments from look_aligns to placements" << endl;
  VecPlacementVec all_placements( n_unipaths );
  for ( size_t i = 0; i < aligns.size(); i++ ) {
    const look_align & a = aligns[i];
    placement p( a.target_id, a.pos2(), a.Pos2(), a.Rc1() );
    all_placements[ a.query_id ].push_back( p );
  }
  Destroy( aligns );
  
  
  // Filter the alignments by seed copy number.
  cout << Date() << ": Filtering and sorting alignments" << endl;
  vec<placement> placements;
  placements.reserve( N );
  size_t n_not_placed = 0, n_uniquely_placed = 0, n_multiply_placed = 0;
  
  for ( size_t i = 0; i < ids.size(); i++ ) {
    size_t uid = ids[i];
    
    size_t n_places = all_placements[uid].size();
    
    if ( n_places == 0 ) { // unplaced seed
      n_not_placed++;
      placements.push_back( placement(-1,-1,-1,-1) );
    }
    else if ( n_places == 1 ) { // singly placed seed (most)
      n_uniquely_placed++;
      placements.push_back( all_placements[uid][0] );
    }
    else { // multiply placed seed
      n_multiply_placed++;
      placements.push_back( placement( -n_places, -1, -1, -1 ) );
    }
  }
  Destroy( all_placements );
  
  
  // Report seed placements to seeds.map.
  out << "TOTAL NUMBER OF SEEDS: " << N << endl;
  out << "Number of seeds with..." << endl;
  out << "... no placement:\t\t" << n_not_placed << endl;
  out << "... a unique placement:\t\t" << n_uniquely_placed << endl;
  out << "... multiple placements:\t" << n_multiply_placed << endl;
  out << endl;
  
  
  // Sort the placements by reference location.
  vec<int> index( N, vec<int>::IDENTITY );
  SortSync( placements, index );
  
  
  // Analyze the sizes of the seeds and of the gaps between them.
  cout << Date() << ": Analyzing seed and gap sizes" << endl;
  vec<size_t> seed_sizes, gap_sizes;
  int last_seed_gID = 0;
  int last_seed_end = -1;
  
  for ( int i = 0; i < N; i++ ) {
    const placement & p = placements[i];
    if ( p.GenomeId() < 0 ) continue;
    
    seed_sizes.push_back( p.Pos() - p.pos() );
    
    // Find the gap between this and the previous seed.
    // Note that we do not include gaps at the starts/ends of reference contigs.
    if ( last_seed_gID == p.GenomeId() && last_seed_end != -1 )
      if ( p.pos() > last_seed_end ) // avoid a gap size of 0
	gap_sizes.push_back( p.pos() - last_seed_end );
    
    last_seed_gID = p.GenomeId();
    last_seed_end = p.Pos();
  }
  
  // Print a report on seed/gap sizes to seeds.map.
  out << "Seed unipath sizes:\t\t" << BigSum( seed_sizes ) / double(N)
      << " (avg)";
  if ( seed_sizes.nonempty( ) ) out << ", " << N50( seed_sizes ) << " (N50)";
  out << endl;
  out << "Sizes of gaps between seeds:\t" << BigSum( gap_sizes ) / double(N)
      << " (avg)";
  if ( gap_sizes.nonempty( ) ) out << ", " << N50( gap_sizes ) << " (N50)";
  out << endl;
  out << endl;
  
  
  
  // Make a chart of seed placements on reference.
  cout << Date() << ": Writing output to file <SUBDIR>/seeds.map" << endl;
  vecKmerPath unipaths( run_dir + "/" + READS + ".unipaths.k" + KS );
  vec< vec<String> > rows;
  vec<String> title1, title2;
  title1.push_back( "Unipath ID", "Seed ID", "Component", "Length", "Ref location" );
  title2.push_back( "----------", "-------", "---------", "------", "------------" );
  rows.push_back( title1, title2 );
  for ( int i = 0; i < N; i++ ) {
    
    int seed_id = index[i];
    int uid = ids[seed_id];
    vec<String> row;
    
    String loc = "unplaced";
    const placement & p = placements[i];
    
    if ( p.GenomeId() >= 0 ) { // singly placed seed
      ostringstream s;
      s << p;
      loc = s.str();
    }
    else if ( p.GenomeId() < -1 ) { // multiply placed seed
      loc = "N PLACEMENTS: " + ToString( -p.GenomeId() );
    }
    
    row.push_back( ToString(uid), ToString(seed_id),
		   ToString( u_components[uid] ),
		   ToString( unipaths[uid].KmerCount() ),
		   loc );
    rows.push_back(row);
  }
  
  
  // Print the chart of seed placements to seeds.map.
  PrintTabular( out, rows, 2, "rrrrl" );
  
  
  
  
  // Cleanup.
  for ( size_t j = 0; j < NUM_THREADS; j++ ) {
    Remove( temp_dir + "/seeds.ids.alt." + ToString(j) );
    Remove( temp_dir + "/seeds.aligns." + ToString(j) );
  }
  cout << Date() << ": Done with MapSeeds" << endl;
}
