///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"

#include "SeqInterval.h"
#include "VecTemplate.h"
#include "VecUtilities.h"
#include "paths/Alignlet.h"
#include "paths/CommonPatherCore.h"
#include "paths/KmerPath.h"
#include "paths/PdfEntry.h"
#include "system/MiscUtil.h"
#include "util/RunCommand.h"

/**
 * RemoveHighCNAligns
 *
 * Remove alignments of reads onto regions of contigs having copy
 * number bigger than PLOIDY_MULTIPLIER * ploidy.
 *
 * REMARK: this actually works by only resetting the indexes of the
 * deleted to -1 (the alignments are not changed).
 *
 * UNIPATHS: relative to run_dir (it needs copy number estimates)
 * ALIGNS_IN: input aligns
 * ALIGNS_OUT: output aligns
 * SCAFFOLDS_IN: assembly contigs
 * RATIO_TO_EXCLUDE: if a read extends in an unique region, do not remove it
 * PLOIDY_MULTIPLIER: see above
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  BeginCommandArguments;
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_String_OrDefault( SUBDIR, "test" );
  CommandArgument_Int_OrDefault( K, 96 );
  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
    "Number of threads to use (use all available processors if set to 0)");
  CommandArgument_String_OrDefault( UNIPATHS, "all_reads" );
  CommandArgument_String_OrDefault( ALIGNS_IN, "scaffold_reads" );
  CommandArgument_String_OrDefault( ALIGNS_OUT, "scaffold_reads_filtered" );
  CommandArgument_String_OrDefault( SCAFFOLDS_IN, "initial_scaffolds" );
  CommandArgument_Double_OrDefault( RATIO_TO_EXCLUDE, 0.5 );
  CommandArgument_Int_OrDefault( PLOIDY_MULTIPLIER, 1 );
  CommandArgument_Bool_OrDefault( FORCE, False );
  EndCommandArguments;

  // Thread control
  
  NUM_THREADS = configNumThreads(NUM_THREADS);

  // Dir and file names.
  String data_dir = PRE + "/" + DATA;
  String run_dir = data_dir + "/" + RUN;
  String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;
  String out_dir = sub_dir + "/" + SCAFFOLDS_IN + ".tmp";
  String out_reldir = "ASSEMBLIES/" + SUBDIR + "/" + SCAFFOLDS_IN + ".tmp";

  String ploidy_file = run_dir + "/ploidy";

  String strK = "k" + ToString( K );
  String uni_base = run_dir + "/" + UNIPATHS; 
  String copynum_file = uni_base + ".unipaths.predicted_count." + strK;
  String unibases_orig_file = uni_base + ".unibases." + strK;

  String aligns_file = sub_dir + "/" + ALIGNS_IN + ".qltoutlet";
  String index_file =  sub_dir + "/" + ALIGNS_IN + ".qltoutlet.index";
  String contigs_orig_file = sub_dir + "/" + SCAFFOLDS_IN + ".contigs.fasta";

  String index_out_file =  sub_dir + "/" + ALIGNS_OUT + ".qltoutlet.index";

  String unibases_file = out_dir + "/unibases.fastb";
  String contigs_file = out_dir + "/contigs.fastb";
  String unipaths_file = out_dir + "/unibases.paths";
  String cgpaths_file = out_dir + "/contigs.paths";
  String cgpathsdb_full_file = cgpaths_file + "db." + strK;
  String cgpathsdb_big_full_file = cgpaths_file + "db_big." + strK;
  
  Mkpath( out_dir );

  // Link original unibases into out_dir.
  String rel_unibases_file;
  AbsToRelPaths(unibases_orig_file, unibases_file, rel_unibases_file);
  SymlinkForce(rel_unibases_file, unibases_file);

  // Generate contigs fastb.
  if ( FORCE || ! IsRegularFile( contigs_file ) ) {
    cout << Date( ) << ": converting contigs fasta to fastb" << endl;
    vec<fastavector> fasta;
    LoadFromFastaFile( contigs_orig_file, fasta );

    // Warning! ambiguous bases will be randomized.
    GenCharToRandomBaseMapper ran;
    vecbvec contigs;
    contigs.reserve( fasta.size( ) ); 
    for (size_t ii=0; ii<fasta.size( ); ii++)
      contigs.push_back( bvec( fasta[ii].begin( ), fasta[ii].end( ), ran ) );

    contigs.WriteAll( contigs_file );
  }

  // Run CommonPather.
  bool unipaths_missing = ! IsRegularFile( unipaths_file + "." + strK );
  bool cgpaths_missing =  ! IsRegularFile( cgpaths_file + "." + strK );
  if ( FORCE || unipaths_missing || cgpaths_missing ) {
    cout << Date( ) << ": running CommonPather" << endl;
    vec<String> fastb_in;
    fastb_in.push_back( contigs_file );
    fastb_in.push_back( unibases_file );
    
    vec<String> paths_out;
    paths_out.push_back( cgpaths_file );
    paths_out.push_back( unipaths_file );
    
    CommonPather( K, fastb_in, paths_out, NUM_THREADS );
  }

  // Run MakeRcDb (logging within).
  if ( FORCE || ( ! IsRegularFile(cgpathsdb_full_file ) 
      && ! IsRegularFile(cgpathsdb_big_full_file)) ) {
    String theCommand
      = "MakeRcDb K=" + ToString( K )
      + " PRE=" + PRE
      + " DATA=" + DATA
      + " RUN=" + RUN
      + " READS=" + out_reldir + "/contigs";
    RunCommand( theCommand );
  }
  
  // Load.
  int ploidy = FirstLineOfFile( ploidy_file ).Int( );
  
  cout << Date( ) << ": loading aligns" << endl;
  vec<alignlet> aligns;
  BinaryRead3( aligns_file, aligns );
  
  cout << Date( ) << ": loading index" << endl;
  vec<int> index;
  BinaryRead3( index_file, index );
  
  cout << Date( ) << ": loading copy numbers" << endl;
  VecPdfEntryVec copynum( copynum_file );

  cout << Date( ) << ": loading contigs db" << endl;
  ForceAssert( IsRegularFile(cgpathsdb_full_file)||IsRegularFile(cgpathsdb_big_full_file) );
  bool big_db = False;
  if ( IsRegularFile(cgpathsdb_big_full_file) ) big_db = True;
  vec<tagged_rpint> cgpathsdb;
  vec<big_tagged_rpint> cgpathsdb_big;
  if(big_db)
    BinaryRead2(cgpathsdb_big_full_file, cgpathsdb_big);
  else
    BinaryRead2(cgpathsdb_full_file, cgpathsdb);
  cout << Date( ) << ": loading contigs unipaths" << endl;
  vecKmerPath cgpaths( cgpaths_file + "." + strK );

  cout << Date( ) << ": loading unipaths" << endl;
  vecKmerPath unipaths( unipaths_file + "." + strK );
  
  // Intervals of coverage on contigs (reserve memory).
  vec<seq_interval> cov;
  size_t n_segs = 0;
  for (size_t ii=0; ii<unipaths.size( ); ii++)
    n_segs += unipaths[ii].NSegments( );
  cov.reserve( n_segs );

  // Place all unipaths onto the contigs (non trivial, but purely mechanical).
  cout << Date( ) << ": identifying high copy number regions" << endl;
  for (size_t unipath_id=0; unipath_id<unipaths.size( ); unipath_id++) {
    const KmerPath &uni = unipaths[unipath_id];
    int uni_klen = uni.TotalLength( );
    
    // Copy number estimate for this unipath (naive: highest probability).
    int uni_cn = copynum[unipath_id][0].first;
    double uni_cnprob = copynum[unipath_id][0].second;
    for (size_t ii=1; ii<copynum[unipath_id].size( ); ii++) {
      if ( copynum[unipath_id][ii].second > uni_cnprob ) {
	uni_cn = copynum[unipath_id][ii].first;
	uni_cnprob = copynum[unipath_id][ii].second;	
      }
    }
    
    // Loop over all segments in the unipath.
    for (int segment_id=0; segment_id<uni.NSegments( ); segment_id++) {
      vec<longlong> locs;
      if(big_db)
	Contains( cgpathsdb_big, uni.Segment( segment_id ), locs );
      else
	Contains( cgpathsdb, uni.Segment( segment_id ), locs );
      longlong seg_start = uni.Segment( segment_id ).Start( );
      longlong seg_len = uni.Segment( segment_id ).Length( );
      
      // Loop over all placements for this segment.
      for (int ii=0; ii<locs.isize( ); ii++) {
	// make it compatable with both tagged_rprint and big_tagged_rprint
	//const tagged_rpint& t_rpint = cgpathsdb[locs[ii]];
	if (big_db && cgpathsdb_big[locs[ii]].Rc( )) continue;
	if (!big_db && cgpathsdb[locs[ii]].Rc( )) continue;

	int contig_id = big_db ? cgpathsdb_big[locs[ii]].ReadId( ): cgpathsdb[locs[ii]].ReadId( );
	int contig_pos =big_db ? cgpathsdb_big[locs[ii]].PathPos(): cgpathsdb[locs[ii]].PathPos( );
	longlong contig_len = cgpaths[contig_id].TotalLength( );
	const KmerPathInterval &kpi = cgpaths[contig_id].Segment( contig_pos );
	ForceAssert( kpi.Overlaps( uni.Segment( segment_id ) ) );

	longlong begin = Max( (longlong)0, seg_start - kpi.Start( ) );
	for (int jj=0; jj<contig_pos; jj++)
	  begin += cgpaths[contig_id].Segment( jj ).Length( );
	longlong end = Min( begin + seg_len, contig_len );
	
	// Add info.
	seq_interval si( uni_cn, contig_id, (int)begin, (int)end );
	cov.push_back( si );
      }
    }
  }
  
  // Compactify cov (merge adjacent intervals with the same copy number).
  sort( cov.begin( ), cov.end( ) );
  {
    vec<seq_interval> new_cov;
    new_cov.reserve( cov.size( ) );
    
    if ( cov.size( ) > 0 ) new_cov.push_back( cov[0] );
    for (size_t ii=1; ii<cov.size( ); ii++) {
      if ( new_cov.back( ).SeqId( ) != cov[ii].SeqId( ) ||
	   new_cov.back( ).IntervalId( ) != cov[ii].IntervalId( ) ||
	   new_cov.back( ).End( ) < cov[ii].Begin( ) ) {
	new_cov.push_back( cov[ii] );
	continue;
      }
      new_cov[new_cov.size( )-1].SetEnd( cov[ii].End( ) );
    }
    
    swap( new_cov, cov );
  }
  
  // Define bad regions: intervals on contig with coverage > ploidy.
  vec<seq_interval> bad;
  bad.reserve( cov.size( ) );
  for (size_t ii=0; ii<cov.size( ); ii++)
    if ( cov[ii].IntervalId( ) > PLOIDY_MULTIPLIER * ploidy )
      bad.push_back( cov[ii] );

  // Translate bad in base coordinates (from kmer coordinates).
  for (size_t ii=0; ii<bad.size( ); ii++)
    bad[ii].SetEnd( bad[ii].End( ) + K - 1 );
  
  // Update index (tally totals).
  cout << Date( ) << ": resetting index of alignments" << endl;
  longlong n_aligned = 0;
  longlong n_reset = 0;
  vec<seq_interval>::iterator it;
  for( size_t ii=0; ii<index.size( ); ii++) {
    if ( index[ii] < 0 ) continue;
    n_aligned++;
    
    const alignlet &al = aligns[ index[ii] ];
    seq_interval si( 0, al.TargetId( ), al.pos2( ), al.Pos2( ) );
    it = lower_bound( bad.begin( ), bad.end( ), si );
    if ( it != bad.end( ) ) {
      double overlap = it->HasAmountOfOverlapWith( si );
      if ( overlap > double( al.Pos2( ) - al.pos2( ) ) * RATIO_TO_EXCLUDE ) {
	index[ii] = -1;
	n_reset++;
	continue;
      }
    }
    if ( it != bad.begin( ) ) {
      it--;
      double overlap = it->HasAmountOfOverlapWith( si );
      if ( overlap > double( al.Pos2( ) - al.pos2( ) ) * RATIO_TO_EXCLUDE ) {
	index[ii] = -1;
	n_reset++;
	continue;
      }
    }
  }

  // Save.
  cout << Date( ) << ": saving" << endl;
  BinaryWrite3( index_out_file, index );

  // Compute and print some stats.
  String cn_threshold = ToString( PLOIDY_MULTIPLIER * ploidy );

  longlong total_cglen = 0;
  for (size_t ii=0; ii<cgpaths.size( ); ii++)
    total_cglen += cgpaths[ii].TotalLength( );
  
  longlong total_badlen = 0;
  for (size_t ii=0; ii<bad.size( ); ii++)
    total_badlen += bad[ii].Length( );

  String str_ratio_cg = "na";
  if ( total_cglen > 0 ) {
    double ratio = 100.0 * SafeQuotient( total_badlen, total_cglen );
    str_ratio_cg = ToString( ratio, 2 ) + "%";
  }

  String str_ratio_al = "na";
  if ( n_aligned > 0 ) {
    double ratio = 100.0 * SafeQuotient( n_reset, n_aligned );
    str_ratio_al = ToString( ratio, 2 ) + "%";
  }
  
  vec< vec<String> > table;
  vec<String> aline;
  
  aline.push_back( "Total bases in input" );
  aline.push_back( ToString( total_cglen ) );
  table.push_back( aline );
  aline.clear( );

  aline.push_back( "Bases having copy number >" + cn_threshold );
  aline.push_back( ToString( total_badlen ) );
  table.push_back( aline );
  aline.clear( );
  
  aline.push_back( "Fraction of high copy number bases" );
  aline.push_back( str_ratio_cg );
  table.push_back( aline );
  aline.clear( );

  aline.push_back( "Aligns in input" );
  aline.push_back( ToString( aligns.size( ) ) );
  table.push_back( aline );
  aline.clear( );

  aline.push_back( "Valid aligns in input" );
  aline.push_back( ToString( n_aligned ) );
  table.push_back( aline );
  aline.clear( );
  
  aline.push_back( "Aligns deleted" );
  aline.push_back( ToString( n_reset ) );
  table.push_back( aline );
  aline.clear( );
  
  aline.push_back( "Fraction deleted" );
  aline.push_back( str_ratio_al );
  table.push_back( aline );
  aline.clear( );
  
  cout << "\n";
  PrintTabular( cout, table, 3, "rl" );
  cout << "\n";

  // Done.
  cout << Date( ) << ": done" << endl;
  
}


