///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "Fastavector.h"
#include "FetchReads.h"
#include "PairsManager.h"
#include "Superb.h"
#include "graph/Digraph.h"
#include "math/NStatsTools.h"
#include "paths/Alignlet.h"
#include "paths/Alignlets2ReadLocs.h"
#include "paths/BuildScaffoldGraph.h"
#include "paths/SaveScaffoldGraph.h"
#include "paths/ScaffoldsUtils.h"
#include "paths/Sepdev.h"
#include "paths/reporting/CLinkBundle.h"
#include "util/RunCommand.h"
#include "util/SearchFastb2Core.h"

/**
 * AlignReadsToContigs
 */
void AlignReadsToContigs( const int K,
			  const String &tmp_dir,
			  const String &reads_fastb_file,
			  const vec<fastavector> &contigs,
			  vec<alignlet> &aligns,
			  vec<int> &index,
			  ostream &log,
			  bool VERBOSE )
{
  log << Date( ) << ": prepare contigs for realignment" << endl;

  temp_file tmp_contigs_fasta( tmp_dir + "/tmp_contigs_fasta_XXXXXX" );
  ofstream cg_out( tmp_contigs_fasta.c_str( ) );
  for (size_t ii=0; ii<contigs.size( ); ii++)
    contigs[ii].Print( cg_out, "contig_" + ToString( ii ) );
  cg_out.close();
    
  temp_file tmp_contigs_fastb( tmp_dir + "/tmp_contigs_fastb_XXXXXX" );
  vecbasevector target;
  FetchReads( target, 0, tmp_contigs_fasta.c_str( ) );
  target.WriteAll( tmp_contigs_fastb );

  log << Date( ) << ": SearchFastb2 - parallel align reads" << endl;
  const int mp = -1;
  const double mf = 0.90;
  vec< triple<int64_t,int64_t,int> > als;
  SearchFastb2( reads_fastb_file, tmp_contigs_fastb, K, &als, 0, mp, mf, VERBOSE );
    
  log << Date( ) << ": load reads (lengths will be needed)" << endl;
  vecbvec reads( reads_fastb_file );

  log << Date( ) << ": fill aligns and index" << endl;
  aligns.clear( );
  aligns.reserve( als.size( ) );
  for ( size_t i = 0; i < als.size( ); i++ ) {
    int id = als[i].first;
    int t = als[i].second;
    int pos2 = ( als[i].third >= 0 ? als[i].third : - als[i].third - 1 );
    int Pos2 = pos2 + reads[id].isize( );
    bool rc = als[i].third >= 0;
    aligns.push( pos2, Pos2, t, target[t].size( ), rc  );
  }

  index.clear( );
  index.resize( reads.size( ), -2 );
  for ( size_t i = 0; i < als.size( ); i++ ) {
    int id = als[i].first;
    if ( index[id] == -2 ) index[id] = i;
    else index[id] = -1;
  }
  
}

/**
 * ReportScaffoldsN50
 */
void ReportScaffoldsN50( const vec<superb> &supers, ostream &out )
{
  vec<int> cg_lens;
  vec<int> s_lens;
  vec<int> gaps_unsigned;
  vec<int> gaps_signed;
  vec<int> devs;
  int n_contigs = 0;
  for (uint ii=0; ii<supers.size( ); ii++)
    n_contigs += supers[ii].Ntigs( );
  cg_lens.reserve( n_contigs );
  s_lens.reserve( supers.size( ) );
  gaps_unsigned.reserve( n_contigs );
  gaps_signed.reserve( n_contigs );
  devs.reserve( n_contigs );
  for (uint ii=0; ii<supers.size( ); ii++) {
    for (int cgpos=0; cgpos<supers[ii].Ntigs( ); cgpos++) {
      cg_lens.push_back( supers[ii].Len( cgpos ) );
      if ( cgpos < supers[ii].Ntigs( )-1 ) {
	gaps_unsigned.push_back( Max( 0, supers[ii].Gap( cgpos ) ) );
	gaps_signed.push_back( supers[ii].Gap( cgpos ) );
	devs.push_back( supers[ii].Dev( cgpos ) );
      }
    }
    s_lens.push_back( supers[ii].ReducedLength( ) );
  }
  
  if ( gaps_unsigned.size( ) < 1 )
    out << "\nNo gaps in this scaffolds - nothing to report.\n";
  else {
    sort( devs.begin( ), devs.end( ) );
    out << "\nsmallest gap:     " << Min( gaps_signed ) << "\n"
	<< "largest gap:      " << Max( gaps_signed ) << "\n"
	<< "gap devs' median: " << Median( devs ) << "\n";
  }
  out << "\n";

  String name = "contig_size";
  PrintBasicNStats( name, cg_lens, out );
  out << "\n";

  name = "scaffold_size";
  PrintBasicNStats( name, s_lens, out );
  out << "\n";
  
}

/**
 * ReportScaffoldsBrief
 */
void ReportScaffoldsBrief( const vec<superb> &supers,
			   const int min_links,
			   const int step,
			   ostream &out )
{
  if ( supers.size( ) < 1 ) {
    out << "no supers found." << endl;
    return;
  }

  int n_contigs = 0;
  for (uint ii=0; ii<supers.size( ); ii++)
    n_contigs += supers[ii].Ntigs( );

  vec<int> slens;
  vec<int> gaps;
  vec<int> devs;
  slens.reserve( supers.size( ) );
  gaps.reserve( n_contigs );
  devs.reserve( n_contigs );
  for (uint ii=0; ii<supers.size( ); ii++) {
    slens.push_back( supers[ii].TrueLength( ) );
    for (int cgpos=0; cgpos<supers[ii].Ntigs( ); cgpos++) {
      if ( cgpos < supers[ii].Ntigs( )-1 ) {
	gaps.push_back( supers[ii].Gap( cgpos ) );
	devs.push_back( supers[ii].Dev( cgpos ) );
      }
    }
  }
  
  String str_gapdev = "na";
  String str_min = "na";
  String str_max = "na";
  if ( gaps.size( ) > 0 ) {
    str_gapdev = ToString( Median( devs ) );
    str_min = ToString( Min( gaps ) );
    str_max = ToString( Max( gaps ) );
  }
  
  out << "ml: " << min_links << "." << step << "   "
      << "ntigs: " << n_contigs << "   "
      << "nsupers: " << supers.size( ) << " (N50: " << N50( slens ) << ")   "
      << "gaps_range: [" << str_min << ", " << str_max << "]   "
      << "median_gap_dev: " << str_gapdev
      << endl;
}

/**
 * SaveInterimScaffolds
 */
void SaveInterimScaffolds( const String &data_dir,
			   const String &out_dir,
			   const PairsManager &pairs,
			   const vec<fastavector> &contigs,
			   const vec<superb> &supers,
			   const vec<alignlet> *aligns,
			   const vec<int> *index )
{
  String lookup_file = data_dir + "/genome.lookup";
  String head_out = out_dir + "/regap";
  String aligns_file = out_dir + "/aligns.qltoutlet";
  String index_file = out_dir + "/aligns.index";
  String contig_graph_file = head_out + ".contig.graph";
  String super_graph_file = head_out + ".super.graph";
  String log_file = head_out + ".log";
  
  Mkpath( out_dir );
  
  // The contig graph (singleton-scaffold graph).
  digraphE<sepdev> unused;
  digraphE<CLinkBundle> graph;
  vec<superb> sing;
  int n_contigs = contigs.size( );
  sing.resize( n_contigs );
  for (int tig=0; tig<n_contigs; tig++) {
    sing[tig].SetNtigs( 1 );
    sing[tig].SetTig( 0, tig );
    sing[tig].SetLen( 0, contigs[tig].size( ) );
  }

  BuildScaffoldGraph( pairs, sing, *aligns, *index, unused, &graph );
  int fd = OpenForWrite( contig_graph_file );
  BinaryWrite( fd, graph );

  // The scaffold graph.
  digraphE<CLinkBundle> supergraph;
  BuildScaffoldGraph( pairs, supers, *aligns, *index, unused, &supergraph );
  fd = OpenForWrite( super_graph_file );
  BinaryWrite( fd, supergraph );	
  
  // Save the usual scaffold assembly files.
  SaveScaffoldAssembly( head_out, supers, contigs );

  // Save aligns and index (optional).
  if ( aligns && index ) {
    BinaryWrite3( aligns_file, *aligns );
    BinaryWrite3( index_file, *index );
  }  

  // Run eval scaffolds (optional).
  if ( IsRegularFile( lookup_file ) ) {
    String theCommand 
      = "EvalScaffolds LOOKUP=" + lookup_file
      + " SCAFFOLDS=" + head_out;
    RunCommand( theCommand );
  }

}

/**
 * UpdateIndexFile
 */
size_t UpdateIndexFile( const vec<int> &orig_index, vec<int> &filt_index )
{
  ForceAssertEq( orig_index.size( ), filt_index.size( ) );
  size_t n_events = 0;
  for (size_t ii=0; ii<orig_index.size( ); ii++) {
    if ( orig_index[ii] < 0 && filt_index[ii] >=0 ) {
      filt_index[ii] = orig_index[ii];
      n_events++;
    }
  }
  return n_events;
}			

/**
 * CheckScaffoldsIntegrity
 */
bool CheckScaffoldsIntegrity( const size_t n_contigs,
			      const vec<superb> &supers )
{
  vec<bool> placed( n_contigs, false );
  for (size_t ii=0; ii<supers.size( ); ii++)
    for (int jj=0; jj<supers[ii].Ntigs( ); jj++)
      placed[ supers[ii].Tig( jj ) ] = true;

  int n_unplaced = 0;
  for (size_t ii=0; ii<n_contigs; ii++)
    if ( ! placed[ii] ) n_unplaced++;
  
  return ( n_unplaced < 1 );
}			

