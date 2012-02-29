///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include <omp.h>

#include "CoreTools.h"
#include "FetchReads.h"
#include "FetchReadsAmb.h"
#include "lookup/LookAlign.h"
#include "paths/reporting/ReftigUtils.h"
#include "util/SearchFastb2Core.h"

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS
// MakeDepend: dependency QueryLookupTable

Bool cmp_target( const look_align& la1, const look_align& la2 )
{    if ( la1.target_id < la2.target_id ) return True;
     if ( la1.target_id > la2.target_id ) return False;
     if ( la1.pos2( ) < la2.pos2( ) ) return True;
     if ( la1.pos2( ) > la2.pos2( ) ) return False;
     if ( la1.query_id < la2.query_id ) return True;
     return False;    }

void GetAlignsFast( const int K, const String& fastb_file, const String& lookup_file, 
		    const String& aligns_file, vec<look_align>& aligns,
		    const Bool USE_CACHE, const String& tmp_dir )
{    if ( !USE_CACHE || !IsRegularFile(aligns_file) ) 
     {    vec< triple<int64_t,int64_t,int> > ALIGNS;
          String gfastb = lookup_file.Before( ".lookup" ) + ".fastb";
          if ( !IsRegularFile(fastb_file) )
          {    vecbasevector seqs;
               FetchReads( seqs, 0, fastb_file.Before( ".fastb" ) + ".fasta" );
               seqs.WriteAll(fastb_file);    }
          vecbasevector seqs(fastb_file), genome(gfastb);
	  size_t nseqs = 0;
	  for (size_t ii=0; ii<seqs.size( ); ii++)
	    if ( (int)seqs[ii].size( ) >= K )
	      nseqs++;
	  cout << Date( ) << " (GAF): calling SearchFastb2 on "
	       << nseqs <<  " sequences" << endl;
          SearchFastb2( fastb_file, gfastb, K, &ALIGNS, 0, -1, 0.5, False );
	  {
	    for ( size_t i = 0; i < ALIGNS.size( ); i++ )
	      seqs[ ALIGNS[i].first ].resize(0);
	    size_t nseqs2 = 0;
	    for (size_t ii=0; ii<seqs.size( ); ii++)
	      if ( (int)seqs[ii].size( ) >= K )
		nseqs2++;
	    seqs.WriteAll( fastb_file + ".modified" );
	    cout << Date( ) << " (GAF): running QueryLookupTable on "
		 << nseqs2 << " queries" << endl;
	  }	  
          int nproc = omp_get_max_threads( ), N = seqs.size( );
          for ( int j = 0; j < nproc; j++ )
          {    Ofstream( idsout, aligns_file + ".ids." + ToString(j) );
               for ( int i = 0; i < N; i++ )
                    if ( i % nproc == j ) idsout << i << "\n";    }
          #pragma omp parallel for
          for ( int j = 0; j < nproc; j++ )
          {    SystemSucceed( "QueryLookupTable K=12 MM=12 MC=0.15 SEQS=" 
                    + fastb_file + ".modified " + ARG(SEQS_IS_FASTB, True) 
                    + ARG(L, lookup_file) + ARG(PARSEABLE, True) 
                    + ARG(SMITH_WAT, True) + ARG(TMP_DIR,tmp_dir)
                    + ARG(SEQS_TO_PROCESS, "@" + aligns_file + ".ids." + ToString(j))
                    + " > " + aligns_file + "." + ToString(j) );    }
          {    Ofstream( out, aligns_file );
               for ( int j = 0; j < nproc; j++ )
               {    vec<look_align> A;
                    LoadLookAligns( aligns_file + "." + ToString(j), A );
                    for ( int i = 0; i < A.isize( ); i++ )
                         A[i].PrintParseable(out);
                    Remove( aligns_file + ".ids." + ToString(j) );
                    Remove( aligns_file + "." + ToString(j) );    }    }
          Remove( fastb_file + ".modified" );
          seqs.ReadAll(fastb_file);
          OfstreamMode( out, aligns_file, ios::app );
          for ( size_t i = 0; i < ALIGNS.size( ); i++ )
          {    int id = ALIGNS[i].first, g = ALIGNS[i].second; 
               int pos = ALIGNS[i].third;
               int xpos = ( pos >= 0 ? pos : -pos-1 );
               avector<int> gaps(1), lengths(1);
               gaps(0) = 0;
               lengths(0) = seqs[id].size( );
               align a( 0, xpos, gaps, lengths );
               look_align la( id, g, seqs[id].size( ), genome[g].size( ), pos < 0,
                    a, 0, 0, 0 );
               la.PrintParseable(out);    }    }
     cout << Date( ) << " (GAF): loading aligns" << endl;
     LoadLookAligns( aligns_file, aligns );    
     sort( aligns.begin( ), aligns.end( ), cmp_target );    }

template<class HYPER_T> void HyperToReftigsCore( const int K, const HYPER_T& h,
     const vec<look_align>& aligns, vec< pair<int,ho_interval> >& reftigs,
     digraph *alignsG )
{
     reftigs.clear( );

     // Define heuristic constants.

     const int delta = 50;
     const int delta2 = 2000;

     // Form the graph G whose vertices are alignments and whose adjacencies 
     // correspond to 'adjacent' alignments.

     vec<int> to_left, to_right;
     h.ToLeft(to_left), h.ToRight(to_right);
     int nA = aligns.size( );
     vec< vec<int> > from(nA), to(nA);
     vec<int> places( h.EdgeObjectCount( ), 0 );
     for ( int i = 0; i < aligns.isize( ); i++ )
          ++places[ aligns[i].query_id ];
     for ( int i1 = 0; i1 < aligns.isize( ); i1++ )
     {    for ( int i2 = i1 + 1; i2 < aligns.isize( ); i2++ )
          {    const look_align &la1 = aligns[i1], &la2 = aligns[i2];
               if ( la2.target_id != la1.target_id ) break;
               if ( la2.pos2( ) >= la1.Pos2( ) + delta2 ) break;
               if ( la1.Fw1( ) != la2.Fw1( ) ) continue;
               int id1 = la1.query_id, id2 = la2.query_id;
               int discrep = Abs( la1.Pos2( ) - la2.pos2( ) - (K-1) );
               Bool OK = False;
               if ( discrep <= delta ) OK = True;
               if ( places[id1] == 1 && places[id2] == 1 && discrep <= delta2 )
                    OK = True;
               if ( !OK ) continue;
               if ( la1.Fw1( ) && to_right[id1] != to_left[id2] ) continue;
               if ( la1.Rc1( ) && to_right[id2] != to_left[id1] ) continue;
               if ( la1.Pos2( ) >= la2.Pos2( ) ) continue;
               from[i1].push_back(i2), to[i2].push_back(i1);    }    }
     for ( int i = 0; i < nA; i++ )
          Sort( to[i] );
     digraph G( from, to );
     if ( alignsG ) *alignsG = G;

     // By construction, G is acylic.  Therefore we can understand its 'extent'
     // on the reference by finding the sources and sinks of its components.

     vec< vec<int> > comp;
     G.Components(comp);
     for ( int i = 0; i < comp.isize( ); i++ )
     {    vec<int> sources, sinks;
          G.SubgraphSources( comp[i], sources );
          G.SubgraphSinks( comp[i], sinks );
          int left = 1000000000, right = 0;
          for ( int j = 0; j < sources.isize( ); j++ )
               left = Min( left, aligns[ sources[j] ].pos2( ) );
          for ( int j = 0; j < sinks.isize( ); j++ )
               right = Max( right, aligns[ sinks[j] ].Pos2( ) );
          reftigs.push_back( make_pair( aligns[ comp[i][0] ].target_id,
               ho_interval( left, right ) ) );    }

     // Clean up answer.

     UniqueSort(reftigs);
     vec<Bool> to_delete( reftigs.size( ), False );
     for ( int i1 = 0; i1 < reftigs.isize( ); i1++ )
     {    for ( int i2 = i1 + 1; i2 < reftigs.isize( ); i2++ )
          {    if ( reftigs[i2].first != reftigs[i1].first ) break;
               if ( reftigs[i2].second.Start( ) > reftigs[i1].second.Stop( ) ) break;
               if ( Subset( reftigs[i1].second, reftigs[i2].second ) ) 
                    to_delete[i1] = True;
               if ( Subset( reftigs[i2].second, reftigs[i1].second ) )
                    to_delete[i2] = True;    }    }
     EraseIf( reftigs, to_delete );    }

template<class HYPER_T>
void GenerateDot( const int ORIGIN,
		  const String &dot_base,
		  const HYPER_T &hyper,
		  const digraph &agraph,
		  const vec<look_align> &aligns,
		  const vec< pair<int,ho_interval> >& reftigs,
		  int VERBOSITY )
{
  ForceAssert( is_sorted( reftigs.begin( ), reftigs.end( ) ) );

  String dot_file = dot_base + ".dot";

  int n_edges = hyper.EdgeObjectCount( );
  vec<double> lengths( n_edges );
  for (int ii=0; ii<n_edges; ii++)
    lengths[ii] = hyper.EdgeObject( ii ).MidLength( );

  // Map edge_id to align_ids.
  vec< vec<int> > to_align( n_edges );
  for (size_t ii=0; ii<aligns.size( ); ii++) {
    int edge_id = aligns[ii].query_id;
    to_align[edge_id].push_back( int( ii ) );
  }

  // Map align_id to reftig_id.  
  vec<int> to_reftig( aligns.size( ), -1 );
  {
    vec< vec<int> > comp;
    agraph.Components( comp );
    
    for (int comp_id=0; comp_id<(int)comp.size( ); comp_id++) {
      int v0 = comp[comp_id][0];
      const look_align &al0 = aligns[v0];

      // Find reftig id for this components.
      int reftig_id = -1;
      for (int rid=0; rid<(int)reftigs.size( ); rid++) {
	if ( al0.target_id != reftigs[rid].first ) continue;
	if ( al0.pos2( ) < reftigs[rid].second.Start( ) ) continue;
	if ( al0.Pos2( ) > reftigs[rid].second.Stop( ) ) continue;
	reftig_id = rid;
	break;
      }
      
      // Tag all aligns in component as belonging to this reftig id.
      for (int ii=0; ii<(int)comp[comp_id].size( ); ii++)
	to_reftig[ comp[comp_id][ii] ] = reftig_id;
    }
  }
  
  // Labels for edges in the dot file.
  vec<String> eLabels( n_edges, "" );
  if ( VERBOSITY > 0 ) {
    for (int edge_id=0; edge_id<n_edges; edge_id++) {
      const vec<int> &align_ids = to_align[edge_id];
      if ( align_ids.size( ) < 1 ) continue;
      for (size_t ii=0; ii<align_ids.size( ); ii++) {
	int al_id = align_ids[ii];
	int r_id = to_reftig[al_id];
	String str_r_id = ToString( r_id );
	if ( r_id < 0 ) continue;
	const look_align &al = aligns[al_id];
	eLabels[edge_id] += " [r" + str_r_id + ( al.Rc1( ) ? "-" : "+" );
	if ( VERBOSITY > 1 ) {
	  String strId = ToString( al.target_id );
	  String strBeg = ToStringAddCommas( ORIGIN + al.pos2( ) );
	  String strEnd = ToStringAddCommas( ORIGIN + al.Pos2( ) );
	  eLabels[edge_id] += ": " + strId + "." + strBeg + "-" + strEnd;
	}
	eLabels[edge_id] += "]";
      }
    }
  }
  
  vec< vec<int> > compE;
  hyper.ComponentsE( compE );
  int n_components = compE.size( );

  // Map component_id to reftig_ids.
  vec< vec<int> > comp_to_reftigs( n_components );
  for (int comp_id=0; comp_id<n_components; comp_id++) {
    vec<int> lmap;
    for (size_t epos=0; epos<compE[comp_id].size( ); epos++) {
      int edge_id = compE[comp_id][epos];
      const vec<int> &align_ids = to_align[edge_id];
      for (size_t ii=0; ii<align_ids.size( ); ii++) {
	int al_id = align_ids[ii];
	int r_id = to_reftig[al_id];
	if ( r_id < 0 ) continue;
	lmap.push_back( r_id );
      }
    }
    sort( lmap.begin( ), lmap.end( ) );
    lmap.erase( unique( lmap.begin( ), lmap.end( ) ), lmap.end( ) );

    comp_to_reftigs[comp_id] = lmap;
  }

  // Contig labels (ie the reftig ranges).
  vec<String> compRange( n_components, "" );
  vec<int> low_rid( n_components, -1 );
  for (int comp_id=0; comp_id<n_components; comp_id++) {
    const vec<int> &local_reftigs= comp_to_reftigs[comp_id];
    if ( local_reftigs.size( ) < 1 ) continue;
    low_rid[comp_id] = Min( local_reftigs );
    for (size_t ii=0; ii<local_reftigs.size( ); ii++) {
      int reftig_id = local_reftigs[ii];
      if ( ii > 0 ) compRange[comp_id] += "  ";
      compRange[comp_id]
	+= "r" + ToString( reftig_id )
	+ " [" + ToStringAddCommas( reftigs[reftig_id].first )
	+ "." + ToStringAddCommas( ORIGIN + reftigs[reftig_id].second.Start( ) )
	+ "-" + ToStringAddCommas( ORIGIN + reftigs[reftig_id].second.Stop( ) )
	+ "]";
    }
  }
  for (int ii=0; ii<n_components; ii++)
    if ( compRange[ii] == "" )
      compRange[ii] = "[unmapped]";

  // Sort components by reftig placement on reference.
  vec<int> comp_ids ( n_components );
  for (int ii=0; ii<n_components; ii++)
    comp_ids[ii] = ii;
  SortSync( low_rid, comp_ids );

  // Run PrettyDOT.
  bool alpha = VERBOSITY == 2 ? False: True;
  ofstream out( dot_file.c_str( ) );
  hyper.PrettyDOT( out, lengths, True, False, True, &comp_ids,
		   alpha, &eLabels, &compRange );
  out.close( );
}

void PrintReftigs ( ostream &out,
		    const int K,
		    const int ORIGIN,
		    const vec< pair<int,ho_interval> > &reftigs,
		    const vecbitvector *amb,
		    const int min_len )
{
  vec< vec<ho_interval> > gaps;
  if ( amb ) {
    gaps.resize( amb->size( ) );
    for ( size_t i = 0; i < amb->size( ); i++ ) {
      for ( uint j = 0; j < (*amb)[i].size( ); j++ ) {
	if ( !(*amb)[i][j] ) continue;
	uint k;
	for ( k = j + 1; k < amb[i].size( ); k++ )
	  if ( !(*amb)[i][k] ) break;
	gaps[i].push( j, k );
	j = k;
      }
    }
  }	  
  for ( int i = 0; i < reftigs.isize( ); i++ ) {
    int t = reftigs[i].first;
    if ( amb && i > 0 && reftigs[i-1].first == t ) {
      vec<ho_interval> g;
      int low = Min( reftigs[i-1].second.Stop( ), 
		     reftigs[i].second.Start( ) ) - K;
      int high = Max( reftigs[i-1].second.Stop( ), 
		      reftigs[i].second.Start( ) ) + K;
      int gsum = 0;
      for ( int j = 0; j < gaps[t].isize( ); j++ ) {
	if ( Meets( gaps[t][j], ho_interval(low, high) ) ) {
	  g.push_back( gaps[t][j] );
	  gsum += gaps[t][j].Length( );
	}
      }
      if ( gsum + 4*K >= (high-low)/4 ) {
	for ( int u = 0; u < g.isize( ); u++ ) {
	  if ( g[u].Length( ) < min_len ) continue;
	  out << t << "." 
	      << ToStringAddCommas( g[u].Start( ) + ORIGIN ) 
	      << "-"
	      << ToStringAddCommas( g[u].Stop( ) + ORIGIN ) 
	      << " (N)" << endl;
	}
      }
    }
    if ( reftigs[i].second.Length( ) >= min_len ) {
      out << t << "." 
	  << ToStringAddCommas( reftigs[i].second.Start( ) + ORIGIN ) 
	  << "-"
	  << ToStringAddCommas( reftigs[i].second.Stop( ) + ORIGIN ) 
	  << endl;
    }
  }
  out << "\n";
  
}

/**
 * Instantiate templates
 */
template void HyperToReftigsCore( const int K, const HyperKmerPath& h,
				  const vec<look_align>& aligns,
				  vec< pair<int,ho_interval> >& reftigs,
				  digraph *alignsG );

template void HyperToReftigsCore( const int K, const HyperBasevector& h,
				  const vec<look_align>& aligns,
				  vec< pair<int,ho_interval> >& reftigs,
				  digraph *alignsG );

template void GenerateDot( const int ORIGIN, const String &dot_base,
			   const HyperKmerPath &hyper, const digraph &agraph,
			   const vec<look_align> &aligns,
			   const vec< pair<int,ho_interval> >& reftigs,
			   int VERBOSITY );

template void GenerateDot( const int ORIGIN, const String &dot_base,
			   const HyperBasevector &hyper, const digraph &agraph,
			   const vec<look_align> &aligns,
			   const vec< pair<int,ho_interval> >& reftigs,
			   int VERBOSITY );
