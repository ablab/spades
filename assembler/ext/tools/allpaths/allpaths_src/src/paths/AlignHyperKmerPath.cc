///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Alignment.h"
#include "Basevector.h"
#include "CommonSemanticTypes.h"
#include "CoreTools.h"
#include "Equiv.h"
#include "graph/DigraphTemplate.h"
#include "lookup/LookAlign.h"
#include "lookup/PerfectLookup.h"
#include "lookup/QueryLookupTableCore.h"
#include "math/Functions.h"
#include "math/HoInterval.h"
#include "pairwise_aligners/AffineRefiner.h"
#include "paths/AlignHyperKmerPath.h"
#include "paths/EvalUtils.h"
#include "paths/HyperKmerPath.h"
#include "paths/KmerBaseBroker.h"
#include "paths/KmerPath.h"

// Function: AlignHyperKmerPath
//
// AlignHyperKmerPath takes a HyperKmerPath h, whose KmerPath edges are assumed
// not to have any gaps, and aligns each of them to a reference genome (using files 
// GENOME.fastb and GENOME.lookup), then removes some edges alignments which appear 
// to be wrong.  The result is returned as a vec<look_align>, in which the query 
// ids are the edge ids in h.  An index (mapping each edge to its alignment(s)
// to the reference) is also generated.  A work directory must 
// be provided for intermediate calculations.
//
// We use two programs to do this, QueryLookupTable and FindAllPerfectAlignments.  
// The first finds alignments, whether perfect or not, but will miss multiple 
// alignments of the same sequence to overlapping places on the reference.  The 
// second is guaranteed to find all perfect alignments.
//
// Output parameters:
//
//    aligns - each element gives the alignment of one edge of the HyperKmerPath h
//             to the reference.
//    aligns_index - for each edge, all alignments of this edge to the reference
//       (the indices in 'aligns' of the alignments of this edge).
//

// There are no clocks in here.  The output is burdensome to general users.

void AlignHyperKmerPath( const HyperKmerPath& h, const KmerBaseBroker* kbb,
     const String& GENOME, const String& tmp_dir, vec<look_align>& aligns,
     vec< vec<int> >& aligns_index, bool filter )
{    
     // Get edges and genome.
               
     vecbasevector edges, genome( GENOME + ".fastb" );
     int K = h.K( );
     int nbases = 0, nedges = h.EdgeObjectCount( );
     for ( int i = 0; i < nedges; i++ )
          nbases += h.EdgeLength(i) + K - 1;
     edges.Reserve( nbases/16 + nedges, nedges );
     for ( int i = 0; i < nedges; i++ )
     {
       SuperBaseVector s = kbb->ToSequence( h.EdgeObject(i) );
          ForceAssertEq( s.size( ), 1 );
          edges.push_back( s.Seq(0) );    }

     // Find perfect alignments, build index.

     ForceAssertLe( 12, K );
     PerfectLookup( 12, edges, GENOME + ".lookup", aligns, FW_OR_RC );
     aligns_index.clear_and_resize(nedges);
     for ( int i = 0; i < aligns.isize( ); i++ )
          aligns_index[ aligns[i].query_id ].push_back(i);

     // Find edges that don't have a perfect alignment, realign using 
     // QueryLookupTable.

     vec<int> unaligned;
     vecbasevector unaligned_edges;
     for ( int i = 0; i < nedges; i++ )
     {    if ( aligns_index[i].empty( ) )
          {    unaligned.push_back(i);
               unaligned_edges.push_back_reserve( edges[i] );    }    }
     if ( unaligned.nonempty( ) )
     {    temp_file un_fastb( tmp_dir + "/unedges.fastb.XXXXXXX" );
          unaligned_edges.WriteAll(un_fastb);
          temp_file qlt( tmp_dir + "/QueryLookupTable.out.XXXXXXX" );
          QueryLookupTableCore( "K=12 L=" + GENOME + ".lookup " + " MM="
               + ToString(K-1) + " SEQS=" + un_fastb + " SMITH_WAT=True"
               + " SEQS_IS_FASTB=True PARSEABLE=True NH=True OUTFILE=" + qlt );
          static vec<look_align> aligns2;
          static vec< vec<int> > aligns_index2;
          LoadLookAligns( qlt, aligns2, aligns_index2, unaligned.size( ) );
          for ( int i = 0; i < aligns2.isize( ); i++ )
          {    look_align& la = aligns2[i];
               basevector e = unaligned_edges[ la.query_id ];
               if ( la.Rc1( ) ) e.ReverseComplement( );
               AffineRefiner AR;
               AR.RefineAlign( la.a, e, genome[la.target_id] );    }
          for ( int i = 0; i < aligns2.isize( ); i++ )
               aligns2[i].query_id = unaligned[ aligns2[i].query_id ];
          for ( int i = 0; i < aligns_index2.isize( ); i++ )
          {    for ( int j = 0; j < aligns_index2[i].isize( ); j++ )
                    aligns_index2[i][j] += aligns.size( );
               aligns_index[ unaligned[i] ] = aligns_index2[i];    }
          aligns.append(aligns2);    }

     // Try to remove incorrect alignments, if asked.

     if ( filter ) {
       vec<Bool> to_remove( aligns.size( ), False );
       for ( int pass = 1; pass <= 2; pass++ )
       {    for ( int v = 0; v < h.N( ); v++ )
            {    const vec<int>& S1 = ( pass == 1 ? h.To(v) : h.From(v) );
                 if ( S1.size( ) != 1 ) continue;
                 int u1, u2;
                 if ( pass == 1 ) u1 = h.EdgeObjectIndexByIndexTo( v, 0 );
                 else u1 = h.EdgeObjectIndexByIndexFrom( v, 0 );
                 if ( aligns_index[u1].size( ) != 1 ) continue;
                 const look_align& la1 = aligns[ aligns_index[u1][0] ];
                 if ( la1.rc1 || !la1.FullLength( ) ) continue;
                 const vec<int>& S2 = ( pass == 1 ? h.From(v) : h.To(v) );
                 for ( int i = 0; i < S2.isize( ); i++ )
                 {    if ( pass == 1 ) u2 = h.EdgeObjectIndexByIndexFrom( v, i );
                      else u2 = h.EdgeObjectIndexByIndexTo( v, i );
                      if ( aligns_index[u2].size( ) <= 1 ) continue;
                      static vec<int> matches;
                      matches.clear( );
                      int min_errs = 1000000;
                      for ( int j = 0; j < aligns_index[u2].isize( ); j++ )
                      {    const look_align& la2 
                                = aligns[ aligns_index[u2][j] ];
                           if ( la2.rc1 || !la2.FullLength( ) ) continue;
                           min_errs = Min( min_errs, la2.Errors( ) );
                           if ( pass == 1 && la1.a.Pos2( ) - la2.a.pos2( ) == K - 1 )
                                matches.push_back(j);
                           if ( pass == 2 && la2.a.Pos2( ) - la1.a.pos2( ) == K - 1 )
                                matches.push_back(j);    }
                      if ( matches.size( ) != 1 ) continue;
                      int keep = aligns_index[u2][ matches[0] ];
                      const look_align& la2 = aligns[keep];
                      if ( la2.Errors( ) > min_errs ) continue;
                      for ( int j = 0; j < aligns_index[u2].isize( ); j++ )
                           if ( aligns_index[u2][j] != keep )
                                to_remove[ aligns_index[u2][j] ] = True;
                      aligns_index[u2].resize(1);
                      aligns_index[u2][0] = keep;    }    }    }
       EraseIf( aligns, to_remove );
     }

     // Rebuild index.

     aligns_index.clear( );
     aligns_index.resize(nedges);
     for ( int i = 0; i < aligns.isize( ); i++ )
          aligns_index[ aligns[i].query_id ].push_back(i);    }

// Function: ReorderToFollowReference
//
// Given edge alignments, try to reorder the <components>
// of a HyperKmerPath so that they follow the reference.  Flip them (both components 
// and aligns) if necessary.

void ReorderToFollowReference( HyperKmerPath& h, vec<look_align>& aligns,
     const vec< vec<int> >& aligns_index )
{    equiv_rel e;
     h.ComponentRelation(e);
     vec<int> reps;
     e.OrbitRepsAlt(reps);

     // Local var: reps_pos
     // For each <component>, (position-on-reference, component-id).
     // ("position" of a component is, here, the leftmost position of its
     // maximal-length edge).
     // Once we get this for each component, we can sort by position on reference.
     // Note that since the reference may consist of multiple unconnected genome parts (e.g. chromosomes),
     // position on reference is a pair (genome-part-id, position-on-that-genome-part).
     vec< pair< pair<genome_part_id_t,genome_part_pos_t>, int > > reps_pos;
     for ( int i = 0; i < reps.isize( ); i++ )
     {
       // Get a list of edges in this component
          static vec<int> o;
          static vec<int> oe;
          e.Orbit( reps[i], o );
          oe.clear( );
          for ( int j = 0; j < o.isize( ); j++ )
          {    int v = o[j];
               for ( int u = 0; u < h.From(v).isize( ); u++ )
                    oe.push_back( h.EdgeObjectIndexByIndexFrom( v, u ) );    }
          UniqueSort(oe);

	  // Get a list of the longest (maximal-length) edges of this component.
	  // We'll align this component to the reference by aligning its longest edges.
	  // There may be several longest edges, and each may have several alignments
	  // to the reference, so there may be several candidate alignments of this
	  // component to the reference.
          nkmers_t maxkmers = 0;
          vec<int> maxkmers_index;
          for ( int j = 0; j < oe.isize( ); j++ )
          {    nkmers_t x = h.EdgeObject( oe[j] ).KmerCount( );
               if ( x > maxkmers )
               {    maxkmers = x;
                    maxkmers_index.clear( );    }
               if ( x >= maxkmers ) maxkmers_index.push_back( oe[j] );    }

	  // Find the "leftmost" alignment of any maximal-length edge.  That will be the
	  // alignment of this component to the reference.  We'll order components
	  // in order of this alignment.
	  
          int infinity = 1000000000;
          genome_part_id_t min_tig = infinity;
	  genome_part_pos_t min_pos = infinity;
	  Bool min_rc = 0;
          for ( int u = 0; u < maxkmers_index.isize( ); u++ )
          {    int M = maxkmers_index[u];
               for ( int j = 0; j < aligns_index[M].isize( ); j++ )
               {    const look_align& la = aligns[ aligns_index[M][j] ];
                    if ( la.target_id > min_tig ) continue;
                    if ( la.target_id < min_tig )
                    {    min_tig = la.target_id;
                         min_pos = la.a.pos2( );
                         min_rc = la.rc1;    }
                    else 
                    {    if ( la.a.pos2( ) < min_pos ) 
                              min_pos = la.a.pos2( );
                              min_rc = la.rc1;    }    }    }

	  // Flip the component (and its alignment to the reference) if necessary
	  // (if the leftmost alignment of a maximal-length edge of the component to the reference
	  // is an rc alignment).
          if ( min_rc ) {
            h.ReverseComponent( reps[i] );
            for ( int edgeIdx = 0; edgeIdx < oe.isize(); ++edgeIdx ) {
              int edgeId = oe[edgeIdx];
              for ( int j = 0; j < aligns_index[edgeId].isize(); j++ ) {
                look_align& la = aligns[ aligns_index[edgeId][j] ];
                la.rc1 = ! la.rc1; 
              }
            }
          }  // if ( min_rc ) ...
          reps_pos.push_back( make_pair( make_pair( min_tig, min_pos ), i ) );
     }  // for each component, get its leftmost alignment to the reference.
     Sort(reps_pos);
     vec<int> new_order;
     for ( int i = 0; i < reps.isize( ); i++ )
          new_order.push_back( reps_pos[i].second );
     h.ReorderComponents(new_order);
}  // ReorderToFollowReference()

struct cmp_target: public binary_function<const look_align&, const look_align&, Bool> {
public:
  Bool operator()( const look_align& la1, const look_align& la2 ) const
  {
    if ( la1.target_id < la2.target_id ) return True;
    if ( la1.target_id > la2.target_id ) return False;
    if ( la1.rc1 < la2.rc1 ) return True;
    if ( la1.rc1 > la2.rc1 ) return False;
    if ( la1.a.pos2( ) < la2.a.pos2( ) ) return True;
    if ( la1.a.pos2( ) > la2.a.pos2( ) ) return False;
    return False;
  }
};

// Function: PrintAlignedHyperKmerPath
//
// Given edge alignments, print a HyperKmerPath, by
// components, with the edge alignments displayed.

void PrintAlignedHyperKmerPath( ostream& out, const HyperKmerPath& h,      
     const KmerBaseBroker* kbb, const vecbasevector& genome, 
     const vec<look_align>& aligns, const vec< vec<int> >& aligns_index,
     Bool print_component_id_line, const vec<TrustedPath>* trusted_pathsp,
     const Bool brief, const Bool diploid )
{
     // Set up map from edges to vertices.

     int nedges = h.EdgeObjectCount( );
     vec<int> to_left_vertex, to_right_vertex;
     h.ToLeft(to_left_vertex), h.ToRight(to_right_vertex);

     // Map trusted paths to components.

     equiv_rel e;
     h.ComponentRelation(e);
     vec<int> reps;
     e.OrbitRepsAlt(reps);

     // Local var: tpc
     // For each component, the <trusted paths> through that component
     // (as indices into *trusted_pathsp).  Trusted paths through components
     // are found by FilterByReference().

     vec< vec<int> > tpc( reps.size( ) );
     if ( trusted_pathsp != 0 )
     {    const vec<TrustedPath>& trusted_paths = *trusted_pathsp;
          for ( int i = 0; i < trusted_paths.isize( ); i++ )
          {    const TrustedPath& p = trusted_paths[i];
               int c = e.ClassId( to_left_vertex[ p.GetAlign(0).query_id ] );
               for ( int j = 0; j < reps.isize( ); j++ )
               {    if ( e.ClassId( reps[j] ) == c )
                    {    tpc[j].push_back(i);
                         break;    }    }    }    }

     // Define component remarks based on the trusted paths.

     vec<String> component_remarks( reps.size( ) );
     vec<Bool> component_remarks_only( reps.size( ), False );
     vec<int> to_left, to_right;
     h.ToLeft(to_left), h.ToRight(to_right);
     if ( trusted_pathsp != 0 )
     {    const vec<TrustedPath>& trusted_paths = *trusted_pathsp;
          for ( int r = 0; r < reps.isize( ); r++ )
          {    vec<int> o;
               e.Orbit( reps[r], o );
               int nsources = 0, nsinks = 0;
               for ( int i = 0; i < o.isize( ); i++ )
               {    if ( h.Source( o[i] ) ) ++nsources;
                    if ( h.Sink( o[i] ) ) ++nsinks;    }

	       // Local var: npaths
	       // The number of trusted paths through the current component.

               int npaths = tpc[r].size( );

	       // Local vars: Info about each trusted path through the current 
               // component.
	       //
	       //   errors - the number of alignment errors in this trusted path
	       //   haploid_errors - the number of alignment errors in edges that 
               //      align to only one strand of a diploid genome
	       //   proper - whether all the path's edges align fully to the 
               //      reference
	       //   source_start - whether the path starts at a <source node>
	       //   sink_stop - whether the path stops at a <sink node>

               vec<int> errors( npaths, 0 ), haploid_errors( npaths, 0 );
               vec<Bool> proper( npaths, True );
               vec<Bool> source_start( npaths, False ), sink_stop( npaths, False );
               for ( int u = 0; u < npaths; u++ )
               {    int i = tpc[r][u];
                    const TrustedPath& p = trusted_paths[i];
                    int naligns = p.GetNumAligns( );
                    int firstVertex = p.GetAlign(0).query_id;
                    int lastVertex = p.GetAlign(naligns-1).query_id;
                    if ( p.IsFw( ) )
                    {    if ( h.Source( to_left[ firstVertex ] ) )
                              source_start[u] = True;
                         if ( h.Sink( to_right[ lastVertex ] ) )
                              sink_stop[u] = True;    }
                    else
                    {    if ( h.Source( to_left[ lastVertex ] ) )
                              source_start[u] = True;
                         if ( h.Sink( to_right[ firstVertex ] ) )
                              sink_stop[u] = True;    }
                    for ( int j = 0; j < naligns; j++ )
                    {    const look_align& la = p.GetAlign(j);
                         errors[u] += la.Errors( );
                         int m = p.GetAlign(j).query_id;
                         if ( aligns_index[m].solo( ) && diploid )
                              haploid_errors[u] += la.Errors( );
                         if ( !la.FullLength( ) ) proper[u] = False;    }    }

	       // Count how many "good" trusted paths there are through this 
               // component.  Good means they (i.e. their edges) align to the 
               // reference with few errors.

               const int max_errors = 10;
               int num_good = 0, num_OK = 0;
               for ( int u = 0; u < npaths; u++ )
               {    if ( proper[u] && source_start[u] && sink_stop[u] ) 
                    {    if ( errors[u] <= max_errors ) ++num_OK;
                         if ( errors[u] == 0 ) ++num_good;    }    }

	       // For each trusted path through this component, construct a string
	       // showing the extent of the reference (i.e. the range of locations 
               // on the reference) covered by the path.  If the path is good or 
               // perfect (in matching the reference), add a remark 
               // highlighting this.

               vec<String> extent(npaths);
               for ( int u = 0; u < npaths; u++ )
               {    int i = tpc[r][u];
                    const TrustedPath& p = trusted_paths[i];
                    int naligns = p.GetNumAligns( );
                    genome_part_id_t tig = p.GetFinishedId();
                    genome_part_pos_t start = p.Begin();
                    genome_part_pos_t stop = p.End();
                    extent[u] = ToString(tig) + "." + ToString(start)
                         + "-" + ToString(stop);    }
               if ( nsources == 1 && nsinks == 1 && num_OK > 0 )
               {    component_remarks[r] = "\nThere is a unique source and a unique "
                         + String("sink and a ") 
                         + ( num_good > 0 ? "perfect" : "good" ) 
                         + " path from source to sink.\n";
                    int nmatches = 0;
                    for ( int u = 0; u < npaths; u++ )
                    {    if ( errors[u] > max_errors || !proper[u] ) continue;
                         if ( !source_start[u] || !sink_stop[u] ) continue;
                         ++nmatches;
                         component_remarks[r] += "Matches reference " + extent[u];
                         int i = tpc[r][u];
                         const TrustedPath& p = trusted_paths[i];
                         genome_part_pos_t start = p.Begin();
                         genome_part_pos_t stop = p.End();
                         component_remarks[r] += " (l=" + ToString(stop-start) + ")";
                         double nerrors = errors[u];
                         if (diploid) nerrors -= double(haploid_errors[u])/2.0;
                         ostringstream eout;
                         eout << setprecision(3) << nerrors;
                         if ( nerrors > 0 )
                         {    component_remarks[r] += ", has " + String(eout.str( ))
                                   + " base errors";    }
                         component_remarks[r] += ".\n";    }
                    if ( brief && !diploid && nmatches == 1 )
                         component_remarks_only[r] = True;
                    if ( brief && diploid && nmatches == 2 )
                         component_remarks_only[r] = True;
                    if ( !diploid && num_good >= 1 ||
                         diploid && num_good >= 2 || 
                         component_remarks_only[r] ) continue;    }

	       // Print out (to a string) the trusted paths through this component,
	       // with info on/features of each trusted path.
	       
               else
               {    ostringstream out;
                    out << "\ntrusted paths:\n";
                    for ( int u = 0; u < npaths; u++ )
                    {    int i = tpc[r][u];
                         const TrustedPath& p = trusted_paths[i];
                         if ( p.End( ) - p.Begin( ) < 1000 ) continue;
                         int n = p.GetNumAligns( );
                         vec<look_align> a;
                         for ( int j = 0; j < n; j++ )
                              a.push_back( p.GetAlign(j) );
                         out << "\n" << i << "[" << errors[u] 
                              << ( proper[u] ? "" : "*" ) << "]: " 
                              << p.GetFinishedId() << "." << p.Begin() << "-" 
                              << p.End();
                         if ( p.GetRc( ) ) 
                         {    a.ReverseMe( );
                              out << "(rc)";    }
                         out << "\n";
                         for ( int j = 0; j < n; j++ )
                         {    if ( j > 0 && j % 4 == 0 ) out << "\n";
                              const look_align& la = a[j];
                              int e = la.query_id;
                              int v = to_left_vertex[e], w = to_right_vertex[e];
                              if ( j == 0 )
                              {    out << v;
                                   if ( h.Source(v) && h.Sink(v) ) 
                                        out << "[source,sink]";
                                   else if ( h.Source(v) ) out << "[source]";
                                   else if ( h.Sink(v) ) out << "[sink]";    }
                              out << " -- " << h.EdgeLength(e) << "[" << la.Errors( ) 
                                   << ( la.FullLength( ) ? "" : "*" ) << "] --> " 
                                   << w;
                              if ( h.Source(w) && h.Sink(w) ) out << "[source,sink]";
                              else if ( h.Source(w) ) out << "[source]";
                              else if ( h.Sink(w) ) out << "[sink]";    }
                         out << "\n";    }    
                    component_remarks[r] = out.str( );    }    }
     }  // if we were given trusted paths

     // Get the sequence of each edge of the component
               
     vecbasevector edges;
     int nbases = 0;
     for ( int pass = 1; pass <= 2; pass++ )
     {    if ( pass == 2 ) edges.Reserve( nbases/16 + nedges, nedges );
          for ( int i = 0; i < nedges; i++ )
          {    const KmerPath& p = h.EdgeObject(i);
               SuperBaseVector s = kbb->ToSequence(p);
               ForceAssertEq( s.size( ), 1 );
               if ( pass == 1 ) nbases += s.Seq(0).size( );
               else edges.push_back( s.Seq(0) );    }    }

     // Define edge comments.

     vec<String> align_descrip(nedges);
     for ( int i = 0; i < nedges; i++ )
     {    ostringstream xout;
          //xout << endl;
          static vec<look_align> these_aligns;
          these_aligns.clear( );
          for ( int j = 0; j < aligns_index[i].isize( ); j++ )
               these_aligns.push_back( aligns[ aligns_index[i][j] ] );
          Sort( these_aligns,  cmp_target() );
          for ( int j = 0; j < these_aligns.isize( ); j++ )
          {    const look_align& la = these_aligns[j];
               int t;
               int delta = 0;
               for ( t = j + 1; t < these_aligns.isize( ); t++ )
               {    const look_align& lat = these_aligns[t];
                    if ( lat.target_id != la.target_id ) break;
                    if ( lat.rc1 != la.rc1 ) break;
                    if ( !la.FullLength( ) || !lat.FullLength( ) ) break;
                    if ( la.Errors( ) > 0 || lat.Errors( ) > 0 ) break;
                    if ( lat.Pos2( ) - lat.pos2( ) != la.Pos2( ) - la.pos2( ) ) 
                         break;
                    delta = these_aligns[t].pos2( ) - these_aligns[t-1].pos2( );
                    if ( delta > 500 ) break;
                    if ( t >= j + 2 && delta
                         != these_aligns[t-1].pos2( ) - these_aligns[t-2].pos2( ) )
                    {    break;    }    }
               xout << "  ";
               xout << ( la.rc1 ? "rc" : "fw" ) << " vs " << la.target_id << ", ";
               if ( t - j > 1 )
               {    xout << "perfect match to " << la.pos2( ) << "-" 
                         << la.Pos2( ) << " + " 
                         << these_aligns[j+1].pos2( ) - these_aligns[j].pos2( )
                         << "n, n = 0,...," << t - j - 1 << "\n";
                    j = t - 1;
                    continue;    }
               static qualvector query_qual;
               if ( la.Errors( ) == 0 && la.FullLength( ) )
               {    xout << "perfect match to " << la.pos2( ) << "-" 
                         << la.Pos2( ) << "\n";
                    continue;    }
               xout << la.mutations << " mismatches/" << la.indels << " indels, ";
               if ( !la.FullLength( ) )
               {    xout << "from " << la.pos1( ) << "-" << la.Pos1( ) 
                         << " (of " << la.query_length << ") ";    }
               xout << "to " << la.pos2( ) << "-" << la.Pos2( ) << "\n";
               la.PrintVisual( xout, edges[la.query_id], query_qual,
                    genome[la.target_id], 0 );    }
          xout << endl;
          String s = xout.str( );
          if ( s.size( ) > 0 ) s.resize( s.size( ) - 1 );
          align_descrip[i] = s;    }

     // Print.

     h.PrintSummaryPlus( out, 0, 0, 0, 0, 0, False, &align_descrip, 
          print_component_id_line, &component_remarks, 
          &component_remarks_only );
     out << "\n";    }

template void digraphE<KmerPath>::ReorderComponents(vec<int> const&);

/**
   Function: TrustedPathsToIndexedAligns
   
   Convert a vec of TrustedPaths into an indexed vec of look_aligns.

   So, for each <component> we have a set of <trusted paths> that it
   captures.   Each path has a set of edges, so make a flat list of
   alignments of those edges, and an index from each edge of its alignments
   to the reference.
*/
void TrustedPathsToIndexedAligns( const vec<TrustedPath>& paths,
                                  const int numEdges,
                                  vec<look_align>& aligns,
                                  vec< vec<int> >& aligns_index, 
                                  const bool removeImproper )
{
  aligns.clear();

  unsigned int numAligns = 0;
  for ( unsigned int i = 0; i < paths.size(); ++i )
    numAligns += paths[i].GetNumAligns();
  aligns.reserve( numAligns );

  for ( unsigned int i = 0; i < paths.size(); ++i )
    copy( paths[i].GetAllAligns().begin(), paths[i].GetAllAligns().end(),
          back_inserter( aligns ) );

  if ( removeImproper )
    aligns.erase( remove_if( aligns.begin(), aligns.end(),
                             not1( mem_fun_ref( &look_align::FullLength ) ) ),
                  aligns.end() );

  aligns_index.clear();
  aligns_index.resize( numEdges );
  for ( int i = 0; i < aligns.isize( ); i++ )
    aligns_index[ aligns[i].query_id ].push_back(i);
}

// Function: FilterByReference
//
// For each <component> of a HyperKmerPath, create <trusted paths> through the component
// that follow the component's graph structure but also align well to the reference.
//
// Find edge alignments and filter them by finding
// the alignments most consistent with the graph structure and the
// finished sequence.  These alignments are then packaged into <trusted paths>.
//
void FilterByReference( const HyperKmerPath& theGraph, 
                        const int K,
                        const vec<look_align>& aligns,
                        const vec< vec<align_id_t> >& aligns_index,
                        vec<TrustedPath>& trustedPaths )
{
  trustedPaths.clear();

  int numEdges = theGraph.EdgeObjectCount( );
  int numVerts = theGraph.N();

  // Set up reverse lookup table so we can easily find the vertices
  // associated with a given edge.
  vec<int> sources( numEdges, -1 );
  vec<int> targets( numEdges, -1 );
  
  for ( int v = 0; v < theGraph.N(); ++v ) {
    for ( int i = 0; i < theGraph.FromEdgeObj( v ).isize(); ++i )
      sources[ theGraph.FromEdgeObj(v)[i] ] = v;
    for ( int i = 0; i < theGraph.ToEdgeObj( v ).isize(); ++i )
      targets[ theGraph.ToEdgeObj(v)[i] ] = v;
  }
   
  // Find the contigs (connected components) of the graph.
  equiv_rel contigs;
  theGraph.ComponentRelation( contigs );
  vec<int> sampleVerts;
  contigs.OrbitRepsAlt( sampleVerts );
  int numContigs = sampleVerts.size();

  // Process the contigs one at a time.
  for ( int contig = 0; contig < numContigs; ++contig ) {

    int sampleVert = sampleVerts[contig];

    vec<int> contigVerts;
    contigs.Orbit( sampleVert, contigVerts );

    // Copy aligns for given contig:
    // alignments of this contig's edges to the reference.
    // Note that each edge may in general align to many places in the reference
    // (even to different <genome parts>).
    // Also, determine the length and id of the longest edge. 
    vec<look_align> rawAligns;

    nbases_t totalEdgeLen = 0;
    unsigned int maxEdgeLen = 0;
    int maxEdgeId = -1;
    for ( int i = 0; i < contigVerts.isize(); ++i ) {
      vec<int> fromEdges( theGraph.FromEdgeObj( contigVerts[i] ) );
      for ( int j = 0; j < fromEdges.isize(); ++j ) {
        int edgeId = fromEdges[j];

        totalEdgeLen += theGraph.EdgeObject( edgeId ).KmerCount() + (K-1);

        vec<align_id_t> alignIdxs = aligns_index[edgeId];

        if ( ! alignIdxs.empty() ) {
          if ( maxEdgeLen < aligns[ alignIdxs[0] ].query_length ) {
            maxEdgeLen = aligns[ alignIdxs[0] ].query_length;
            maxEdgeId = edgeId;
          }
          
          for ( int k = 0; k < alignIdxs.isize(); ++k ) 
            rawAligns.push_back( aligns[ alignIdxs[k] ] );
        }
      }
    }

    sort( rawAligns.begin(), rawAligns.end(),
          order_lookalign_TargetBegin() );

    vec<Bool> isSeed( rawAligns.size(), false );

    // Flag perfect, full-length alignments and all alignments to the
    // longest edge as possible seeds.
    for ( align_id_t thisIdx = 0; thisIdx < rawAligns.isize(); ++thisIdx ) {
      look_align& anAlign = rawAligns[thisIdx];
      if ( anAlign.query_id == maxEdgeId ||
           anAlign.Errors() == 0 && 
           anAlign.FullLength() ) 
        isSeed[ thisIdx ] = true;
    }
    
    const int nullSeed = -1;
    vec<align_id_t> fromSeed( rawAligns.size(), nullSeed );

    // Start at each seed align, and work forwards, flagging aligns
    // as trusted if they match the graph structure, stopping when you
    // hit an already-trusted align.
    for ( unsigned int seedIdx = 0; seedIdx < isSeed.size(); ++seedIdx ) {
      // If this is not a seed, skip it.
      if ( ! isSeed[seedIdx] )
        continue;
      
      // If this align has already been claimed by some other seed, skip it.
      if ( fromSeed[seedIdx] != nullSeed )
        continue;
      
      // Claim this align for this seed.
      fromSeed[seedIdx] = seedIdx;

      look_align* pLastTrustedAlign = &rawAligns[seedIdx];
      int lastTrustedEdgeIdx = pLastTrustedAlign->query_id;
      genome_part_id_t targetId = pLastTrustedAlign->target_id;

      // What will be the pos2 (start on the reference) of the next align?
      // It will be at the end of the last trusted align
      // (with an adjustment for the fact that edges are KmerPaths in kmer space,
      // so KmerPaths from adjacent edges overlap by K-1 bases).
      genome_part_pos_t pos2Target = pLastTrustedAlign->Pos2() - (int)(K-1);
      
      int nextVertex;
      vec<int> validEdges;
      if ( pLastTrustedAlign->rc1 ) {
        nextVertex = sources[ lastTrustedEdgeIdx ];
        validEdges = theGraph.ToEdgeObj( nextVertex );
      } else {
        nextVertex = targets[ lastTrustedEdgeIdx ];
        validEdges = theGraph.FromEdgeObj( nextVertex );
      }
      sort( validEdges.begin(), validEdges.end() );  // so we can binary_search it

      // Find edge alignments to the reference, that align an edge of this contig
      // to the endpoint of the current trusted path.
      // Remember that rawAligns is sorted by the genome part and within a genome part
      // by the start of the alignment on that genome part.
      for ( unsigned int forwIdx = seedIdx+1; forwIdx < rawAligns.size(); ++forwIdx ) {
        // Skip already claimed aligns.
        if ( fromSeed[forwIdx] >= 0 )
          continue;

        look_align* pThisAlign = &rawAligns[forwIdx];

	// If we're getting alignments to a different genome part, stop --
	// because of how we sorted rawAligns, we won't see any more alignments
	// to our genome part further down in the array.
        if ( longlong(pThisAlign->target_id) != longlong(targetId) )
          break;

	// If we're getting alignments that start earlier than the end of the
	// current trusted path, keep looking for alignments that start
	// further to the right on the genome part.
        else if ( pThisAlign->pos2() < pos2Target )
          continue;

	// If we're already getting alignments that start further to the right
	// on the genome part than we need, stop -- because of how we sorted rawAligns,
	// we won't find any more alignments to the end of the current trusted path. 
        else if ( pThisAlign->pos2() > pos2Target )
          break;

        else if ( pThisAlign->pos2() == pos2Target ) {
	  // Aha! found an alignment of an edge to the point on the reference
	  // where the current trusted path ends.

	  // First of all, find _all_ such alignments of edges.  They will all be
	  // together in rawAligns, because of how we sorted rawAligns.
	  // We will then pick one of these edges to extend our trusted trusted path
	  // through the contig.
          unsigned int candidatesBegin = forwIdx;
          unsigned int candidatesEnd = forwIdx+1;
          while ( candidatesEnd < rawAligns.size() &&
                  rawAligns[candidatesEnd].pos2() == pos2Target )
            ++candidatesEnd;
          
          Float lowestErrorRate = -1;
          unsigned int bestAlignIdx = candidatesBegin;

          for ( unsigned int candidateIdx = candidatesBegin;
                candidateIdx < candidatesEnd; 
                ++candidateIdx ) {
            pThisAlign = &rawAligns[candidateIdx];
            
            if ( pThisAlign->rc1 == pLastTrustedAlign->rc1 &&
                 pThisAlign->pos1() == 0 &&
                 binary_search( validEdges.begin(), validEdges.end(), pThisAlign->query_id ) ) {

	      // So, this candidate edge aligns to the reference in the same orientation as the
	      // last edge of this trusted trusted path, and the alignment of the edge to the reference
	      // starts at the beginning of the edge, and this edge happens to be a successor in the graph
	      // of our last edge.

              int numErrors = pThisAlign->Errors();

              if ( ! pThisAlign->FullLength() )
                numErrors += (int) pThisAlign->query_length -
                  ( pThisAlign->Pos1() - pThisAlign->pos1() );

              Float errorRate = Float(numErrors)/Float(pThisAlign->query_length);

              if ( abs(lowestErrorRate) < 1e-6 &&
                   abs(errorRate) < 1e-6 ) {
		/*
                cout << "Warning: duplicate edges found out of vertex " 
                     << nextVertex << ":" << endl;

                cout << "  -> " << targets[ pThisAlign->query_id ] << " "
                     << "\"" << BaseAlpha( pThisAlign->query_id ) << "\" "
                     << theGraph.EdgeObject( pThisAlign->query_id ) << endl;

                look_align* pOtherAlign = &rawAligns[bestAlignIdx];
                cout << "  -> " << targets[ pOtherAlign->query_id ] << " "
                     << "\"" << BaseAlpha( pOtherAlign->query_id ) << "\" "
                     << theGraph.EdgeObject( pOtherAlign->query_id ) << endl;
		*/
              }

              if ( lowestErrorRate < Float(0) ||
                   lowestErrorRate > errorRate ) {
                lowestErrorRate = errorRate;
                bestAlignIdx = candidateIdx;
              }
            }  // if this candidate edge is a successor in the graph of our last edge
          }  // for each candidate edge (an edge aligned to the end of the trusted-path-so-far).

          // If lowestErrorRate is still less than zero, we found no
          // valid align to a valid edge, so we're done.
          if ( lowestErrorRate < Float(0) )
            break;

          // Otherwise, bestAlignIdx is the align we'll follow.
          pThisAlign = &rawAligns[bestAlignIdx];
          fromSeed[ bestAlignIdx ] = seedIdx;

          // The for loop will increment this.
          forwIdx = candidatesEnd-1;

          // Set up the criteria for the next align.
          lastTrustedEdgeIdx = pThisAlign->query_id;
          pLastTrustedAlign = pThisAlign;
          if ( pLastTrustedAlign->rc1 ) {
            int nextVertex = sources[ lastTrustedEdgeIdx ];
            validEdges = theGraph.ToEdgeObj( nextVertex );
          } else {
            int nextVertex = targets[ lastTrustedEdgeIdx ];
            validEdges = theGraph.FromEdgeObj( nextVertex );
          }
          sort( validEdges.begin(), validEdges.end() );

          pos2Target = pLastTrustedAlign->Pos2() - (int)(K-1);
        }  // else if ( pThisAlign->pos2() == pos2Target ) ...
      }  // for ( unsigned int forwIdx = seedIdx+1; forwIdx < rawAligns.size(); ++forwIdx ) 
    }  // from each unclaimed seed alignment, work forwards

    SortSync( rawAligns, fromSeed, order_lookalign_TargetEnd() );
      
    set<align_id_t> seedsBackwalked;

    // Now work backwards from the seed (the alignment of an edge to the reference).
    // A trusted path is then the concatenation of walking backwards from the seed
    // as far as we can, and walking forward as far as we can.
    for ( unsigned int alignIdx = 0; alignIdx < fromSeed.size(); ++alignIdx ) {
      align_id_t seedIdx = fromSeed[alignIdx];

      // Skip aligns that aren't in a path.
      if ( seedIdx < 0 )
        continue;

      // Skip already processed seeds.
      if ( seedsBackwalked.count( seedIdx ) )
        continue;
      seedsBackwalked.insert( seedIdx );
      
      look_align* pLastTrustedAlign = &rawAligns[alignIdx];
      unsigned int lastTrustedEdgeIdx = pLastTrustedAlign->query_id;

      // Note that valid edge determination is reversed here.
      int nextVertex;
      vec<int> validEdges;
      if ( pLastTrustedAlign->rc1 ) {
        nextVertex = targets[ lastTrustedEdgeIdx ];
        validEdges = theGraph.FromEdgeObj( nextVertex );
      } else {
        nextVertex = sources[ lastTrustedEdgeIdx ];
        validEdges = theGraph.ToEdgeObj( nextVertex );
      }
      sort( validEdges.begin(), validEdges.end() );
      
      int Pos2Target = pLastTrustedAlign->pos2() + (int)(K-1);
      int targetId = pLastTrustedAlign->target_id;

      for ( int backIdx = alignIdx-1; backIdx >= 0; --backIdx ) {
        if ( fromSeed[backIdx] >= 0 )
          continue;
        
        look_align* pThisAlign = &rawAligns[backIdx];
        if ( pThisAlign->target_id != targetId )
          break;

        else if ( pThisAlign->Pos2() > Pos2Target )
          continue;

        else if ( pThisAlign->Pos2() < Pos2Target )
          break;

        else if ( pThisAlign->Pos2() == Pos2Target ) {
          int candidatesEnd = backIdx+1;
          int candidatesBegin = backIdx-1;
          while ( candidatesBegin >= 0 &&
                  rawAligns[candidatesBegin].Pos2() == Pos2Target )
            --candidatesBegin;
          ++candidatesBegin;

          Float lowestErrorRate = -1;
          unsigned int bestAlignIdx = candidatesBegin;

          for ( int candidateIdx = candidatesBegin;
                candidateIdx < candidatesEnd; 
                ++candidateIdx ) {
            pThisAlign = &rawAligns[candidateIdx];
            
            if ( pThisAlign->rc1 == pLastTrustedAlign->rc1 &&
                 pThisAlign->Pos1() == (int)pThisAlign->query_length &&
                 binary_search( validEdges.begin(), validEdges.end(), pThisAlign->query_id ) ) {

              int numErrors = pThisAlign->Errors();

              if ( ! pThisAlign->FullLength() )
                numErrors += (int) pThisAlign->query_length -
                  ( pThisAlign->Pos1() - pThisAlign->pos1() );

              Float errorRate = Float(numErrors)/Float(pThisAlign->query_length);

              if ( lowestErrorRate == Float(0) &&
                   errorRate == Float(0) ) {
                cout << "Warning: duplicate edges found out of vertex " 
                     << nextVertex << ":" << endl;

                cout << "  " << BaseAlpha( pThisAlign->query_id ) << " "
                     << theGraph.EdgeObject( pThisAlign->query_id ) << endl;

                look_align* pOtherAlign = &rawAligns[bestAlignIdx];
                cout << "  " << BaseAlpha( pOtherAlign->query_id ) << " "
                     << theGraph.EdgeObject( pOtherAlign->query_id ) << endl;
              }

              if ( lowestErrorRate < Float(0) ||
                   lowestErrorRate > errorRate ) {
                lowestErrorRate = errorRate;
                bestAlignIdx = candidateIdx;
              }
            }
          }

          if ( lowestErrorRate < Float(0) )
            break;
            
          pThisAlign = &rawAligns[bestAlignIdx];
          fromSeed[ bestAlignIdx ] = seedIdx;

          backIdx = candidatesBegin;

          lastTrustedEdgeIdx = pThisAlign->query_id;
          pLastTrustedAlign = pThisAlign;
          // Note that valid edge determination is reversed here.
          if ( pLastTrustedAlign->rc1 ) {
            int nextVertex = targets[ lastTrustedEdgeIdx ];
            validEdges = theGraph.FromEdgeObj( nextVertex );
          } else {
            int nextVertex = sources[ lastTrustedEdgeIdx ];
            validEdges = theGraph.ToEdgeObj( nextVertex );
          }
          sort( validEdges.begin(), validEdges.end() );

          Pos2Target = pLastTrustedAlign->pos2() + (int)(K-1);
        }
      }
    }
    
    vec<int> vertexIdsInPath;
    vec<look_align> alignsInPath;
    
    set<align_id_t> seedsSaved;
    seedsSaved.insert( nullSeed );

    for ( unsigned int i = 0; i < rawAligns.size(); ++i ) {
      int seedIdx = fromSeed[i];

      if ( seedsSaved.count( seedIdx ) )
        continue;

      seedsSaved.insert( seedIdx );

      const look_align& firstAlign = rawAligns[i];

      if ( firstAlign.rc1 )
        vertexIdsInPath.push_back( targets[ firstAlign.query_id ] );
      else
        vertexIdsInPath.push_back( sources[ firstAlign.query_id ] );

      for ( unsigned int j = i; j < rawAligns.size(); ++j ) {
        if ( fromSeed[j] != seedIdx ) 
          continue;
        
        const look_align& thisAlign = rawAligns[j];

        alignsInPath.push_back( thisAlign );

        if ( thisAlign.rc1 )
          vertexIdsInPath.push_back( sources[ thisAlign.query_id ] );
        else
          vertexIdsInPath.push_back( targets[ thisAlign.query_id ] );
      }
        
      trustedPaths.push_back( TrustedPath( contig, totalEdgeLen, vertexIdsInPath, alignsInPath ) );
      vertexIdsInPath.clear();
      alignsInPath.clear();
    }  // for all rawAligns
  }  // for each contig
  
  sort( trustedPaths.begin(), trustedPaths.end() );
}  // FilterByReference()

// Constructor for TrustedPath.

TrustedPath::TrustedPath( int contig,
                          int contigTotalEdgeLength,
                          const vec<int>& vertexIds,
                          const vec<look_align>& aligns )
  : m_contig( contig ),
    m_vertexIds( vertexIds ),
    m_aligns( aligns )
{
  ForceAssertEq( m_vertexIds.size() - 1, m_aligns.size() );
  
  map<int,int> maxAlignLengthPerEdge;
  for ( unsigned int i = 0; i < aligns.size(); ++i ) {
    int edgeId = aligns[i].query_id;
    int alignLen = aligns[i].Pos1() - aligns[i].pos1();
    
    map<int,int>::iterator found = maxAlignLengthPerEdge.find( edgeId );
    if ( found == maxAlignLengthPerEdge.end() )
      maxAlignLengthPerEdge.insert( make_pair( edgeId, alignLen ) );
    else
      if ( found->second < alignLen )
        found->second = alignLen;
  }
  
  int edgeLengthAligned = 0;
  for ( map<int,int>::iterator iEdge = maxAlignLengthPerEdge.begin();
        iEdge != maxAlignLengthPerEdge.end(); ++iEdge )
    edgeLengthAligned += iEdge->second;
  
  m_fractionEdgeLengthAligned = (float)edgeLengthAligned / (float)contigTotalEdgeLength;
}

void TrustedPath::PrintSummary( ostream& out ) const 
{
  int oldPrecision = out.precision();
  out.precision(3);
  out << "contig " << GetContig() << ": " 
      << Begin() << "-" << End() << " on " << GetFinishedId() 
      << " (" << End()-Begin() << " bases" 
      << ", " << setprecision(3) << GetFractionEdgeLengthAligned() * Float(100.0) << "%)"
      << endl;
  out.precision(oldPrecision);
}

void TrustedPath::TestValid( ) const {
  ForceAssert( !m_aligns.empty() );
  ForceAssertEq( m_vertexIds.size() - 1, m_aligns.size() );
  // check that all edges of a path align to the same genome part,
  // to the same strand of that genome part, and follow each other
  // (are adjacent to each other) on that strand.
  
  genome_part_id_t genomePart = -1;
  Bool rc = False;
  for ( int i = 0; i < m_aligns.isize(); i++ ) {
    const look_align& a = m_aligns[i];
    // check that everything aligns to the same genome part
    if ( genomePart == -1 ) {
      genomePart = a.target_id;
      rc = a.rc1;
    } else {
      ForceAssertEq( genomePart, a.target_id );
      ForceAssertEq( rc, a.rc1 );
    }

  }
}

/**
   Method: TestValid

   Test the validity of this TrustedPath -- apply additional tests
   beyond TestValid(), that are not available without having the
   HyperKmerPath.
*/
void TrustedPath::TestValid( const HyperKmerPath& h ) const {
  TestValid( );

  int K = h.K();
  for ( int i = 0; i < m_aligns.isize(); i++ ) {
    if ( i < m_aligns.isize()-1 ) {
      // check that the next edge starts where this edge stops on the reference
      ForceAssertEq( m_aligns[i].Pos2() + 1, m_aligns[i+1].pos2() + K );
    }
  }
  
  vec<int> edgeSources, edgeTargets;
  h.ToLeft( edgeSources );
  h.ToRight( edgeTargets );
  for ( int i = 0; i < m_aligns.isize(); i++ ) {
    const look_align& edgeAlign = m_aligns[i];
    int e = int(edgeAlign.query_id);
    if ( !edgeAlign.rc1 ) {
      ForceAssertEq( m_vertexIds[i], edgeSources[e] );
      ForceAssertEq( m_vertexIds[i+1], edgeTargets[e] );
    } else {
      ForceAssertEq( m_vertexIds[i], edgeTargets[e] );
      ForceAssertEq( m_vertexIds[i+1], edgeSources[e] );
    }
  }
}

// Output operator for TrustedPath.

ostream& operator<< ( ostream& out, const TrustedPath& path ) 
{
  path.PrintSummary( out );

  int numAligns = path.GetNumAligns();
  for ( int i = 0; i < numAligns; ++i ) {
    const look_align& thisAlign = path.GetAlign(i);
    if ( i == 0 && ! thisAlign.FullLength() )
      out << "...";
    else
      out << path.GetVertexIdBefore( i );
    out << endl;
    out << "  " << thisAlign.pos2() << "-" << thisAlign.Pos2() 
        << " " << BaseAlpha( thisAlign.query_id ) << ( thisAlign.rc1 ? "-" : "+" ) 
        << " (" << thisAlign.Pos1() - thisAlign.pos1() << " of " << thisAlign.query_length << " bases"
        << ", " << thisAlign.mutations << " mutations" 
        << ", " << thisAlign.indels << " indels";
    if ( ! thisAlign.FullLength() )
      out << ", incomplete";
    out << ")"
        << endl;
  }
  const look_align& lastAlign = path.GetAlign( numAligns-1 );
  if ( ! lastAlign.FullLength() )
    out << "...";
  else
    out << path.GetLastVertexId();
  out << endl;

  return out;
}

// A TrustedPath dominates another's aligns if both paths are from
// the same contig and the finished sequence covered by this path
// extends past that covered by the other path in one or both
// directions.

bool TrustedPath::DominatesAlignsOf( const TrustedPath& other ) const {

  if ( this->m_contig != other.m_contig ) 
    return false;

  if ( this->GetFinishedId() != other.GetFinishedId() )
    return false;

  if ( this->Begin() <= other.Begin() && this->End() > other.End() ||
       this->Begin() < other.Begin() && this->End() >= other.End() )
    return true;
  
  return false;
}


const vec<int>& 
TrustedPath::GetUniquePerfectEdgeIds() const {

  if ( m_uniquePerfectEdgeIds.empty() ) {
    for ( unsigned int i = 0; i < this->m_aligns.size(); ++i ) {
      const look_align& thisAlign = this->m_aligns[i];
      if ( thisAlign.indels == 0 && 
           thisAlign.mutations == 0 &&
           thisAlign.FullLength() )
        m_uniquePerfectEdgeIds.push_back( thisAlign.query_id );
    }
    UniqueSort( m_uniquePerfectEdgeIds );
  }
  
  return m_uniquePerfectEdgeIds;
}
    

// A TrustedPath dominates another's edges if the perfect, full-length
// alignments in the other path involve a strict subset of the edges
// involved in the perfect, full-length alignments of the dominating
// path.

bool TrustedPath::DominatesEdgesOf( const TrustedPath& other ) const {

  if ( this->m_contig != other.m_contig ) 
    return false;

  const vec<int>& thisEdgeIds = this->GetUniquePerfectEdgeIds();
  const vec<int>& otherEdgeIds = other.GetUniquePerfectEdgeIds();

  if ( thisEdgeIds.size() <= otherEdgeIds.size() )
    return false;

  vec<int>::const_iterator iThisEdge = thisEdgeIds.begin();
  vec<int>::const_iterator iOtherEdge = otherEdgeIds.begin();

  while ( iThisEdge != thisEdgeIds.end() &&
          iOtherEdge != otherEdgeIds.end() ) {
    if ( *iThisEdge < *iOtherEdge )
      ++iThisEdge;
    else if ( *iThisEdge == *iOtherEdge ) {
      ++iThisEdge;
      ++iOtherEdge;
    }
    else 
      return false;
  }

  return ( iOtherEdge == otherEdgeIds.end() );
}


// A TrustedPath dominates another's vertices if the vertices of the
// other path are a strict subset of the vertices of the dominating
// path.

bool TrustedPath::DominatesVerticesOf( const TrustedPath& other ) const {

  if ( this->m_contig != other.m_contig ) 
    return false;

  vec<int> thisVertexIds = this->m_vertexIds;
  if ( ! this->GetAllAligns().front().FullLength() )
    thisVertexIds.erase( thisVertexIds.begin() );
  if ( ! this->GetAllAligns().back().FullLength() )
    thisVertexIds.erase( thisVertexIds.end() - 1 );
  UniqueSort( thisVertexIds );

  vec<int> otherVertexIds = other.m_vertexIds;
  if ( ! other.GetAllAligns().front().FullLength() )
    otherVertexIds.erase( otherVertexIds.begin() );
  if ( ! other.GetAllAligns().back().FullLength() )
    otherVertexIds.erase( otherVertexIds.end() - 1 );
  UniqueSort( otherVertexIds );

  if ( thisVertexIds.size() <= otherVertexIds.size() )
    return false;

  vec<int>::iterator iThisVertex = thisVertexIds.begin();
  vec<int>::iterator iOtherVertex = otherVertexIds.begin();

  while ( iThisVertex != thisVertexIds.end() &&
          iOtherVertex != otherVertexIds.end() ) {
    if ( *iThisVertex < *iOtherVertex )
      ++iThisVertex;
    else if ( *iThisVertex == *iOtherVertex ) {
      ++iThisVertex;
      ++iOtherVertex;
    }
    else 
      return false;
  }

  return ( iOtherVertex == otherVertexIds.end() );
}


struct order_TrustedPath_byContig 
  : public binary_function<TrustedPath,TrustedPath,bool> 
{
  bool operator() ( const TrustedPath& lhs, const TrustedPath& rhs ) const 
  {
    return ( lhs.GetContig() < rhs.GetContig() );
  }
};


struct order_TrustedPath_byContigAndGenomePart
  : public binary_function<TrustedPath,TrustedPath,bool> 
{
  bool operator() ( const TrustedPath& lhs, const TrustedPath& rhs ) const 
  {
    return
      lhs.GetContig() < rhs.GetContig() ||
      lhs.GetContig() == rhs.GetContig() && lhs.GetFinishedId() < rhs.GetFinishedId() ;
  }
};


template <typename DominanceFunctor>
void FilterPathsByDominance( vec<TrustedPath>& trustedPaths, DominanceFunctor functor )
{
  cout << Date() << ": sorting " << trustedPaths.isize() << " trusted paths..." << endl;
  Sort( trustedPaths, order_TrustedPath_byContig() );

  vector<bool> isDominated( trustedPaths.size(), false );
  cout << Date() << ": filtering paths by dominance..." << endl;

  unsigned int contigBegin = 0;
  while ( contigBegin < trustedPaths.size() ) {
    unsigned int contigEnd = contigBegin + 1;
    while ( contigEnd < trustedPaths.size() &&
            trustedPaths[contigEnd].GetContig() == trustedPaths[contigBegin].GetContig() )
      ++contigEnd;
    // cout << Date() << ": contig " << trustedPaths[contigBegin].GetContig() 
    //      << " : " << (contigEnd - contigBegin) << " trusted paths." << endl;
    int nfunctorCalls = 0;
    int ndominatedMarked = 0;
    for ( unsigned int i = contigBegin; i < contigEnd-1; ++i )
      for ( unsigned int j = i+1; j < contigEnd && !isDominated[i]; ++j )
	if ( !isDominated[j] ) {
	  nfunctorCalls++;
	  if ( functor( trustedPaths[i], trustedPaths[j] ) ) {
	    isDominated[ j ] = true;
	    ndominatedMarked++;
	  } 
	  else  {
	    nfunctorCalls++;
	    if ( functor( trustedPaths[j], trustedPaths[i] ) ) {
	      isDominated[ i ] = true;
	      ndominatedMarked++;
	    }
	  }
	}
    cout << Date() << ":   --> done with contig " << trustedPaths[contigBegin].GetContig() << " : " << (contigEnd - contigBegin) << " trusted paths; made " <<
      nfunctorCalls << " functor calls; marked " << ndominatedMarked << " dominated paths ( " << (contigEnd - contigBegin - ndominatedMarked) << " left." << endl;
    contigBegin = contigEnd;
  }
  

  vec<TrustedPath> savedPaths;

  for ( unsigned int i = 0; i < trustedPaths.size(); ++i )
    if ( ! isDominated[i] )
      savedPaths.push_back( trustedPaths[i] );
    
  trustedPaths.swap( savedPaths );
  Sort( trustedPaths );
  cout << Date() << "   --> after filtering: " << trustedPaths.size() << " trusted paths left." << endl;
}

// Type: begend_t
// A trusted path id, together with "begin" or "end" marker.  Used in
// FilterPathsByAlignDominance().
typedef pair< tpid_t, Bool > begend_t;

/**
   Functor: begend_cmp

   Used to sort begin/end markers in FilterPathsByAlignDominance().
*/
template <int rightEndPlus = 0>
struct begend_cmp {
private:
  const vec<TrustedPath>& trustedPaths;
public:
  begend_cmp( const vec<TrustedPath>& _trustedPaths):
    trustedPaths(_trustedPaths) { }
  
  bool operator() ( const begend_t& be1, const begend_t& be2 ) const {
    const TrustedPath& p1 = trustedPaths[be1.first];
    const TrustedPath& p2 = trustedPaths[be2.first];
    Bool be1isBeg = be1.second;
    Bool be2isBeg = be2.second;
    genome_part_pos_t be1pos = be1isBeg ? p1.Begin() : ( p1.End() + rightEndPlus );
    genome_part_pos_t be2pos = be2isBeg ? p2.Begin() : ( p2.End() + rightEndPlus );
    // first of all, sort by position on reference
    if ( be1pos < be2pos ) return True;
    if ( be1pos > be2pos ) return False;
    // if position is the same:
    //   ... put begins first, so that the begin marker of a path is always before
    // the end marker even if they're at the same base
    if ( be1isBeg && !be2isBeg ) return True;
    if ( !be1isBeg && be2isBeg ) return False;
    // if they're both begin markers at the same position, put the marker
    // whose end is further to the right first.
    if ( be1isBeg ) {
      if ( p1.End() > p2.End() ) return True;
      if ( p1.End() < p2.End() ) return False;
    } 
    // if all else fails, sort by the path id
    return be1.first < be2.first;
  }
};


genome_part_pos_t BegEndPos( const begend_t& be, const vec<TrustedPath>& trustedPaths ) {
  const TrustedPath& tp = trustedPaths[ be.first ];
  Bool isBeg = be.second;
  return ( isBeg ? tp.Begin() : tp.End()+1 );
}

/**
   Functor: begend_cmp2

   Used to sort begin/end markers in FilterPathsByAlignDominance().
*/
struct begend_cmp2 {
private:
  const vec<TrustedPath>& trustedPaths;
  
public:
  begend_cmp2( const vec<TrustedPath>& _trustedPaths ):
    trustedPaths(_trustedPaths) { }
  
  bool operator() ( const begend_t& be1, const begend_t& be2 ) const {
    const TrustedPath& p1 = trustedPaths[be1.first];
    const TrustedPath& p2 = trustedPaths[be2.first];
    Bool be1isBeg = be1.second;
    Bool be2isBeg = be2.second;
    genome_part_pos_t be1pos = BegEndPos( be1, trustedPaths );
    genome_part_pos_t be2pos = BegEndPos( be2, trustedPaths );

    // first of all, sort by position on reference
    if ( be1pos < be2pos ) return True;
    if ( be1pos > be2pos ) return False;
    // if position is the same:
    //   ... put begins first, so that the begin marker of a path is always before
    // the end marker even if they're at the same base
    if ( be1isBeg && !be2isBeg ) return True;
    if ( !be1isBeg && be2isBeg ) return False;
    // if they're both begin markers at the same position, put the marker
    // whose end is further to the right first.
    if ( be1isBeg ) {
      if ( p1.End() > p2.End() ) return True;
      if ( p1.End() < p2.End() ) return False;
    } 
    // if all else fails, sort by the path id
    return be1.first < be2.first;
  }
};


/**
   Function: FilterPathsByAlignDominance

   Remove paths that cover finished sequence that is covered by some
   other, longer path from that contig.

   Algorithm:

   Stratify paths by contig and (within the contig) by the genome part to which they align.
   For paths within a contig that align to the same genome part:
   Sort the begins & ends of paths along the genome part.
   (So, for each path we create a begin marker and an end marker, and sort the combined vector of all markers
   by position on the reference).
   Go from left to right in this sorted vector, keeping track of the path that extends furthest
   to the right of the paths covering the current position.
   If we see a begin marker, this path is either dominated by the current-rightmost-path or becomes
   the current-rightmost-path.
*/
void FilterPathsByAlignDominance( vec<TrustedPath>& trustedPaths )
{
  
  vec<Bool> isDominated( trustedPaths.size(), False );
  
  Sort( trustedPaths, order_TrustedPath_byContigAndGenomePart() );

  int sectionBegin = 0;
  while ( sectionBegin < trustedPaths.isize() ) {
    int sectionEnd = sectionBegin + 1;
    while ( sectionEnd < trustedPaths.isize() &&
            trustedPaths[sectionEnd].GetContig() == trustedPaths[sectionBegin].GetContig() &&
	    trustedPaths[sectionEnd].GetFinishedId() == trustedPaths[sectionBegin].GetFinishedId() )
      ++sectionEnd;

    static vec< begend_t > begEnds;
    begEnds.clear();
    for ( tpid_t i = sectionBegin; i < sectionEnd; i++ ) {
      begEnds.push_back( make_pair( i, True ) );
      begEnds.push_back( make_pair( i, False ) );
    }
    Sort( begEnds, begend_cmp<0>(trustedPaths) );
    tpid_t rightmostEndingPathId = -1;
    
    for ( int i = 0; i < begEnds.isize(); i++ ) {
      const begend_t& begEnd = begEnds[i];
      tpid_t pathId = begEnd.first;
      Bool isBeg = begEnd.second;
      if ( isBeg ) {
	if ( rightmostEndingPathId != -1  &&
	     ( trustedPaths[ pathId ].End() < trustedPaths[ rightmostEndingPathId ].End() ||
	       trustedPaths[ pathId ].End() == trustedPaths[ rightmostEndingPathId ].End() &&
	       trustedPaths[ pathId ].Begin() > trustedPaths[ rightmostEndingPathId ].Begin() ) ) {
	  isDominated[ pathId ] = True;
	} else
	  rightmostEndingPathId = pathId;
      } else {  // if this is the endpoint of a path
	if ( pathId == rightmostEndingPathId )
	  rightmostEndingPathId = -1;
      }
    }
    
    sectionBegin = sectionEnd;
  }

  EraseIf( trustedPaths, isDominated );
}

/**
   Function: FilterPathsByEdgeDominance

   Remove paths that are edge-dominated by some other path from that contig.

   Algorithm:

   Path A can dominate path B only if they share an edge, so, first built a per-edge
   index of all paths that contain that edge, then for each edge look at all paths that
   contain that edge and see which of them dominates which.  If path A is dominated,
   don't bother checking if it dominates some other path B -- since domination is transitive
   anything dominated by A will be dominated by B.
*/
void FilterPathsByEdgeDominance( vec<TrustedPath>& trustedPaths, int numEdges )
{
  
  vec<Bool> isDominated( trustedPaths.size(), False );
  
  Sort( trustedPaths, order_TrustedPath_byContig() );  // redundant

  // Local var: uwpathsAtEdge
  // For each edge, the trusted paths that contain that edge.
  vec< vec< tpid_t > > uwpathsAtEdge(numEdges+1);

  tpid_t sectionBegin = 0;
  while ( sectionBegin < trustedPaths.isize() ) {
    tpid_t sectionEnd = sectionBegin + 1;
    while ( sectionEnd < trustedPaths.isize() &&
            trustedPaths[sectionEnd].GetContig() == trustedPaths[sectionBegin].GetContig() )
      ++sectionEnd;
    // cout << Date() << ": contig " << trustedPaths[sectionBegin].GetContig() 
    //      << " : " << (sectionEnd - sectionBegin) << " trusted paths." << endl;

    
    vec< int > componentEdges;
    for ( tpid_t i = sectionBegin; i < sectionEnd; i++ ) {
      const vec<int>& edgeIds = trustedPaths[i].GetUniquePerfectEdgeIds();
      for ( int j = 0; j < edgeIds.isize(); j++ ) {
	uwpathsAtEdge[ edgeIds[j] ].push_back( i );
	componentEdges.push_back( edgeIds[j] );
      }
    }
    UniqueSort( componentEdges );
    for ( int cei = 0; cei < componentEdges.isize(); cei++ ) {
      const vec<tpid_t>& pathsAtEdge = uwpathsAtEdge[ componentEdges[cei] ];
      for ( int i = 0; i < pathsAtEdge.isize()-1; i++ )
	for ( int j = i+1; j < pathsAtEdge.isize() && !isDominated[ pathsAtEdge[i] ]; j++ )
	  if ( !isDominated[ pathsAtEdge[j] ] ) {
	    if ( trustedPaths[ pathsAtEdge[i] ].DominatesEdgesOf( trustedPaths[ pathsAtEdge[j] ] ) )
	      isDominated[ pathsAtEdge[j] ] = True;
	    else if ( trustedPaths[ pathsAtEdge[j] ].DominatesEdgesOf( trustedPaths[ pathsAtEdge[i] ] ) )
	      isDominated[ pathsAtEdge[i] ] = True;
	  }
    }

    sectionBegin = sectionEnd;
  }
  EraseIf( trustedPaths, isDominated );
  
}


void FilterPathsByVertexDominance( vec<TrustedPath>& trustedPaths )
{
  FilterPathsByDominance( trustedPaths, mem_fun_ref( &TrustedPath::DominatesVerticesOf ) );
}


void FilterPathsByEdgeCoverage( vec<TrustedPath>& trustedPaths, 
                                const Float minFraction )
{
  vec<TrustedPath> savedPaths;

  for ( unsigned int i = 0; i < trustedPaths.size(); ++i )
    if ( ! ( trustedPaths[i].GetFractionEdgeLengthAligned() < minFraction ) )
      savedPaths.push_back( trustedPaths[i] );

  trustedPaths.swap( savedPaths );
}


void FilterPathsByLength( vec<TrustedPath>& trustedPaths, 
                          const int minLength, 
                          const int minPercent )
{
  vec<TrustedPath> savedPaths;

  for ( unsigned int i = 0; i < trustedPaths.size(); ++i )
    if ( trustedPaths[i].Length() >= minLength &&
         trustedPaths[i].GetFractionEdgeLengthAligned() >= Float(minPercent)/Float(100) )
      savedPaths.push_back( trustedPaths[i] );

  trustedPaths.swap( savedPaths );
}

// Function: ReportMisassemblies
//
// ReportMisassemblies.  Look for putative misassemblies:
//
// 1. Report edges that have no end-to-end alignment (or that have no alignment 
// at all).
//
// 2. For each component C, consider the parts of the genome that are covered by it.
// Divide these parts into their connected components, and consider those components
// that include a uniquely anchored edge of length >= MIN_LEN.  If the separation
// between two of these components is >= MIN_SEP, report C.  A given component is
// reported at most once.
void ReportMisassemblies( ostream& out, const HyperKmerPath& h,      
     const vec<look_align>& aligns, const vec< vec<int> >& aligns_index,
     const int MIN_LEN, const int MIN_SEP, const int MIN_LEN_NO_REPORT )
{    vec<String> events;

     // Check for edges that are not aligned or not aligned end-to-end.

     for ( int v = 0; v < h.N( ); v++ )
     {    for ( int j = 0; j < h.From(v).isize( ); j++ )
          {    int w = h.From(v)[j];
               int e = h.EdgeObjectIndexByIndexFrom( v, j );
               if ( h.EdgeLength(e) < MIN_LEN_NO_REPORT ) continue;
               if ( aligns_index[e].empty( ) )
               {    events.push_back( "An edge (" + BaseAlpha(e) + ")"
                         + " from vertex " + ToString(v) + " to " + ToString(w)
                         + " of length " + ToString( h.EdgeLength(e) )
                         + " has no consistent alignment." );
                    continue;    }
               Bool have_full_length = False;
               for ( int r = 0; r < aligns_index[e].isize( ); r++ )
               {    const look_align& la = aligns[ aligns_index[e][r] ];
                    if ( la.FullLength( ) ) have_full_length = True;    }
               if ( !have_full_length )
               {    events.push_back( "An edge (" + BaseAlpha(e) + ")"
                         + " from vertex " + ToString(v) + " to " + ToString(w) 
                         + " of length " + ToString( h.EdgeLength(e) ) + " has no "
                         + "consistent end-to-end alignment." );    }    }    }
     
     // Get number of <genome contigs> (that is known but not accessible here).
     // Similarly, get genome contig lengths (actually lower bounds).

     int genome_contigs = 0;
     for ( int i = 0; i < aligns.isize( ); i++ )
          genome_contigs = Max( genome_contigs, aligns[i].target_id + 1 );
     vec<int> genome_lengths(genome_contigs, 0);
     for ( int i = 0; i < aligns.isize( ); i++ )
     {    int t = aligns[i].target_id;
          genome_lengths[t] = Max( genome_lengths[t], aligns[i].Pos2( ) );    }

     // Find the components of the HyperKmerPath and representative vertices
     // from each.

     equiv_rel e;
     h.ComponentRelation(e);
     vec<int> reps, o;
     e.OrbitRepsAlt(reps);

     // Go through the components.

     for ( int i = 0; i < reps.isize( ); i++ )
     {    
          // Let o = {vertices in the component}.  Compute what is covered by the 
          // edges in the component, and also what is covered by long edges having
          // only a single alignment.

          e.Orbit( reps[i], o );
          vec< vec<ho_interval> > 
               covered(genome_contigs), covered_uniquely(genome_contigs);
          for ( int m = 0; m < o.isize( ); m++ )
          {    int v = o[m];
               for ( int j = 0; j < h.From(v).isize( ); j++ )
               {    int e = h.EdgeObjectIndexByIndexFrom( v, j );
                    for ( int r = 0; r < aligns_index[e].isize( ); r++ )
                    {    const look_align& la = aligns[ aligns_index[e][r] ];
                         ho_interval ho( la.pos2( ), la.Pos2( ) );
                         int t = la.target_id;
                         covered[t].push_back(ho);
                         if ( aligns_index[e].size( ) == 1 
                              && h.EdgeLength(e) >= MIN_LEN )
                         {    covered_uniquely[t].push_back(ho);    }    }    }    }

          // Find the connected components of the covered parts of the genome.
          // Restrict attention to those which subsume a long edge having a
          // unique alignment.  Look for two such components which are either on 
          // the same genome contig, but separated by at least MIN_SEP, or on 
          // different genome contigs.

          vec< vec<ho_interval> > cov2(genome_contigs);
          for ( int t = 0; t < genome_contigs; t++ )
          {    vec<ho_interval> cov;
               ExtractGivenCoverage( genome_lengths[t], 1, covered[t], cov );
               for ( int u = 0; u < cov.isize( ); u++ )
               {    for ( int x = 0; x < covered_uniquely[t].isize( ); x++ )
                    {    if ( Subset( covered_uniquely[t][x], cov[u] ) )
                         {    cov2[t].push_back( cov[u] );
                              break;    }    }    }
               for ( int j = 1; j < cov2[t].isize( ); j++ )
               {    const ho_interval &h1 = cov2[t][j-1], &h2 = cov2[t][j];
                    if ( h2.Start( ) >= h1.Stop( ) + MIN_SEP )
                    {    events.push_back( "Contig " + ToString(i) 
                              + " may have misjoin between "
                              + ToString(t) + "." + ToString( h1.Start( ) )
                              + "-" + ToString( h1.Stop( ) ) + " and "
                              + ToString(t) + "." + ToString( h2.Start( ) )
                              + "-" + ToString( h2.Stop( ) ) + "." );
                         goto next_component;    }    }
	  }
          next_component: continue;    }

     // Report the events.

     out << "\nPutative misassemblies detected:";
     if ( events.empty( ) ) 
     {    out << " none.\n";
          return;    }
     out << "\n";
     for ( int i = 0; i < events.isize( ); i++ )
          out << i+1 << ". " << events[i] << "\n";
     out << "\n";    }

void AlignAndPrintHyperKmerPath( ostream& out, const HyperKmerPath& h, 
     const KmerBaseBroker* kbb, const String& GENOME, const String& tmp_dir, 
     Bool print_component_id_line, bool filter )
{    static vec<look_align> aligns;
     static vec< vec<int> > aligns_index;
     static vec<TrustedPath> trusted_paths;
     AlignHyperKmerPath( h, kbb, GENOME, tmp_dir, aligns, aligns_index );
     FilterAligns( h, aligns, aligns_index, trusted_paths );
     vecbasevector genome( GENOME + ".fastb" );
     PrintAlignedHyperKmerPath( cout, h, kbb, genome, aligns,
          aligns_index, print_component_id_line, &trusted_paths );    }

/**
   Function: AlignTrustedPathsByGraph

   Find mutual alignments of <trusted paths> to each other, where the assembly graph
   implies a particular alignment of the two paths relative for each other (for example,
   they share an edge, or meet at a vertex).

   ( Not implemented yet ).
*/




// Local term: contig
// Here, it means "a connected component".
