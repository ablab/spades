/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include <strstream>

#include "Equiv.h"
#include "PrintAlignment.h"
#include "VecUtilities.h"
#include "feudal/IncrementalWriter.h"
#include "graph/DigraphTemplate.h"
#include "math/Functions.h"
#include "pairwise_aligners/SmithWatBandedA.h"
#include "paths/HyperKmerPath.h"
#include "paths/KmerBaseBroker.h"

// Template instantiations:

template digraphE<KmerPath>::digraphE( const digraphE& g, int n );
template digraphE<KmerPath>::digraphE( const digraphE& g, const vec< vec<int> >& C );
template digraphE<KmerPath>::digraphE(vec<KmerPath> const&,
     const ConstructorBehavior );
template digraphE<KmerPath>::digraphE(vec<KmerPath> const&, const equiv_rel&);
template digraphE<KmerPath>::digraphE( const vec< vec<int> >&, const vec< vec<int> >&,
          const vec<KmerPath>&, const vec< vec<int> >&, const vec< vec<int> >& );
template
     void digraphE<KmerPath>::SplitEdge(int, int, KmerPath const&, KmerPath const&);
template void digraphE<KmerPath>::Reverse();
template void digraphE<int>::Reverse();
template void digraphE<KmerPath>::ReverseComponent( int v );
template void digraphE<KmerPath>::ComponentEdges( vec< vec<int> >& edges ) const;
template void digraphE<KmerPath>::ComponentsE( vec< vec<int> >& comp ) const;
template void digraphE<KmerPath>::Glue(EmbeddedSubPath<KmerPath> const&,
     EmbeddedSubPath<KmerPath> const&, vec<int> const&, vec<int> const&,
     digraphE<KmerPath> const&);
template digraphE<KmerPath>::digraphE(digraphE<KmerPath> const&,
     vec<int> const&);
template void digraphE<KmerPath>::Used( vec<Bool>& used ) const;
template void digraphE<KmerPath>::DeleteEdges(const vec<int> &);
template void digraphE<KmerPath>::DeleteEdges(const vec<int>&,const vec<int>&);
template void digraphE<KmerPath>::DeleteEdgeTo(int,int);
template void digraphE<KmerPath>::DeleteEdgesAtVertex(int);
template void digraphE<KmerPath>::RemoveDuplicateEdges();
template void digraphE<KmerPath>::RemoveDeadEdgeObjects();
template void digraphE<KmerPath>::PopBubbles( const vec<int>& bubble_vs );
template void BinaryWrite(int, const digraphE<KmerPath>&);
template void BinaryRead(int, digraphE<KmerPath>&);
template Bool digraphE<KmerPath>::IsComplete( const vec<int>& vertices,
     const vec<int>& edges ) const;
template void digraphE<KmerPath>::DualComponentRelation(
     equiv_rel& e, const vec<Bool>& exclude ) const;
template digraphE<KmerPath>::digraphE(digraphE<KmerPath> const&, equiv_rel const&);
template void digraphE<KmerPath>::Initialize(vec<digraphE<KmerPath> > const&,
     vec<pair<pair<int, int>, pair<int, int> > > const&);
template void digraphE<KmerPath>::ToLeft(vec<int>&) const;
template void digraphE<KmerPath>::ToRight(vec<int>&) const;
template void digraphE<KmerPath>::PrettyDOT( ostream& out, const vec<double>& lengths, Bool label_contigs, Bool label_vertices, Bool label_edges, const vec<int>* componentsToPrint, const Bool, const vec<String> * , const vec<String> *, const vec<int> * ) const;
template void digraphE<KmerPath>::DumpGraphML( const String& ) const;





void HyperKmerPath::ReduceLoops( )
{    for ( int v = 0; v < N( ); v++ )
     {    if ( To(v).size( ) != 2 || From(v).size( ) != 1 ) continue;
          int w = From(v)[0];
          if ( v == w ) continue;
          int u;
          if ( To(v)[0] == w ) u = To(v)[1];
          else if ( To(v)[1] == w ) u = To(v)[0];
          else continue;
          if ( u == v || u == w ) continue;
          if ( From(w).size( ) != 2 || To(w).size( ) != 1 ) continue;
          int x;
          if ( From(w)[0] == v ) x = From(w)[1];
          else x = From(w)[0];
          if ( x == u || x == v || x == w ) continue;

          KmerPath vv = EdgeObjectByIndexFrom( v, 0 );
          vv.Append( EdgeObjectByIndexTo( v, ( To(v)[0] == w ? 0 : 1 ) ) );
          KmerPath vx = EdgeObjectByIndexFrom( v, 0 );
          vx.Append( EdgeObjectByIndexFrom( w, ( From(w)[0] == x ? 0 : 1 ) ) );
          int vx_id = EdgeObjectIndexByIndexTo( v, ( To(v)[0] == w ? 0 : 1 ) );
          EdgeObjectMutable(vx_id) = vx;
          FromMutable(w).clear( ), ToMutable(w).clear( );
          FromEdgeObjMutable(w).clear( ), ToEdgeObjMutable(w).clear( );
          int vv_id = EdgeObjectIndexByIndexFrom( v, 0 );
          EdgeObjectMutable(vv_id) = vv;
          FromMutable(v)[0] = v;
          if ( To(v)[0] == w )
          {    ToMutable(v)[0] = u;
               ToEdgeObjMutable(v)[0] = ToEdgeObj(v)[1];    }
          FromEdgeObjMutable(v)[0] = vv_id;
          ToMutable(v).resize(1), ToEdgeObjMutable( )[v].resize(1);
          ToMutable(v).push_back(v);
          ToEdgeObjMutable( )[v].push_back(vv_id);
          FromMutable(v).push_back(x);
          FromEdgeObjMutable( )[v].push_back(vx_id);
          for ( size_t j = 0; j < To(x).size(); j++ )
          {    if ( To(x)[j] == w )
               {    ToMutable(x)[j] = v;
                    ToEdgeObjMutable(x)[j] = vx_id;
                    break;    }    }
          SortSync( FromMutable(v), FromEdgeObjMutable(v) );
          SortSync( ToMutable(v), ToEdgeObjMutable(v) );
          SortSync( FromMutable(w), FromEdgeObjMutable(w) );
          SortSync( ToMutable(w), ToEdgeObjMutable(w) );
          SortSync( ToMutable(x), ToEdgeObjMutable(x) );    }    }



void HyperKmerPath::DumpFasta( const String& fn, const KmerBaseBroker& kbb ) const
{    Ofstream( out, fn );
     vec<int> to_left, to_right;
     ToLeft(to_left), ToRight(to_right);
     for ( int e = 0; e < EdgeObjectCount( ); e++ )
     {    int v = to_left[e], w = to_right[e];
	  out << ">edge_" << BaseAlpha(e) << " " << v << ":" << w << "\n";
          kbb.ToSequence( EdgeObject(e) ).PrintN(out);    }    }

void HyperKmerPath::DumpFastb( const String& fn, const KmerBaseBroker& kbb ) const
{    IncrementalWriter<basevector> basesOut( fn.c_str( ) );
     for ( int e = 0; e < EdgeObjectCount( ); e++ )
     {    basesOut.add( kbb.Seq( EdgeObject(e) ) );    }
     basesOut.close();    }


void HyperKmerPath::SetToDisjointUnionOf( const vec<HyperKmerPath>& v )
{    for ( size_t i = 0; i < v.size(); i++ )
          ForceAssertEq( K( ), v[i].K( ) );
     Clear( );
     int nvert = 0, nedge = 0;
     for ( size_t i = 0; i < v.size(); i++ )
     {    nvert += v[i].N( );
          nedge += v[i].EdgeObjectCount( );    }
     FromMutable( ).reserve(nvert), ToMutable( ).reserve(nvert);
     FromEdgeObjMutable( ).reserve(nvert), ToEdgeObjMutable( ).reserve(nvert);
     EdgesMutable( ).reserve(nedge);
     int vcount = 0, ecount = 0;
     for ( size_t i = 0; i < v.size(); i++ )
     {    const HyperKmerPath& h = v[i];
          EdgesMutable( ).append( h.Edges( ) );
          vec< vec<int> > from = h.From( ), to = h.To( );
          vec< vec<int> > frome = h.FromEdgeObj( ), toe = h.ToEdgeObj( );
          for ( int j = 0; j < h.N( ); j++ )
          {    for ( size_t u = 0; u < from[j].size(); u++ )
               {    from[j][u] += vcount;
                    frome[j][u] += ecount;    }
               for ( size_t u = 0; u < to[j].size(); u++ )
               {    to[j][u] += vcount;
                    toe[j][u] += ecount;    }    }
          FromMutable( ).append(from), ToMutable( ).append(to);
          FromEdgeObjMutable( ).append(frome), ToEdgeObjMutable( ).append(toe);
          vcount += h.N( );
          ecount += h.EdgeObjectCount( );    }    }


void BinaryWrite( int fd, const HyperKmerPath& h )
{    WriteBytes( fd, &h.K_, sizeof(int) );
     BinaryWrite( fd, (const digraphE<KmerPath>&) h );    }

void BinaryRead( int fd, HyperKmerPath& h )
{    ReadBytes( fd, &h.K_, sizeof(int) );
     BinaryRead( fd, (digraphE<KmerPath>&) h );    }

void HyperKmerPath::ReverseComponent( int x )
{    equiv_rel e( N( ) );
     for ( int v = 0; v < N( ); v++ )
     {    for ( size_t i = 0; i < From(v).size(); i++ )
          {    int w = From(v)[i];
               e.Join( v, w );    }    }
     vec<int> o;
     e.Orbit( x, o );
     for ( size_t j = 0; j < o.size(); j++ )
     {    int i = o[j];
          for ( size_t u = 0; u < From(i).size(); u++ )
          {    int e = EdgeObjectIndexByIndexFrom( i, u );
               KmerPath& p = EdgeObjectMutable(e);
               p.Reverse( );    }    }
     digraphE<KmerPath>::ReverseComponent(x);    }

void HyperKmerPath::Reverse( )
{    for ( int i = 0; i < EdgeObjectCount( ); i++ )
          EdgeObjectMutable(i).Reverse( );
     digraphE<KmerPath>::Reverse( );    }

String EdgeLabel( int v, int w, const KmerPath& p, int K )
{    float len = 0, dev = 0;
     for ( int u = 0; u < p.NSegments( ); u++ )
     {    const KmerPathInterval& x = p.Segment(u);
          if ( x.isSeq( ) ) len += x.Length( );
          else
          {    len += float( x.Maximum( ) + x.Minimum( ) ) / 2.0;
               dev += float( x.Maximum( ) - x.Minimum( ) ) / 2.0;    }    }
     ostrstream out;
     out << v << " --- " << setprecision(9) << len << " +/- " << dev
          << " --> " << w << ends;
     return out.str( );    }




     


void HyperKmerPath::TestValid( ) const
{    digraphE<KmerPath>::TestValid( );
     for ( int v = 0; v < N( ); v++ )
     {    for ( size_t j = 0; j < From(v).size(); j++ )
          {    int w = From(v)[j];
               const KmerPath& p = EdgeObjectByIndexFrom( v, j );
               ForceAssertGt( p.NSegments( ), 0 );    }    }    }

void HyperKmerPath::PrintSummary( ostream& out ) const
{    for ( int v = 0; v < N( ); v++ )
     {    for ( size_t j = 0; j < From(v).size(); j++ )
          {    int w = From(v)[j];
               const KmerPath& p = EdgeObjectByIndexFrom( v, j );
               out << EdgeLabel( v, w, p, K( ) ) << "\n";    }    }    }

void HyperKmerPath::PrintSummaryDOT( ostream& out ) const
{    vec< vec<String> > edge_labels;
     edge_labels.resize( N( ) );
     for ( int v = 0; v < N( ); v++ )
     {    const vec<int>& from = From(v);
          edge_labels[v].resize( from.size( ) );
          for ( size_t j = 0; j < from.size(); j++ )
          {    const KmerPath& p = EdgeObjectByIndexFrom( v, j );
               int w = from[j];
               edge_labels[v][j] = EdgeLabel( v, w, p, K( ) );    }    }
     DOT( out, edge_labels );    }

void HyperKmerPath::PrintSummaryDOTAlt( ostream& out ) const
{    vec< vec<String> > edge_labels;
     edge_labels.resize( N( ) );
     for ( int v = 0; v < N( ); v++ )
     {    const vec<int>& from = From(v);
          edge_labels[v].resize( from.size( ) );
          for ( size_t j = 0; j < from.size(); j++ )
          {    int id = EdgeObjectIndexByIndexFrom( v, j );
               edge_labels[v][j] = ToString(id);    }    }
     DOT( out, edge_labels );    }

void HyperKmerPath::PrintSummaryDOT0( ostream& out ) const
{    DOT( out );    }

// Wrapper for the templatized digraphE function PrettyDOT.

void HyperKmerPath::PrintSummaryDOT0w( ostream& out, Bool label_contigs,
     Bool label_vertices, Bool label_edges, const vec<int>* componentsToPrint,
     const Bool edge_labels_base_alpha, const vec<String>* edge_labels_extra,
     const vec<String>* label_contigs_extra, const vec<int> *verticesToPrint ) const
{    vec<double> lengths( EdgeObjectCount( ) );
     for ( int i = 0; i < EdgeObjectCount( ); i++ )
         lengths[i] = EdgeObject(i).MidLength( );
     PrettyDOT( out, lengths, label_contigs, label_vertices, label_edges,
		componentsToPrint, edge_labels_base_alpha, edge_labels_extra,
		label_contigs_extra, verticesToPrint );    }

void HyperKmerPath::PrintSummaryPlus( ostream& out,
     const void * scratch1, KmerBaseBroker* kbb,
     void * scratch2, void * scratch3,
     int scratch4, Bool print_kmer_ranges,
     const vec<String>* edge_remarks, Bool print_component_id_line,
     const vec<String>* component_remarks,
     const vec<Bool>* component_remarks_only ) const
{
     // Define equivalence relation whose orbits are the components.

     equiv_rel e( N( ) );
     for ( int v = 0; v < N( ); v++ )
     {    for ( size_t j = 0; j < From(v).size(); j++ )
               e.Join( v, From(v)[j] );    }

     // Set up indices.
     vec< vec<int> > to_comp, to_comp_pos;

     // Go through the components.  There are two passes.  In the first pass, we
     // just fill in the indices.

     for ( int pass = 1; pass <= 2; pass++ )
     {
     int count = 0;
     for ( int x = 0; x < N( ); x++ )
     {    if ( e.Representative(x) )
          {    if ( pass == 2 && print_component_id_line )
                    out << "\ncomponent " << count << "\n";
               if ( pass == 2 && component_remarks != 0 )
                    out << (*component_remarks)[count];
               if ( component_remarks_only != 0 && (*component_remarks_only)[count] )
               {    ++count;
                    continue;    }
               ++count;
               vec<int> o;
               e.Orbit( x, o );
               Sort(o);
               vec<float> pos( o.size( ) );
               vec<Bool> placed( o.size( ), False );
               pos[0] = 0.0, placed[0] = True;
               while( Sum(placed) < (int)o.size() )
               {    for ( size_t i1 = 0; i1 < o.size(); i1++ )
                    {    int v = o[i1];
                         for ( size_t j = 0; j < From(v).size(); j++ )
                         {    int w = From(v)[j];
                              int i2 = BinPosition( o, w );
                              if ( !( placed[i1] ^ placed[i2] ) ) continue;
                              const KmerPath& p = EdgeObjectByIndexFrom( v, j );
                              if ( placed[i1] ) pos[i2] = pos[i1] + p.MidLength( );
                              else pos[i1] = pos[i2] - p.MidLength( );
                              placed[i1] = placed[i2] = True;    }    }    }
               vec<float> opos( o.size( ) );
               for ( size_t i = 0; i < o.size(); i++ )
                    opos[i] = pos[i];
               SortSync( opos, o );
               int edgeid = 0;
               for ( size_t i = 0; i < o.size(); i++ )
               {    int v = o[i];
                    vec<int> f = From(v);
                    vec<int> ind( f.size( ) );
                    for ( size_t j = 0; j < ind.size(); j++ )
                         ind[j] = j;
                    vec<float> flen( f.size( ) );
                    for ( size_t j = 0; j < f.size(); j++ )
                         flen[j] = EdgeObjectByIndexFrom( v, j ).MidLength( );
                    SortSync( flen, ind );
                    for ( size_t j = 0; j < ind.size(); j++ )
                    {    ++edgeid;
                         int w = f[ ind[j] ];
                         const KmerPath& p = EdgeObjectByIndexFrom( v, ind[j] );

                         // Record component and edge indices.

                         int e = EdgeObjectIndexByIndexFrom( v, ind[j] );
                         if ( pass == 1 ) continue;

                         // Print edge.

                         out << "\n" << "[";
                         if (print_component_id_line)
                              out << count-1 << "." << edgeid-1
                                  << " = " << BaseAlpha(e) << "] ";
                         else out << BaseAlpha(e) << "] ";
                         out << EdgeLabel( v, w, p, K( ) ) << "\n";
                         if (print_kmer_ranges) PrintFolded( out, p );
                         if ( To(v).empty( ) ) out << "[" << v << " is source]\n";
                         if ( From(w).empty( ) ) out << "[" << w << " is sink]\n";

                         // Check for bubble: v != w and there are exactly two edges
                         // from v to w, and both are gap-free.  Then print the
                         // alignment of the two edges.

                         vec<int> js;
                         for ( size_t j2 = 0; j2 < ind.size(); j2++ )
                         {    int w2 = f[ ind[j2] ];
                              if ( w2 == w ) js.push_back(j2);    }
                         if ( v != w && js.size() == 2 && js[0] == (int)j && kbb != 0 )
                         {    int j2 = js[1];
                              const KmerPath& p2
                                   = EdgeObjectByIndexFrom( v, ind[j2] );
                              if ( p.GapFree( ) && p2.GapFree( ) )
                              {    basevector b = kbb->ToSequence(p).Seq(0);
                                   basevector b2 = kbb->ToSequence(p2).Seq(0);
                                   align a;
                                   int offset = 0, bandwidth = 20, errors;
                                   SmithWatBandedA( b, b2, offset,
                                        bandwidth, a, errors );
                                   out << "\nalignment of the two paths from "
                                        << v << " to " << w << "\n";
                                   PrintVisualAlignment(
                                        True, out, b, b2, a );    }    }

                         if ( edge_remarks != 0 ) out << (*edge_remarks)[e];
		         }    }    }    }    }
     flush(out);    }

void HyperKmerPath::PrintSummaryPlusAlt( ostream& out,
     const vec<String>* edge_remarks ) const
{
     // Compute the components.

     vec< vec<int> > comps;
     Components(comps);

     // Go through the components.  There are two passes.  In the first pass, we
     // just fill in the indices.

     for ( int pass = 1; pass <= 2; pass++ )
     {    for ( size_t ic = 0; ic < comps.size(); ic++ )
          {    if ( pass == 2 ) out << "\ncomponent " << ic << "\n";
               vec<int> o = comps[ic];
               vec<float> pos( o.size( ) );
               vec<Bool> placed( o.size( ), False );
               pos[0] = 0.0, placed[0] = True;
               while( Sum(placed) < (int)o.size() )
               {    for ( size_t i1 = 0; i1 < o.size(); i1++ )
                    {    int v = o[i1];
                         for ( size_t j = 0; j < From(v).size(); j++ )
                         {    int w = From(v)[j];
                              int i2 = BinPosition( o, w );
                              if ( !( placed[i1] ^ placed[i2] ) ) continue;
                              const KmerPath& p = EdgeObjectByIndexFrom( v, j );
                              if ( placed[i1] ) pos[i2] = pos[i1] + p.MidLength( );
                              else pos[i1] = pos[i2] - p.MidLength( );
                              placed[i1] = placed[i2] = True;    }    }    }
               vec<float> opos( o.size( ) );
               for ( size_t i = 0; i < o.size(); i++ )
                    opos[i] = pos[i];
               SortSync( opos, o );
               int edgeid = 0;
               for ( size_t i = 0; i < o.size(); i++ )
               {    int v = o[i];
                    vec<int> f = From(v);
                    vec<int> ind( f.size( ) );
                    for ( size_t j = 0; j < ind.size(); j++ )
                         ind[j] = j;
                    vec<float> flen( f.size( ) );
                    for ( size_t j = 0; j < f.size(); j++ )
                         flen[j] = EdgeObjectByIndexFrom( v, j ).MidLength( );
                    SortSync( flen, ind );
                    for ( size_t j = 0; j < ind.size(); j++ )
                    {    ++edgeid;
                         int w = f[ ind[j] ];
                         const KmerPath& p = EdgeObjectByIndexFrom( v, ind[j] );

                         // Record component and edge indices.

                         int e = EdgeObjectIndexByIndexFrom( v, ind[j] );
                         if ( pass == 1 ) continue;

                         // Print edge.

                         out << "\n" << "[" << ic << "." << edgeid-1
                             << " = " << BaseAlpha(e) << "] ";
                         out << EdgeLabel( v, w, p, K( ) ) << "\n";
                         if ( To(v).empty( ) ) out << "[" << v << " is source]\n";
                         if ( From(w).empty( ) )
                              out << "[" << w << " is sink]\n";
                         if ( edge_remarks != 0 )
                              out << (*edge_remarks)[e];    }    }    }    }
     flush(out);    }

void HyperKmerPath::ContractEmptyEdges( ) {
  bool contracted = false;
  for( int v = 0; v < N(); v += (!contracted) ) {
    contracted = false;
    for( int j = From(v).size()-1; j>=0; j-- )
      if( EdgeObjectByIndexFrom(v, j).IsEmpty() ) {
	ContractEdgeFrom(v,j);
	contracted = true;
	break;
      }
  }
}




void HyperKmerPath::CompressEdgeObjects( )
{    for ( int i = 0; i < EdgeObjectCount( ); i++ )
     {    KmerPath& p = EdgeObjectMutable(i);
          KmerPath pnew;
	  pnew.Reserve( p.NSegments( ) );
          for ( int j = 0; j < p.NSegments( ); j++ )
               pnew.AddSegment( p.Segment(j) );
          p = pnew;    }    }







HyperKmerPath::HyperKmerPath( const String& filename )
{    int fd = OpenForRead(filename);
     BinaryRead( fd, *this );
     close(fd);    }

void HyperKmerPath::RemoveSmallComponents( int min_kmers )
{    vec< vec<int> > comps;
     Components(comps);
     vec<int> o, keep;
     for ( size_t i = 0; i < comps.size(); i++ )
     {    const vec<int>& o = comps[i];
          int nkmers = 0;
          for ( size_t j = 0; j < o.size(); j++ )
          {    int v = o[j];
               for ( size_t t = 0; t < From(v).size(); t++ )
                    nkmers += EdgeObjectByIndexFrom( v, t ).KmerCount( );    }
          if ( nkmers < min_kmers ) continue;
          keep.append(o);    }
     HyperKmerPath h( *this, keep );
     *this = h;    }


// Remove the loop subgraph from this HKP, leaivng only non-loop edges.
void HyperKmerPath::RemoveLoopSubgraph()
{
  
  // Define reverse complementation.
  vec<int> h_to_rcx( EdgeObjectCount() );
  vec<KmerPath> edges;
  vec<int> ids( EdgeObjectCount(), vec<int>::IDENTITY );
  for ( int i = 0; i < EdgeObjectCount(); i++ ) {
    KmerPath p = EdgeObject(i), q;
    p.Canonicalize( );
    edges.push_back(p);
  }
  SortSync( edges, ids );
  for ( int i = 0; i < EdgeObjectCount(); i++ ) {
    KmerPath q = EdgeObject(i);
    q.Reverse( );
    q.Canonicalize( );
    h_to_rcx[i] = ids[ BinPosition( edges, q ) ];
  }
  
  // Find the loop subgraph.
  vec<int> loop_edges;
  LoopSubgraph(loop_edges);
  cout << Date( ) << ": deleting loop subgraph" << endl;
  vec<int> loops_to_delete;
  for ( size_t i = 0; i < loop_edges.size(); i++ )
    loops_to_delete.push_back( loop_edges[i] );
  size_t nloops = loops_to_delete.size( );
  for ( size_t i = 0; i < nloops; i++ )
    loops_to_delete.push_back( h_to_rcx[ loops_to_delete[i] ] );
  UniqueSort(loops_to_delete);
  
  // Delete the loop subgraph.
  DeleteEdges(loops_to_delete);
}





// Remove the loop subgraph, and also remove branches as much as possible.
// The resulting HyperKmerPath should be clean and acyclic.
void
HyperKmerPath::MakeAcyclic( )
{
  // Identify loop subgraph.  
  
  vec<int> loop_edges, to_split;
  LoopSubgraph(loop_edges);

  // Delete loop edges, and for every vertex incident
  // upon a loop edge, separate the residual in edges from the out edges.

  for ( int v = 0; v < N( ); v++ )
    for ( int j = 0; j < From(v).isize( ); j++ ) {
      int e = EdgeObjectIndexByIndexFrom( v, j );
      if ( BinMember( loop_edges, e ) )
	to_split.push_back( v, From(v)[j] );
    }
  
  UniqueSort(to_split);
  DeleteEdges(loop_edges);
  for ( int i = 0; i < to_split.isize( ); i++ ) {
    int v = to_split[i];
    if ( From(v).nonempty( ) && To(v).nonempty( ) ) {
      AddVertices(1);
      TransferEdges( v, N( ) - 1, True );
    }
  }

  // Identify branch edges.  These are edges having two completely separate
  // sequences exiting, each of length at least 1 kb.  Or the other direction.
  // There is one pass for each direction.  The comments below apply to pass 1.

  vec<int> branches;
  int min_len = 1000;
  for ( int pass = 1; pass <= 2; pass++ ) {
    Bool fw = ( pass == 1 );
    vec<int> dist_to_end;
    DistancesToEnd( *this, &KmerPath::KmerCount, min_len, fw, dist_to_end );
    for ( int v = 0; v < N( ); v++ ) {
      const vec<int>& FTv = ( fw ? From(v) : To(v) );
      for ( int j = 0; j < FTv.isize( ); j++ ) {
	int w = FTv[j]; 
	int e = ( fw ? EdgeObjectIndexByIndexFrom( v, j )
		  : EdgeObjectIndexByIndexTo  ( v, j ) );
	
	// Now we have an edge e from v to w.  Define an equivalence 
	// relation on the edges the emanate from w.
	
	const vec<int>& FTw = ( fw ? From(w) : To(w) );
	int n = FTw.size( );
	if ( n < 2 ) continue; // optimization
	if ( n == 2 && FTw[0] == FTw[1] ) continue; // optimization
	equiv_rel E(n);
	vec< vec<int> > suc(n);
	for ( int i = 0; i < n; i++ )
	  if (fw) GetSuccessors1( FTw[i], suc[i] );
	  else GetPredecessors1( FTw[i], suc[i] );
	for ( int i1 = 0; i1 < n; i1++ )
	  for ( int i2 = i1+1; i2 < n; i2++ )
	    if ( Meet( suc[i1], suc[i2] ) )
	      E.Join( i1, i2 );
	
	vec<int> reps;
	E.OrbitRepsAlt(reps);
	if ( reps.size( ) < 2 ) continue; // optimization
	
	// Compute exit length.
	
	vec<int> dist( reps.size( ), 0 );
	for ( int i = 0; i < n; i++ ) {
	  int r = BinPosition( reps, E.ClassId(i) );
	  dist[r] = Max( dist[r],
			 ( fw ?  EdgeObjectByIndexFrom( w, i )
			   : EdgeObjectByIndexTo  ( w, i ) ).KmerCount( )
			 + dist_to_end[ FTw[i] ] );
	}
	ReverseSort(dist);
	if ( dist[1] >= min_len ) branches.push_back(e);
      }
    }
  }
  UniqueSort(branches);

     // Remove the branch edges from the graph and spread out any edges that 
     // touch them.

     vec<int> to_left, to_right;
     ToLeft(to_left), ToRight(to_right);
     int N0 = N( );
     for ( int i = 0; i < branches.isize( ); i++ ) 
     {
          // First find the edge.  Its endpoints may have moved.
    
          int v = to_left[ branches[i] ], w = to_right[ branches[i] ];
          int l = EdgeObjectIndexToFromIndex( v, branches[i] );
          int r = EdgeObjectIndexToToIndex  ( w, branches[i] );
          if ( l < 0 && r < 0 ) continue; // already liberated
          if ( l < 0 )
          {    for ( int x = N0; x < N( ); x++ )
	       {    for ( int j = 0; j < From(x).isize( ); j++ )
	            {    if ( EdgeObjectIndexByIndexFrom( x, j ) == branches[i] ) 
                         {    v = x;
	                      goto found_left;    }    }    }    }
          found_left:
          if ( r < 0 )
          {    for ( int x = N0; x < N( ); x++ )
	       {    for ( int j = 0; j < To(x).isize( ); j++ )
	            {    if ( EdgeObjectIndexByIndexTo( x, j ) == branches[i] ) 
                         {    w = x;
	                      goto found_right;    }    }    }    }
          found_right:
    
          // Now liberate the edge.
    
          LiberateEdge( branches[i], v, w );    }

     // Clean up the graph.

  RemoveEdgelessVertices();
  RemoveUnneededVertices();
}



void HyperKmerPath::MakeEdgeDatabase( vec<tagged_rpint>& edgedb ) const
{    edgedb.clear( );
     edgedb.reserve( EdgeObjectCount( ) );
     for ( int i = 0; i < EdgeObjectCount( ); i++ )
          EdgeObject(i).AppendToDatabase( edgedb, i );
     Prepare(edgedb);    }



bool ComponentsAreIsomorphic( const HyperKmerPath& hkp1, int seed_vx_1,
			      const HyperKmerPath& hkp2, int seed_vx_2,
			      vec<int>* p_iso ) {
  vec<int> local_iso;
  if( p_iso == NULL ) p_iso = &local_iso;

  p_iso->resize( hkp1.N(), -1 );
  // (*p_iso)[i]=-1 means hkp1 vx i unused; j>=0 means matches hkp2 vx j.
  vec<Bool> used2( hkp2.N(), false );

  // The hypothesis <i,j> means vertices i and j ought to match
  vec< pair<int,int> > hypotheses;
  // The total number of hypotheses will be two per edge in the
  // connected component (but generally not all on the stack at once).
  hypotheses.reserve( 2 * hkp1.EdgeObjectCount() );
  // Initial hypothesis: the given seed vertices match.
  hypotheses.push_back( make_pair(seed_vx_1,seed_vx_2) );

  while( ! hypotheses.empty() ) {

    int v1 = hypotheses.back().first, v2=hypotheses.back().second;
    hypotheses.pop_back();

    if( (*p_iso)[v1]==v2 )                     // already knew this hypothesis
      continue;
    else if( used2[v2] || (*p_iso)[v1]!=-1 )   // knew a contradictory hypothesis
      return false;

    // Otherwise:
    // (1) Assume this hypothesis is true
    (*p_iso)[v1] = v2;
    used2[v2] = true;

    // (2) check that the edges in/out of these vertices match up
    vec<int> from_v1 = hkp1.From(v1), to_v1 = hkp1.To(v1);
    vec<int> from_v2 = hkp2.From(v2), to_v2 = hkp2.To(v2);

    if( from_v1.size() != from_v2.size() || to_v1.size() != to_v2.size() )
      return false;

    vec<KmerPath> from_v1_edges( from_v1.size() );
    vec<KmerPath> from_v2_edges( from_v2.size() );
    vec<KmerPath> to_v1_edges( to_v1.size() );
    vec<KmerPath> to_v2_edges( to_v2.size() );

    for(size_t i=0; i<from_v1.size(); i++)
      from_v1_edges[i] = hkp1.EdgeObjectByIndexFrom(v1,i);
    for(size_t i=0; i<from_v2.size(); i++)
      from_v2_edges[i] = hkp2.EdgeObjectByIndexFrom(v2,i);

    for(size_t j=0; j<to_v1.size(); j++)
      to_v1_edges[j] = hkp1.EdgeObjectByIndexTo(v1,j);
    for(size_t j=0; j<to_v2.size(); j++)
      to_v2_edges[j] = hkp2.EdgeObjectByIndexTo(v2,j);

    SortSync( from_v1_edges, from_v1 );
    SortSync( from_v2_edges, from_v2 );
    SortSync( to_v1_edges, to_v1 );
    SortSync( to_v2_edges, to_v2 );

    if( from_v1_edges != from_v2_edges || to_v2_edges != to_v2_edges )
      return false;

    ForceAssert( from_v1_edges.UniqueOrdered() );
    ForceAssert( from_v2_edges.UniqueOrdered() );
    ForceAssert( to_v1_edges.UniqueOrdered() );
    ForceAssert( to_v2_edges.UniqueOrdered() );

    for( size_t i=0; i < from_v1.size(); i++ )
      hypotheses.push_back( make_pair(from_v1[i],from_v2[i]) );
    for( size_t j=0; j < to_v1.size(); j++ )
      hypotheses.push_back( make_pair(to_v1[j],to_v2[j]) );


    // (3) formulate all the corollary hypotheses imposed by edge matching
    for( size_t i=0; i < from_v1.size(); i++ )
      hypotheses.push_back( make_pair(from_v1[i],from_v2[i]) );
    for( size_t j=0; j < to_v1.size(); j++ )
      hypotheses.push_back( make_pair(to_v1[j],to_v2[j]) );

  } // end of processing of hypotheses.back()

  // If no hypotheses gave rise to a contradiction,
  // then we have an isomorphism.  Hooray!
  return true;
}

void HyperKmerPath::Zipper( )
{
     // Iterate until no further edits are possible.

     while(1)
     {    int edits = 0;

          // Go through two passes, one for each orientation.

          for ( int pass = 1; pass <= 2; pass++ )
          {
               // Go through each vertex v.

               for ( int v = 0; v < N( ); v++ )
               {    restart:

                    // Go through the edges e1 exiting v.  Edges from v to v are
                    // ignored.

                    for ( size_t j1 = 0; j1 < From(v).size(); j1++ )
                    {    int w1 = From(v)[j1];
                         if ( w1 == v ) continue;
                         const KmerPath& e1 = EdgeObjectByIndexFrom( v, j1 );

                         // Now go through the other edges e2 exiting v.  Again
                         // edges from v to v are ignored.

                         for ( size_t j2 = j1+1; j2 < From(v).size(); j2++ )
                         {    int w2 = From(v)[j2];
                              if ( w2 == v ) continue;
                              const KmerPath& e2 = EdgeObjectByIndexFrom( v, j2 );

                              // Find the maximum perfect match between e1 and e2,
                              // starting at the beginning.  A match of at least one
                              // kmer is required.

                              if ( e1.Start(0) != e2.Start(0) ) continue;
                              KmerPathLoc loc1(e1, 0), loc2(e2, 0);
                              ScanRightPerfectMatch( loc1, loc2 );

                              // Handle the case where the two edges are the same.

                              if ( loc1.atEnd( ) && loc2.atEnd( ) )
                              {    DeleteEdgeFrom( v, j2 );
                                   if ( w1 != w2 ) TransferEdges( w1, w2 );
                                   ++edits;
                                   goto restart;    }

                              // Handle the case where the match goes all the way
                              // to the end on just one edge.

                              if ( loc1.atEnd( ) )
                              {    KmerPath before, after;
                                   e2.CopySubpath( e2.Begin( ), loc2, before );
                                   e2.CopySubpathNoFirstKmer(loc2, e2.End( ), after);
                                   DeleteEdgeFrom( v, j2 ); DeleteEdgeFrom( v, j1 );
                                   int x = N( );
                                   AddVertices(1);
                                   AddEdge( v, x, before );
                                   AddEdge( x, w2, after );
                                   TransferEdges( w1, x );
                                   ++edits;
                                   goto restart;    }
                              if ( loc2.atEnd( ) )
                              {    KmerPath before, after;
                                   e1.CopySubpath( e1.Begin( ), loc1, before );
                                   e1.CopySubpathNoFirstKmer(loc1, e1.End( ), after);
                                   DeleteEdgeFrom( v, j2 ); DeleteEdgeFrom( v, j1 );
                                   int x = N( );
                                   AddVertices(1);
                                   AddEdge( v, x, before );
                                   AddEdge( x, w1, after );
                                   TransferEdges( w2, x );
                                   ++edits;
                                   goto restart;    }

                              // Join the two edges up to the point where they
                              // disagree.

                              ++edits;
                              KmerPath before, after1, after2;
                              e1.CopySubpath( e1.Begin( ), loc1, before );
                              e1.CopySubpathNoFirstKmer( loc1, e1.End( ), after1 );
                              e2.CopySubpathNoFirstKmer( loc2, e2.End( ), after2 );
                              DeleteEdgeFrom( v, j2 );
                              DeleteEdgeFrom( v, j1 );
                              int x = N( );
                              AddVertices(1);
                              AddEdge( v, x, before );
                              AddEdge( x, w1, after1 );
                              AddEdge( x, w2, after2 );
                              goto restart;    }    }    }
               Reverse( );    }
         if ( edits == 0 ) break;
         ReduceLoops( );
         CompressEdgeObjects( );
	 RemoveDuplicateEdges( ); // the zippering sometimes leaves two vertices with multiple copies of the exact same edge between then
	 RemoveUnneededVertices( );
         RemoveDeadEdgeObjects( );
         RemoveEdgelessVertices( );    }    }

// MethodDecl: CanonicalizeEdges
// Canonicalize the KmerPaths on all the edges.
void HyperKmerPath::CanonicalizeEdges() {
  vec<KmerPath>& edges = EdgesMutable();
  for ( size_t i = 0; i < edges.size(); i++ )
    edges[i].Canonicalize();
}


/**
   Method: FindIsomorphicComponents

   Test whether the HyperKmerPath has two isomorphic connected components.

   Return the equivalence relation of the components, and the list of ids of those components
   that have an isomorphic partner in the graph.

*/
Bool HyperKmerPath::FindIsomorphicComponents( equiv_rel& componentRelation,
					      vec<int>& isomorphicComponentReps,
					      Bool stopIfFoundOne ) const {
  // Determine the connected components
  componentRelation.Initialize( N() );
  ComponentRelation( componentRelation );

  // For each connected component compute a "signature" -- a vector of
  // integers such that isomorphic components will have equal signatures.
  // We'll then test components with equal signatures for isomorphism.

  // A component's signature is a vector of integers that looks like this:
  // #vertices #edges sorted-list-of-vertex-out-degrees sorted-list-of-vertex-in-degrees sorted-list-of-edge-lengths

  typedef vec<int> component_signature_t;

  vec< component_signature_t > componentSignatures;

  vec<int> componentReps;
  componentRelation.OrbitRepsAlt( componentReps );
  for ( size_t i = 0 ; i < componentReps.size(); i++ ) {
    int thisComponentRep = componentReps[i];
    vec<int> thisComponentVerts;
    componentRelation.Orbit( thisComponentRep, thisComponentVerts );
    component_signature_t thisComponentSig;
    thisComponentSig.push_back( thisComponentVerts.size() );

    int thisComponentNumEdges = 0;
    for ( size_t j = 0; j < thisComponentVerts.size(); j++ )
      thisComponentNumEdges += From(thisComponentVerts[j]).size();

    thisComponentSig.push_back( thisComponentNumEdges );

    vec<int> thisComponentOutDegrees;
    for ( size_t j = 0; j < thisComponentVerts.size(); j++ )
      thisComponentOutDegrees.push_back( From(thisComponentVerts[j]).size() );
    Sort( thisComponentOutDegrees );
    thisComponentSig.append( thisComponentOutDegrees );

    vec<int> thisComponentInDegrees;
    for ( size_t j = 0; j < thisComponentVerts.size(); j++ )
      thisComponentInDegrees.push_back( To(thisComponentVerts[j]).size() );
    Sort( thisComponentInDegrees );
    thisComponentSig.append( thisComponentInDegrees );

    vec<int> thisComponentEdgeLens;
    for ( size_t j = 0; j < thisComponentVerts.size(); j++ ) {
      const vec<int>& fromThisVrtx = FromEdgeObj(thisComponentVerts[j]);
      for ( size_t k = 0; k < fromThisVrtx.size(); k++ )
	thisComponentEdgeLens.push_back( EdgeObject( fromThisVrtx[k] ).KmerCount() );
    }

    Sort( thisComponentEdgeLens );
    thisComponentSig.append( thisComponentEdgeLens );

    componentSignatures.push_back( thisComponentSig );
  }

  SortSync( componentSignatures, componentReps );

  Bool foundIsomorphicComponents = False;
  isomorphicComponentReps.clear();
  for ( size_t i = 0; i+1 < componentSignatures.size(); i++ ) {
    if ( componentSignatures[i] == componentSignatures[i+1] &&
	 ComponentsAreIsomorphic( *this, componentReps[i], *this, componentReps[i+1] ) ) {
      foundIsomorphicComponents = True;
      isomorphicComponentReps.push_back( componentReps[i] );
      isomorphicComponentReps.push_back( componentReps[i+1] );
      if ( stopIfFoundOne )
	return True;
    }
  }
  UniqueSort( isomorphicComponentReps );

  return foundIsomorphicComponents;
}

int HyperKmerPath::EdgeN50( ) const
{    vec<int> X( EdgeObjectCount( ) );
     for ( size_t i = 0; i < X.size(); i++ )
          X[i] = EdgeLength(i);
     Sort(X);
     return ( X.nonempty( ) ? N50(X) : 0 );    }

     // If two components of the graph are reverse complements of each other,
     // delete one of them.

void HyperKmerPath::DeleteReverseComplementComponents( )
{    vec<int> to_rc( EdgeObjectCount( ), -1 );
     {    vec<KmerPath> edges;
          vec<int> ids( EdgeObjectCount( ), vec<int>::IDENTITY );
          for ( int i = 0; i < EdgeObjectCount( ); i++ )
          {    KmerPath p = EdgeObject(i), q;
               p.Canonicalize( );
               edges.push_back(p);    }
          SortSync( edges, ids );
          for ( int i = 0; i < EdgeObjectCount( ); i++ )
          {    KmerPath q = EdgeObject(i);
               q.Reverse( );
               q.Canonicalize( );
               int pos = BinPosition( edges, q );
               if ( pos >= 0 ) to_rc[i] = ids[pos];    }    }
     vec< vec<int> > components;
     ComponentEdges(components);
     Sort(components);
     vec<int> rc_to_delete;
     for ( size_t i = 0; i < components.size(); i++ )
     {    vec<int> rc;
          for ( size_t j = 0; j < components[i].size(); j++ )
               rc.push_back( to_rc[ components[i][j] ] );
          Sort(rc);
          int p = BinPosition( components, rc );
          if ( p > (int)i )
          {    for ( size_t j = 0; j < components[p].size(); j++ )
                    rc_to_delete.push_back( components[p][j] );    }    }
     Sort(rc_to_delete);
     DeleteEdges(rc_to_delete);    }
