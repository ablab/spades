///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include <strstream>

#include "CoreTools.h"
#include "Equiv.h"
#include "feudal/IncrementalWriter.h"
#include "graph/DigraphTemplate.h"
#include "paths/HyperBasevector.h"
#include "paths/KmerBaseBroker.h"

// Template instantiations:

template digraphE<basevector>::digraphE(vec<basevector> const&, const ConstructorBehavior );
template void digraphE<basevector>::Clear( );
template void digraphE<basevector>::ComponentEdges( vec< vec<int> >& edges ) const;
template void digraphE<basevector>::DeleteEdgesAtVertex(int);
template void digraphE<basevector>::RemoveEdgelessVertices( );
template void digraphE<basevector>::RemoveDeadEdgeObjects();
template void digraphE<basevector>::ToLeft(vec<int>&) const;
template void digraphE<basevector>::ToRight(vec<int>&) const;
template void digraphE<basevector>::DumpGraphML( const String& ) const;
template void digraphE<basevector>::PrettyDOT( ostream& out, const vec<double>& lengths, Bool label_contigs, Bool label_vertices, Bool label_edges, const vec<int>* componentsToPrint, const Bool, const vec<String> * , const vec<String> *, const vec<int> * ) const;




void HyperBasevector::SetToDisjointUnionOf( const vec<HyperBasevector>& v )
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
     {    const HyperBasevector& h = v[i];
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

void HyperBasevector::Reverse( )
{    for ( int i = 0; i < EdgeObjectCount( ); i++ )
          EdgeObjectMutable(i).ReverseComplement( );
     digraphE<basevector>::Reverse( );    }

void HyperBasevector::PrintSummaryDOT0w( ostream& out, Bool label_contigs,
     Bool label_vertices, Bool label_edges, const vec<int>* componentsToPrint,
     const Bool edge_labels_base_alpha, const vec<String> *label_edges_extra,
     const vec<String> *label_contigs_extra, const vec<int>* verticesToPrint ) const
{    vec<double> lengths( EdgeObjectCount( ) );
     for ( int i = 0; i < EdgeObjectCount( ); i++ )
         lengths[i] = EdgeObject(i).size( );
     PrettyDOT( out, lengths, label_contigs, label_vertices, label_edges,
		componentsToPrint, edge_labels_base_alpha,
		label_edges_extra, label_contigs_extra, verticesToPrint ); }

vec<basevector> Convert( const vec<KmerPath>& in, const KmerBaseBroker& kbb )
{    vec<basevector> out;
     for ( size_t j = 0; j < in.size(); j++ )
     {    SuperBaseVector s = kbb.ToSequence( in[j] );
          ForceAssertEq( s.size( ), 1 );
          out.push_back( s.Seq(0) );    }
     return out;    }

HyperBasevector::HyperBasevector( const String& filename )
{    int fd = OpenForRead(filename);
     BinaryRead( fd, *this );
     close(fd);    }

HyperBasevector::HyperBasevector( const HyperKmerPath& h, const KmerBaseBroker& kbb )
     : digraphE<basevector>( h.From( ), h.To( ), Convert( h.Edges( ), kbb ),
			     h.ToEdgeObj( ), h.FromEdgeObj( ) ),
       K_( h.K() )
{
  TestValid();
}

void HyperBasevector::Initialize( const HyperKmerPath& h, const KmerBaseBroker& kbb )
{    (digraphE<basevector>&)(*this) =
          digraphE<basevector>( h.From( ), h.To( ), Convert( h.Edges( ), kbb ),
          h.ToEdgeObj( ), h.FromEdgeObj( ) );
     K_ = h.K( );
     TestValid();
}

void BinaryWrite( int fd, const HyperBasevector& h )
{    WriteBytes( fd, &h.K_, sizeof(int) );
     BinaryWrite( fd, (const digraphE<basevector>&) h );    }

void BinaryRead( int fd, HyperBasevector& h )
{    ReadBytes( fd, &h.K_, sizeof(int) );
     BinaryRead( fd, (digraphE<basevector>&) h );    }

size_t HyperBasevector::writeBinary( BinaryWriter& writer ) const
{
    size_t len = writer.write(K_);
    return len+writer.write(static_cast<digraphE<basevector> const&>(*this));
}

void HyperBasevector::readBinary( BinaryReader& reader )
{
    reader.read(&K_);
    reader.read(static_cast<digraphE<basevector>*>(this));
}

size_t HyperBasevector::externalSizeof()
{
    return 0;
}

// Check that the edge adjacencies in this HyperBasevector make sense.
// If they don't, this causes a FatalErr.
void
HyperBasevector::TestValid( ) const
{
  const String err = "This HyperBasevector is internally inconsistent.  It might have been created from a HyperKmerPath and a KmerBaseBroker that did not match.";
  
  // For every vertex, determine the (K-1)-mer adjacency implied by that vertex.
  for ( int v = 0; v < N(); v++ ) {
    const vec<int> & to = ToEdgeObj(v), from = FromEdgeObj(v);
    if ( to.empty() && from.empty() ) continue;
    basevector adj;
    if ( to.empty() ) adj.SetToSubOf( EdgeObject( from[0] ), 0, K_-1 );
    else adj.SetToSubOf( EdgeObject( to[0] ), EdgeObject( to[0] ).size() - K_ + 1, K_-1 );
    
    // Verify that this (K-1)-mer appears in every edge to/from this vertex.
    for ( size_t i = 0; i < to.size(); i++ )
      if ( adj != basevector( EdgeObject( to[i] ), EdgeObject( to[i] ).size() - K_ + 1, K_-1 ) )
	FatalErr( err );
    
    for ( size_t i = 0; i < from.size(); i++ )
      if ( adj != basevector( EdgeObject( from[i] ), 0, K_-1 ) )
	FatalErr( err );
  }
}


void HyperBasevector::RemoveSmallComponents( int min_kmers )
{    vec< vec<int> > comps;
     Components(comps);
     vec<int> o, keep;
     for ( size_t i = 0; i < comps.size(); i++ )
     {    const vec<int>& o = comps[i];
          int nkmers = 0;
          for ( size_t j = 0; j < o.size(); j++ )
          {    int v = o[j];
               for ( size_t t = 0; t < From(v).size(); t++ )
                    nkmers += EdgeObjectByIndexFrom( v, t ).size() - K( ) + 1;    }
          if ( nkmers < min_kmers ) continue;
          keep.append(o);    }
     HyperBasevector h( *this, keep );
     *this = h;    }

void HyperBasevector::LowerK( int newK )
{    ForceAssertLe( newK, K( ) );
     for ( int i = 0; i < EdgeObjectCount( ); i++ )
     {    ForceAssertGt( EdgeObject(i).size( ), 0u );
          EdgeObjectMutable(i).resize( EdgeObject(i).size( ) + newK - K( ) );    }
     K_ = newK;    }

void HyperBasevector::GenerateVecbasevector( vecbvec &bases ) const
{
  bases.clear( );

  // Reserve memory.
  longlong n_reads = 0;
  longlong n_bases = 0;
  vec< vec<int> > components;
  this->ComponentEdges( components );
  for (size_t ii=0; ii<components.size(); ii++) {
    n_reads += components[ii].size();
    for (size_t jj=0; jj<components[ii].size(); jj++)
      n_bases += this->EdgeObject( components[ii][jj] ).size( );
  }
  bases.Reserve( n_bases / 16 + n_reads, n_reads );

  // Fill bases.
  for (size_t ii=0; ii<components.size(); ii++)
    for (size_t jj=0; jj<components[ii].size(); jj++)
      bases.push_back( this->EdgeObject( components[ii][jj] ) );
}

void
HyperBasevector::DumpFasta( const String& fn ) const
{
  Ofstream( out, fn );
  vec<int> to_left, to_right;
  ToLeft(to_left), ToRight(to_right);
  for ( int e = 0; e < EdgeObjectCount( ); e++ ) {
    int v = to_left[e], w = to_right[e];
    out << ">edge_" << BaseAlpha(e) << " " << v << ":" << w << "\n";
    EdgeObject(e).Print(out);
  }
}

void
HyperBasevector::DumpFastb( const String& fn ) const
{
  IncrementalWriter<basevector> basesOut( fn.c_str( ) );
  for ( int e = 0; e < EdgeObjectCount( ); e++ ) {
    basesOut.add( EdgeObject(e) );
  }
  basesOut.close();
}

Bool operator==( const HyperBasevector& h1, const HyperBasevector& h2 )
{    if ( h1.K( ) != h2.K( ) ) return False;
     return (const digraphE<basevector>&) h1
          == (const digraphE<basevector>&) h2;    }

void HyperBasevector::RemoveUnneededVertices( )
{    for ( int i = 0; i < N( ); i++ )
     {    if ( From(i).size( ) == 1 && To(i).size( ) == 1 && From(i)[0] != i )
          {    basevector p = EdgeObjectByIndexTo( i, 0 );
               p.resize( p.size( ) - K( ) + 1 );
               p.append( EdgeObjectByIndexFrom( i, 0 ) );
               JoinEdges( i, p );    }    }
     RemoveEdgelessVertices( );    }

basevector HyperBasevector::EdgePathToBases( const vec<int>& e )
{    ForceAssert( e.nonempty( ) );
     basevector b = EdgeObject( e[0] );
     for ( int j = 1; j < e.isize( ); j++ )
     {    b.resize( b.isize( ) - ( K( ) - 1 ) );
          b = Cat( b, EdgeObject( e[j] ) );    }
     return b;    }
