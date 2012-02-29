/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "paths/SeqOnHyper.h"

#include "Basevector.h"
#include "CoreTools.h"
#include "PrintAlignment.h"
#include "paths/HyperBasevector.h"
#include "paths/HyperKmerPath.h"


void
SeqOnHyper::TestValid( const HyperBasevector& hb, const vecbasevector & reads,
		       const vecbasevector & edges,
		       const vec<int>& to_left, const vec<int>& to_right ) const
{    ForceAssertGe( Id1( ), 0 );
     ForceAssertGt( Len1( ), 0 );
     ForceAssert( parts_.nonempty( ) );
     ForceAssertEq( parts_.front( ).pos1( ), 0 );
     int K = hb.K( );
     for ( int i = 1; i < parts_.isize( ); i++ )
     {    ForceAssertEq( Pos1(i-1) - pos1(i), K - 1 );
          ForceAssertEq( to_right[ Id2(i-1) ], to_left[ Id2(i) ] );    }
     for ( int i = 1; i < parts_.isize( ) - 1; i++ )
     {    ForceAssertEq( pos2(i), 0 );
          ForceAssertEq( Pos2(i), hb.EdgeLength( Id2(i) ) );    }
     if ( parts_.size( ) > 1 )
     {    ForceAssertEq( Pos2(0), hb.EdgeLength( Id2(0) ) );
          ForceAssertEq( parts_.back( ).pos2( ), 0 );    }
     
     // Compare the basevector of the read sequence and the basevector of
     // the HyperBasevector where it is supposedly aligned.
     for ( int i = 0; i < parts_.isize(); i++ ) {
       const PartialSeqOnHyper& p = parts_[i];
       basevector r = reads[ Id1( ) ];
       if ( Rc1() ) r.ReverseComplement( );
       basevector bv1( r,                p.pos1(), p.Len() );
       basevector bv2( edges[ p.Id2() ], p.pos2(), p.Len() );
       ForceAssert( bv1 == bv2 );
     }
}

void SeqOnHyper::Print( ostream& out, const vec<int>& to_left, 
     const vec<int>& to_right ) const
{    out << Id1( ) << ( Rc1( ) ? "rc" : "fw" );
     for ( int i = 0; i < parts_.isize( ); i++ )
     {    const PartialSeqOnHyper& p = parts_[i];
          longlong e = p.Id2( );
          cout << " " << e << "[" << to_left[e] << "->" << to_right[e] << "]."
               << p.pos2( ) << "-" << p.Pos2( );    }
     cout << "\n";    }

void SeqOnHyper::PrintVisual( ostream& out, const vecbasevector& reads,
			      const vecbasevector& edges) const
{    basevector r = reads[Id1()];
     r.Print(cout);
     if (Rc1())
       r.ReverseComplement();
     r.Print(cout);
     for ( int i = 0; i < parts_.isize( ); i++ )
     {    const PartialSeqOnHyper& p = parts_[i];
          longlong e = p.Id2( );
	  PrintBlanks(cout, p.pos1());
	  edges[e].PrintBases(out, p.pos2(), p.Len(), False, UINT_MAX);
     }
     cout << "\n";    }

void SeqOnHyper::PrintStats( ostream& out, const vecbasevector& reads,
			      const vecbasevector& edges) const
{    basevector r = reads[Id1()];
     cout << "READ: "<< Id1() << "  ";
     if (Rc1()) {
       r.ReverseComplement();
       cout << "RC" << endl;
     } else
       cout << endl;
     for ( int i = 0; i < parts_.isize( ); i++ )
     {    const PartialSeqOnHyper& p = parts_[i];
          longlong e = p.Id2( );
	  out << "Edge: " << e << "  Len: " << edges[e].size();
	  out << "  Start: " << p.pos2() << "  End: " << p.Pos2();
	  out << "  Overlap: " << p.Len() << endl;
     }
}



// Stuff for CompressSeqOnHyper

longlong CompressedSeqOnHyper::ambig_count = 0;

// Not clever or fast.  Does this really not exist globally?
bool MatchUntilEnd( const basevector& s, const basevector& t,
		    int s_offset = 0, int t_offset = 0, 
		    bool rc1 = false ) {
  if( rc1 ) {
    for(uint i=0; i + s_offset < s.size() && i + t_offset < t.size(); i++)
      if( 3 - s[s.size()-1 - i - s_offset] != t[i + t_offset] ) return false;
    return true;
  }
  else {
    for(uint i=0; i + s_offset < s.size() && i + t_offset < t.size(); i++)
      if( s[i + s_offset] != t[i + t_offset] ) return false;
    return true;
  }
}

// Starting from position bases[pos], what edges out of vx agree perfectly?
// If rc, work as if passed the rc of the basevector.
void CompressedSeqOnHyper::AgreeingFromEdges( const HyperBasevector& hbv,
					      const int vx,
					      const int pos, const bool rc,
					      vec<int>& ans ) const {
  ans.clear();
  for(int j = 0; j < hbv.From(vx).isize(); j++)
    if( MatchUntilEnd( *mp_bases, hbv.EdgeObjectByIndexFrom(vx,j), pos, 0, rc ) ) {
      ans.push_back( hbv.EdgeObjectIndexByIndexFrom(vx,j) );
    }
}

void CompressedSeqOnHyper::Compress( const SeqOnHyper& salign,
				     const HyperBasevector& hbv_aligned,
				     const basevector* bv_aligned,
				     const vec<int>& to_right_vertex ) {

  // Set the things which can be set without looking at the HBV:
  mp_bases = bv_aligned;
  if( mp_ambig ) {
    delete mp_ambig;
    mp_ambig = NULL;
  }
  SetIdRc1( salign.Id1(), salign.Rc1() );

  int num_edges = salign.N();

  if( num_edges == 0 )
    return;

  m_first_edge_alignment = salign.Part(0);
  m_vx0right = to_right_vertex[ m_first_edge_alignment.Id2() ];

  if( num_edges == 1 )
    return;

  // Otherwise, we need to trace the basevector through the graph,
  // and record which choice to make if there's ambiguity.
  // Also, assert if the basevector doesn't actually align.

  int K = hbv_aligned.K();
  int pos = m_first_edge_alignment.Pos1() - K + 1; // I think
  int current_edge, current_vx = to_right_vertex[m_first_edge_alignment.Id2()];
  int vx_id = m_vx0right;
  int edge_id;
  vec<int> from_edges;
  vec<int> ambig;

  for(int i=1; i < num_edges; i++) {
    // These asserts are no longer valid - need to figure out what should replace them.
    //    ForceAssertEq( pos, salign.pos1(i) ); 
    //    ForceAssertLe( pos, bv_aligned->isize() - K ); // must be at least K bases
    edge_id = salign.Id2(i);
    // vx_id should now be the left vertex of edge_id

    AgreeingFromEdges( hbv_aligned, vx_id, pos, Rc1(), from_edges );

    if( ! Member( from_edges, edge_id ) ) { // uh oh
      cout << "SeqOnHyper doesn't align correctly to HyberBasevector!"
	   << "\nwhole basevector: " << mp_bases->ToString()
	   << "\npos=" << pos << ", " << (Rc1() ? "rc" : "fw")
	   << "\nedge " << edge_id << ": " 
	   << hbv_aligned.EdgeObject(edge_id).ToString()
	   << endl;
      FatalErr( "SeqOnHyper doesn't align correctly to HyperBasevector!" );
    }

    if( from_edges.size() != 1 )  // we have ambiguity, captain
      ambig.push_back(edge_id);
      
    // move to the right end of that edge.
    vx_id = to_right_vertex[edge_id];
    pos = salign.Pos1(i) - K + 1;
  }

  // If there was any ambiguity, allocate a vector to store it:
  if( ! ambig.empty() ) {
    mp_ambig = new vec<int>(ambig);
    IncrementAmbigCount();
  }

}

void CompressedSeqOnHyper::Compress(Bool rc1, longlong id1,
				    longlong id2, int pos1, int pos2, int len, 
				    const basevector* bv_aligned,
				    const vec<int>& to_right_vertex ) {
  mp_bases = bv_aligned;

  if( mp_ambig ) {
    delete mp_ambig;
    mp_ambig = NULL;
  }

  SetIdRc1(id1, rc1 );

  m_first_edge_alignment = PartialSeqOnHyper(id2, pos1, pos2, len );

  m_vx0right = to_right_vertex[ id2 ];
}

void CompressedSeqOnHyper::DecompressInto( SeqOnHyper& salign,
			   const HyperBasevector& hbv ) const {

  vec<PartialSeqOnHyper> parts( 1, m_first_edge_alignment );
  PartialSeqOnHyper partial;
  int vx_id = m_vx0right;
  int ambig_index = 0;
  int pos1;
  longlong id2;
  vec<int> from_edges;

  while( parts.back().Pos1() != mp_bases->isize() ) {

    pos1 = parts.back().Pos1() - hbv.K() + 1;

    AgreeingFromEdges( hbv, vx_id, pos1, Rc1(), from_edges );

    if( from_edges.size() == 1 )
      id2 = from_edges[0];
    else if ( from_edges.size() == 0 ) {
      FatalErr( "CompressedSeqOnHyper Decompression Failure!"
		<< "  No edge matches basevector!\n"
		<< "vertex=" << vx_id << ", pos=" << pos1
		<< " in bases=" << mp_bases->ToString() );
    }
    else { // ambiguity
      if( mp_ambig == NULL || ambig_index >= mp_ambig->isize() ) {
	FatalErr( "CompressedSeqOnHyper Decompression Failure!"
		  << "Found ambiguity and I don't know what to do!\n"
		  << "vertex=" << vx_id << ", pos=" << pos1
		  << " in bases=" << mp_bases->ToString() );
      }
      id2 = (*mp_ambig)[ambig_index++];
    }
    
    partial.Setpos1( pos1 );
    partial.Setpos2( 0 );
    partial.SetId2( id2 );
    partial.SetLen( min( hbv.EdgeObject(id2).isize(),
			 mp_bases->isize() - pos1 ) );
    parts.push_back(partial);

    // update vx_id
    vx_id = hbv.From(vx_id)[ Position(hbv.FromEdgeObj(vx_id), (int)id2) ];
  }

  salign.Set( Rc1(), Id1(), parts );
}


void BinaryWrite( int fd, const SeqOnHyper& s )
{    WriteBytes( fd, &s.rc1_, sizeof(s.rc1_) );
     WriteBytes( fd, &s.id1_, sizeof(s.id1_) );
     BinaryWrite( fd, s.parts_ );    }

void BinaryRead( int fd, SeqOnHyper& s )
{    ReadBytes( fd, &s.rc1_, sizeof(s.rc1_) );
     ReadBytes( fd, &s.id1_, sizeof(s.id1_) );
     BinaryRead( fd, s.parts_ );    }

void BinaryWrite( const String& filename, const vec<CompressedSeqOnHyper>& csaligns,
		  const HyperBasevector& hbv) {
  int fd = OpenForWrite(filename);
  int n = csaligns.size();
  WriteBytes( fd, &n, sizeof(int) );
  SeqOnHyper salign;
  for (int i = 0; i < n; i++) {
    csaligns[i].DecompressInto( salign, hbv);
    BinaryWrite(fd, salign);
  }
  Close(fd);
}

  
void BinaryRead( const String& filename, vec<CompressedSeqOnHyper>& csaligns,
		 const HyperKmerPath& h, const HyperBasevector& hbv,
		 const vecbasevector& reads ) {

  csaligns.clear();
  vec<int> to_right_vertex( h.EdgeObjectCount( ) );
  h.ToRight(to_right_vertex);
  
  int fd = OpenForRead(filename);
  int n;
  ReadBytes( fd, &n, sizeof(int) );
  csaligns.reserve(n);
  SeqOnHyper salign;
  for ( int i = 0; i < n; i++ ) {
    BinaryRead(fd, salign);
    csaligns.push_back(CompressedSeqOnHyper(salign, hbv, &reads[salign.Id1()],
					    to_right_vertex)); 
  }
  close(fd);
}
