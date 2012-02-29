/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// A SeqOnHyper is a full-length perfect placement of a basevector on a
// HyperKmerPath.  It is represented as a sequence id, an orientation of it,
// and a sequence of triples:
// [position range on sequence (or rc), edge, position range on edge],
// where the edges are consecutive in the graph.

#ifndef SEQ_ON_HYPER_H
#define SEQ_ON_HYPER_H

#include "Basevector.h"
#include "CoreTools.h"
#include "CommonSemanticTypes.h"
#include "paths/HyperBasevector.h"
#include "paths/HyperKmerPath.h"

struct partial_perfect_align {
  longlong read_id;
  int edge_id;
  int read_length;
  int pos1;
  int pos2;
  int length;
};

class PartialSeqOnHyper {

     public:

     PartialSeqOnHyper( ) { }
     PartialSeqOnHyper( longlong id2, int pos1, int pos2, int len )
          : id2_(id2), pos1_(pos1), pos2_(pos2), len_(len)
     {    ForceAssertGe( id2, 0 );
          ForceAssertGe( pos1_, 0 );
          ForceAssertGe( pos2_, 0 );
          ForceAssertGt( len_, 0 );    }

     int pos1( ) const { return pos1_; }
     int Pos1( ) const { return pos1_ + len_; }
     int pos2( ) const { return pos2_; }
     int Pos2( ) const { return pos2_ + len_; }
     longlong Id2( ) const { return id2_; }
     int Len( ) const { return len_; }
     
     void Setpos1( int pos1 ) { pos1_ = pos1; }
     void Setpos2( int pos2 ) { pos2_ = pos2; }
     void SetId2( longlong id2 ) { id2_ = id2; }
     void SetLen( int len ) { len_ = len; }

     friend Bool operator==( 
          const PartialSeqOnHyper& p1, const PartialSeqOnHyper& p2 )
     {    return p1.id2_ == p2.id2_ && p1.pos1_ == p2.pos1_
               && p1.pos2_ == p2.pos2_ && p1.len_ == p2.len_;    }

     friend Bool operator<( 
          const PartialSeqOnHyper& p1, const PartialSeqOnHyper& p2 )
     {    if ( p1.id2_ < p2.id2_ ) return True;
          if ( p1.id2_ > p2.id2_ ) return False;
          if ( p1.pos1_ < p2.pos1_ ) return True;
          if ( p1.pos1_ > p2.pos1_ ) return False;
          if ( p1.pos2_ < p2.pos2_ ) return True;
          if ( p1.pos2_ > p2.pos2_ ) return False;
          if ( p1.len_ < p2.len_ ) return True;
          return False;    }

     private:

     longlong id2_;
     int pos1_;
     int pos2_;
     int len_;

};  // class PartialSeqOnHyper

class SeqOnHyper {

     public:

     SeqOnHyper( ) { }
     SeqOnHyper( Bool rc1, longlong id1, const vec<PartialSeqOnHyper>& parts )
          : rc1_(rc1), id1_(id1), parts_(parts) { }

     void Set( Bool rc1, longlong id1, const vec<PartialSeqOnHyper>& parts ) {
       rc1_ = rc1;
       id1_ = id1;
       parts_ = parts;
     }

     Bool Rc1( ) const { return rc1_; }
     int Len1( ) const { return parts_.back( ).Pos1( ); }
     int pos1( int i ) const { return parts_[i].pos1( ); }
     int Pos1( int i ) const { return parts_[i].Pos1( ); }
     int pos2( int i ) const { return parts_[i].pos2( ); }
     int Pos2( int i ) const { return parts_[i].Pos2( ); }
     longlong Id1( ) const { return id1_; }
     int Id2( int i ) const { return parts_[i].Id2( ); }
     int N( ) const { return parts_.size( ); }

     void SetRc1( Bool rc1 ) { rc1_ = rc1; }
     void SetId1( longlong id1 ) { id1_ = id1; }

     PartialSeqOnHyper& Part( int i ) { return parts_[i]; }
     const PartialSeqOnHyper& Part( int i ) const { return parts_[i]; }

  void TestValid( const HyperBasevector& hb, const vecbasevector & reads,
		  const vecbasevector & edges,
		  const vec<int>& to_left, const vec<int>& to_right ) const;

     void Print( ostream& out, const vec<int>& to_left, 
          const vec<int>& to_right ) const;

     void PrintVisual( ostream& out, const vecbasevector& reads,
		       const vecbasevector& edges) const;

     void PrintStats( ostream& out, const vecbasevector& reads,
		      const vecbasevector& edges) const;


     friend Bool operator==( const SeqOnHyper& s1, const SeqOnHyper& s2 )
     {    return s1.rc1_ == s2.rc1_ && s1.id1_ == s2.id1_ 
               && s1.parts_ == s2.parts_;    }

     friend Bool operator<( const SeqOnHyper& s1, const SeqOnHyper& s2 )
     {    if ( s1.rc1_ < s2.rc1_ ) return True;
          if ( s1.rc1_ > s2.rc1_ ) return False;
          if ( s1.id1_ < s2.id1_ ) return True;
          if ( s1.id1_ > s2.id1_ ) return False;
          if ( s1.parts_ < s2.parts_ ) return True;
          return False;    }

     friend ostream& operator<<( ostream& out, const SeqOnHyper& s ) {
       out << "read " << s.Id1() << (s.Rc1() ? "rc" : "fw") << " aligns  /";
       for(int i=0; i<s.N(); i++)
         out << "[" << s.pos1(i) << "-" << s.Pos1(i) 
	     << ") to edge " << s.Id2(i)
	     << " [" << s.pos2(i) << "-" << s.Pos2(i) << ")/";
       return out;
     }

     friend void BinaryWrite( int fd, const SeqOnHyper& s );
     friend void BinaryRead( int fd, SeqOnHyper& s );
     friend void BinaryWrite( int fd, const vec<SeqOnHyper>& s )
     {    BinaryWriteComplex( fd, s );    }
     friend void BinaryRead( int fd, vec<SeqOnHyper>& s )
     {    BinaryReadComplex( fd, s );    }

     private:

     Bool rc1_;
     longlong id1_;
     vec<PartialSeqOnHyper> parts_;
};  // class SeqOnHyper

/**
   Class: CompressedSeqOnHyper

   Alignment of a read to a HyperBasevector, represented in a compressed
   form.  The read aligns to a partial path through the HyperBasevector.
   We store the list of edges traversed, and the position on the first
   edge.

   Original comment:
   
   Given both the basevector and HyperBasevector that a SeqOnHyper
   refers to, we can compress the alignment object significantly.
   Usually it is enough to record the point on the HBv where the
   alignment begins.  If the HBv is badly behaved, we may need more
   data to avoid ambiguity; that data lives in a vec<int> which we
   allocate only if needed, so in the common (unambiguous) case, we
   only pay one pointer.
  
   As designed, the object remembers the pointer to the basevector,
   but the caller must pass in the HBv to decompress.
*/
class CompressedSeqOnHyper {

public:

  CompressedSeqOnHyper() : mp_ambig(NULL) {}
  CompressedSeqOnHyper( const SeqOnHyper& salign,
			const HyperBasevector& hbv_aligned,
			const basevector* bv_aligned,
			const vec<int>& to_right_vertex ) 
    : mp_ambig(NULL) 
  { 
    Compress( salign, hbv_aligned, bv_aligned, to_right_vertex );
  }

  // Constructor without SeqOnHyper - special case, on single edge
  CompressedSeqOnHyper( Bool rc1, longlong id1,
			int id2, int pos1, int pos2, int len, 
			const basevector* bv_aligned,
			const vec<int>& to_right_vertex) 
    : mp_bases(bv_aligned),
      mp_ambig(NULL),      
      m_first_edge_alignment(id2, pos1, pos2, len ),
      m_vx0right(to_right_vertex[id2])
  { 
    SetIdRc1(id1, rc1 );;
  }

  CompressedSeqOnHyper( const CompressedSeqOnHyper& other )
    : mp_bases              ( other.mp_bases ),
      mp_ambig              ( NULL ),
      m_first_edge_alignment( other.m_first_edge_alignment ),
      m_vx0right            ( other.m_vx0right ),
      m_idrc1               ( other.m_idrc1 )
  { if( other.mp_ambig != NULL ) mp_ambig = new vec<int>(*other.mp_ambig); }

  CompressedSeqOnHyper& operator=(const CompressedSeqOnHyper& rhs) {
    mp_bases = rhs.mp_bases;
    m_first_edge_alignment = rhs.m_first_edge_alignment;
    m_vx0right = rhs.m_vx0right;
    m_idrc1 = rhs.m_idrc1;
    if( rhs.mp_ambig == NULL ) {
      delete mp_ambig; 
      mp_ambig = NULL;
    }
    else if ( mp_ambig == NULL )
      mp_ambig = new vec<int>(*(rhs.mp_ambig));
    else
      (*mp_ambig) = *(rhs.mp_ambig);
    return *this;
  }

  ~CompressedSeqOnHyper() { delete mp_ambig; }

  friend bool operator==( const CompressedSeqOnHyper& lhs,
			  const CompressedSeqOnHyper& rhs ) {
    return ( lhs.mp_bases == rhs.mp_bases      // must point to same basevector!
	     && ( lhs.mp_ambig == rhs.mp_ambig // ie both NULL
		  || ( lhs.mp_ambig != NULL && rhs.mp_ambig != NULL // or both not NULL and...
		       && *(lhs.mp_ambig) == *(rhs.mp_ambig) ) )  // equal
	     && lhs.m_first_edge_alignment == rhs.m_first_edge_alignment
	     && lhs.m_vx0right == rhs.m_vx0right
	     && lhs.m_idrc1 == rhs.m_idrc1 );
  }

  // Compress from SeqOnHyper -- what is called by the constructor.
  void Compress( const SeqOnHyper& salign,
		 const HyperBasevector& hbv_aligned,
		 const basevector* bv_aligned,
		 const vec<int>& to_right_vertex );

  // Compress without SeqOnHyper - special case, on single edge
  void Compress(Bool rc1, longlong id1,
		longlong id2, int pos1, int pos2, int len, 
		const basevector* bv_aligned,
		const vec<int>& to_right_vertex );

  // Decompress into the offered SeqOnHyper.
  void DecompressInto( SeqOnHyper& salign,
		       const HyperBasevector& hbv_aligned ) const;

  // Data from the SeqOnHyper available without decompressing:
  longlong Id1() const { return ((m_idrc1 < 0) ? -m_idrc1-1 : m_idrc1); }
  Bool Rc1() const { return (m_idrc1 < 0); }
  const PartialSeqOnHyper& Part0() const { return m_first_edge_alignment; }

  // Does the basevector fall on a single edge?
  // Equivalent to decompressing and checking N()==1, but fast.
  bool SingleEdge() const 
  { return ( mp_bases->isize() == m_first_edge_alignment.Len() ); }

  // Just for fun, let's have a way to keep track of how much
  // ambiguity we've seen during compression:
private:
  static longlong ambig_count;
  static void IncrementAmbigCount() { ambig_count++; }
public:
  static void ClearAmbigCount() { ambig_count = 0; }
  static longlong AmbigCount() { return ambig_count; }


  // Custom versions of BinaryWrite and BinaryRead for vec<CompressedSeqOnHyper>
  // Actually stored on disk as vec<SeqOnHyper> and can be reloaded as such.
  // - Cannot store CompressedSeqOnHyper directly as it contains pointers to the
  // associated HyperBaseVector. You must supply BinaryRead with exactly the same
  // HyperKmerPath, HyperBaseVector and read set or it will complain loudly!
  friend void BinaryWrite( const String& filename, const vec<CompressedSeqOnHyper>& csaligns,
			   const HyperBasevector& hbv);
  friend void BinaryRead( const String& filename, vec<CompressedSeqOnHyper>& csaligns,
			  const HyperKmerPath& h, const HyperBasevector& hbv,
			  const vecbasevector& reads );

private:

  void SetIdRc1( longlong id1, bool rc1 ) { m_idrc1 = (rc1 ? -id1-1 : id1 ); }

  // Helper function: given a vertex in a HyperBasevector 
  // and a position on our basevector, what are all the edges
  // out from that vertex which agree with the basevector?
  // We're happiest if there is only one, but can't guarantee it.

  void AgreeingFromEdges( const HyperBasevector& hbv,
			  const int vx,
			  const int pos, const bool rc,
			  vec<int>& ans ) const;

private:

  const basevector* mp_bases;  // the basevector* passed in
  vec<int>* mp_ambig;          // NULL unless there is ambiguity

  PartialSeqOnHyper m_first_edge_alignment; // SeqOnHyper::Part(0)
  int m_vx0right;                           // vx at right end of edge 0

  longlong m_idrc1;   // holds SeqOnHyper fields id1_ and rc1_

};

#endif
