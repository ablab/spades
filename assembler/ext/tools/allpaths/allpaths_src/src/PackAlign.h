///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// This file defines three classes for storing alignments between two DNA sequences:
// class align: stores a forward alignment
// class packalign: stores a forward alignment in a compressed form
// class flatalign: stores an alignment having no indels. [NOT DEFINED YET.]
// Added: class placement_mark.

// Class packalign is optimized to reduced memory usage to a bare minimum.
// A packalign stores the starting position of the alignment on both sequences 
// (pos1 and pos2), plus a sequence of integers 
//             (length, gap, length, ..., length, gap, length),
// where a "length" is the length of an aligning portion of the sequences
// (possibly with mismatches), and a "gap" is a positive integer if there's a gap on 
// the first sequence, else negative.  For backward compatibility we also allow
// a "gap" before the first length [denoted gap(0)], but it is almost always zero.

// Example:
//
//    Consider an alignment between the sequences ACGGTACGTTACTATTT and
//    AAAACGCTGTTTTAGTAT, given by the following picture:
//
//                       ACGGTACG TT ACTATTT
//                    AAAACGCT  GTTTTAGTAT
//
// Then the packalign object would have pos1 = 0, pos2 = 3, 
// (length 5, gap -2, length 1, gap 1, length 2, gap 1, length 5).

// (So length(0) = 5, gap(1) = -2, length(1) = 1, gap(2) = 1, etc.)

// The number of lengths is called nblocks.

// Externally, pos1, pos2, length(i), and gap(i) are handled as int's.
// Also, nblocks is an unsigned short.

// If we wish to declare a packalign "void", we set nblocks to 0 (see Kill( )).

// The memory optimization of packalign is achieved using different storage
// structures, depending on the sizes of pos1, pos2, and the other components
// of a packalign.  At present three are implemented:

// Type 0.  pos1 <= 1023 and pos2 <= 1023 and nblocks <= 6 and length(i) <= 1023
//          for all i and |gap(i)| <= 2 for all i.
//          storage: 12 bytes

// Type 1.  not type 0 and pos1 <= 4095 and pos2 <= 4095 and length(i) <= 4095 for 
//          all i and |gap(i)| <= 8 for all i.
//          storage: 12 + (2 * nblocks) bytes + overhead for one memory allocation

// Type 2.  not type 0 and not type 1.
//          storage: 16 + (8 * nblocks) bytes + overhead for one memory allocation
//
// Technical note.  Actually, any alignment with nblocks = 0 or gaps(0) != 0 is
// treated as type 2.  These alignments are only allowed for backward compatibility.

// As can be seen, packalign is memory-optimized for DNA sequences of length
// <= 1023 and alignments between them which have only a handful of gaps.

// Type 0 storage details:
//
// control -    3 bits (value = 0)
// pos1 -      10 bits
// pos2 -      10 bits
// nblocks -    3 bits
// gap(1) -     2 bits
// gap(2) -     2 bits
// gap(3) -     2 bits
//
// gap(4) -     2 bits
// length(0) - 10 bits
// length(1) - 10 bits
// length(2) - 10 bits
//
// gap(5) -     2 bits
// length(3) - 10 bits
// length(4) - 10 bits
// length(5) - 10 bits

// Type 1 storage details:
//
// control -    3 bits (value = 1)
// unused -     5 bits
// pos1 -      12 bits
// pos2 -      12 bits
// pointer -   64 bits
//
// at pointer: 
// nblocks -    16 bits
// gap(0) -     4  bits
// length(0) -  12 bits
// gap(1) -     4  bits
// length(1) -  12 bits
// ...

// Type 2 storage details:
//
// control -     3 bits (value = 2)
// unused -     13 bits
// nblocks -    16 bits
// pointer -    64 bits
//
// at pointer:
// pos1 -       32 bits
// pos2 -       32 bits
// length(0) -  32 bits
// gap(0) -     32 bits
// length(1) -  32 bits
// gap(1) -     32 bits
// ...

// To reduce running time, we also define a class align, which is the same as
// class packalign, but is not compressed.  Note that it keeps vectors of gaps
// and lengths, but that only the first nblocks of them are relevant.

#ifndef PACKALIGN
#define PACKALIGN

#include "Basevector.h"
#include "CommonSemanticTypes.h"
#include "CoreTools.h"
#include "Fastavector.h"
#include "Qualvector.h"
#include "ShortVector.h"
#include "feudal/MasterVec.h"
#include "feudal/SerfVec.h"
#include "math/HoInterval.h"
#include "pairwise_aligners/Mutmer.h"

const int Bits2  = 3, Bits3  = 7, Bits4 = 15, Bits10 = 1023, Bits12 = 4095, 
  Bits16 = 65535;

class short_pointer_or_words {
 public:
  union {
    unsigned short* p;
    unsigned int x[2];    
  };
};

class int_pointer_or_words {
 public:
  union {
    unsigned int* p;
    unsigned int x[2];    
  };
};

class align;

class packalign {

 public:

  // ========================================================================
  // Here are the packalign members which are MOST important for external use:
  // the constructor, a Set operator, and the reverse operation (Unpack).Note 
  // that gaps and and lengths are assumed to have the same size.
  // ========================================================================

  packalign( int pos1, int pos2, const avector<int>& gaps, 
             const avector<int>& lengths );

  packalign( const align& a );

  void Set( int pos1, int pos2, const avector<int>& gaps, 
            const avector<int>& lengths, int nblocks = -1 );

  void Unpack( int& pos1, int& pos2, avector<int>& gaps, 
               avector<int>& lengths ) const;

  void Unpack( int& pos1, int& pos2, avector<int>& gaps, 
               avector<int>& lengths, int& nblocks ) const;

  // ==========================================================================
  // The following members serve up interesting information about the
  // packalign, for convenience.  We let pos1 and Pos1 give the starting and
  // ending positions on the first sequence; pos2 and Pos2 give the starting
  // and ending positions on the second sequence.
  // ==========================================================================

  Bool Dead( ) const;

  int pos1( ) const;
  int Pos1( ) const;
  int pos2( ) const;
  int Pos2( ) const;

  int Offset( ) const { return pos1( ) - pos2( ); }

  int extent1( ) const { return Pos1( ) - pos1( ); }
  int extent2( ) const { return Pos2( ) - pos2( ); }
  ho_interval Extent1( ) const { return ho_interval( pos1( ), Pos1( ) ); }
  ho_interval Extent2( ) const { return ho_interval( pos2( ), Pos2( ) ); }

  int StartOnQuery() const { return pos1(); }
  int EndOnQuery() const { return Pos1(); }
  int StartOnTarget() const { return pos2(); }
  int EndOnTarget() const { return Pos2(); }

  // ========================================================================
  // The following members change a packalign.  Flip switches the role of the two
  // sequences.  Reverse produces the alignment that one would have if one 
  // reversed both sequences; it takes as input the lengths of the two sequences.
  // =======================================================================

  void Kill( );

  void Flip( );

  packalign Reverse( int b1_len, int b2_len ) const;
  void ReverseThis( int b1_len, int b2_len );

  void SetToFlipOf( const packalign& p );
  void SetToFlipOf( align a );
  void SetToReverseFlipOf( const packalign& p, int b1_len, int b2_len );
  void SetToReverseFlipOf( align a, int b1_len, int b2_len );

  // =========================================================================
  // The remaining public members are boring utilities:
  // =========================================================================

  packalign( ) { word_[0] = 7u << 29; } // type 7 means UNINITIALIZED

  packalign( const packalign& p );

  packalign& operator=( const packalign& p );

  ~packalign( );

  friend Bool operator==( const packalign& x, const packalign& y );

  // =========================================================================
  // Everything else is private!
  // =========================================================================

 private:

  void ConstructorCore( int pos1, int pos2,
                        const avector<int>& gaps, const avector<int>& lengths,
                        int nblocks = -1 );

  int Control( ) const { return (word_[0] >> 29) & Bits3; }

  void ConvertToType0( int pos1, int pos2, 
       const avector<int>& gaps, const avector<int>& lengths, int nblocks = -1 );

  void Unpack0( int& pos1, int& pos2, avector<int>& gaps, 
                avector<int>& lengths, int& n ) const;

  void DeleteType0( ) { }

  void ConvertToType1( int pos1, int pos2, 
       const avector<int>& gaps, const avector<int>& lengths, int nblocks = -1 );

  void Unpack1( int& pos1, int& pos2, avector<int>& gaps, 
                avector<int>& lengths, int& n ) const;

  void DeleteType1( ) 
  {
    short_pointer_or_words pw;
    pw.x[0] = word_[1];
    pw.x[1] = word_[2];
    delete [ ] pw.p;    
  }

  void ConvertToType2( int pos1, int pos2, 
                       const avector<int>& gaps, const avector<int>& lengths, 
                       int nblocks = -1 );

  void Unpack2( int& pos1, int& pos2, avector<int>& gaps, 
                avector<int>& lengths, int& n ) const;

  void DeleteType2( ) 
  {
    int_pointer_or_words pw;
    pw.x[0] = word_[1];
    pw.x[1] = word_[2];
    delete [ ] pw.p;    
  }

  unsigned int word_[3];

};

/**
    Class: align

    Alignment of query to target, possibly with gaps.
    Models type concept Align.
 */
class align {

 public:

  align( ): pos1_(0), pos2_(0), nblocks_(0)
  {
    gaps_.Setsize(0);
    lengths_.Setsize(0);
  }

  align( int pos1, int pos2, const avector<int>& gaps, 
         const avector<int>& lengths )
    : pos1_(pos1), pos2_(pos2), nblocks_( gaps.length ), gaps_(gaps), 
      lengths_(lengths) { }
  
    // copy ctor
  align( const align & a ) {
    pos1_ = a.pos1_;
    pos2_ = a.pos2_;
    nblocks_ = a.nblocks_;
    if ( (int) gaps_.length < nblocks_ ) {
      gaps_.Setsize(nblocks_);
      lengths_.Setsize(nblocks_);   
    }
    memcpy( gaps_.x, a.gaps_.x, nblocks_ * sizeof(int) );
    memcpy( lengths_.x, a.lengths_.x, nblocks_ * sizeof(int) );
  }

  void Set( int p1, int p2, const avector<int>& g, const avector<int>& l,
     int nblocks = -1 )
  {
    pos1_ = p1;
    pos2_ = p2;
    if ( nblocks < 0 ) nblocks_ = g.length;
    else nblocks_ = nblocks;
    if ( (int) gaps_.length < nblocks_ ) {
      gaps_.Setsize(nblocks_);
      lengths_.Setsize(nblocks_);    
    }
    memcpy( gaps_.x, g.x, nblocks_ * sizeof(int) );
    memcpy( lengths_.x, l.x, nblocks_ * sizeof(int) );    
  }

  align& operator=( const align& a )
  {
    pos1_ = a.pos1_;
    pos2_ = a.pos2_;
    nblocks_ = a.nblocks_;
    if ( (int) gaps_.length < nblocks_ ) {
      gaps_.Setsize(nblocks_);
      lengths_.Setsize(nblocks_);   
    }
    memcpy( gaps_.x, a.gaps_.x, nblocks_ * sizeof(int) );
    memcpy( lengths_.x, a.lengths_.x, nblocks_ * sizeof(int) );
    return *this;    
  }

  align( const packalign& p )
  {    p.Unpack( pos1_, pos2_, gaps_, lengths_, nblocks_ );    }

  int BinaryWrite(int fd) const {
    int written = SafeWrite(fd, & pos1_, sizeof(pos1_));
    written += SafeWrite(fd, &pos2_, sizeof(pos2_));
    written += SafeWrite(fd, &nblocks_, sizeof(nblocks_));
    written += SafeWrite(fd, gaps_.x, sizeof(gaps_.x[0])*nblocks_);
    written += SafeWrite(fd, lengths_.x, sizeof(lengths_.x[0]) *nblocks_);
    return written;
  }

  int BinaryRead(int fd) {
    int bytes_read = read(fd, &pos1_, sizeof(pos1_));
    bytes_read += read(fd, &pos2_, sizeof(pos2_));
    bytes_read += read(fd, &nblocks_, sizeof(nblocks_));
    gaps_.resize(nblocks_);
    lengths_.resize(nblocks_);
    bytes_read += read(fd, gaps_.x, sizeof(gaps_.x[0])*nblocks_);
    bytes_read += read(fd, lengths_.x, sizeof(lengths_.x[0]) *nblocks_);
    return bytes_read;
  }

  /// Set start position of the alignment on the query sequence; cryptic legacy interface
  void Setpos1( int p1 ) { pos1_ = p1; }
  /// Set start position of the alignment on the query sequence
  void SetStartOnQuery( int p1 ) { pos1_ = p1; }
  /// Set start position of the alignment on the target.
  void Setpos2( int p2 ) { pos2_ = p2; }
  /// Set start position of the alignment on the target.
  void SetStartOnTarget( int p2 ) { pos2_ = p2; }

  Bool FullLength( const int query_length )
  {    return pos1( ) == 0 && Pos1( ) == query_length;    }

  void UnpackFrom( const packalign& p );

  /// Start position of the alignment on the query sequence; cryptic legacy interface
  int pos1( ) const { return pos1_; }
  /// Start position of the alignment on the query sequence
  int StartOnQuery( ) const { return pos1_; }
  /// Start position of the alignment on the target sequence; cryptic legacy interface
  int pos2( ) const { return pos2_; }
  /// Start position of the alignment on the target sequence
  int StartOnTarget( ) const { return pos2_; }

  int Offset( ) const { return pos1_ - pos2_; }

  ///End position of the alignment on the query sequence
  int Pos1( ) const
  {
    int p1 = pos1( );
    for ( int j = 0; j < nblocks_; j++ ) {
      if ( gaps_(j) < 0 ) p1 -= gaps_(j);
      p1 += lengths_(j);    
    }
    return p1;    
  }

  ///End position of the alignment on the query sequence
  int EndOnQuery( ) const
  {
    int p1 = pos1( );
    for ( int j = 0; j < nblocks_; j++ ) {
      if ( gaps_(j) < 0 ) p1 -= gaps_(j);
      p1 += lengths_(j);    
    }
    return p1;    
  }

  ///End position of the alignment on the target sequence
  int Pos2( ) const
  {
    int p2 = pos2( );
    for ( int j = 0; j < nblocks_; j++ ) {   
      if ( gaps_(j) > 0 ) p2 += gaps_(j);
      p2 += lengths_(j);    
    }
    return p2;    
  }

  ///End position of the alignment on the target sequence
  int EndOnTarget( ) const
  {
    int p2 = pos2( );
    for ( int j = 0; j < nblocks_; j++ ) {   
      if ( gaps_(j) > 0 ) p2 += gaps_(j);
      p2 += lengths_(j);    
    }
    return p2;    
  }

  /// Advance start position of the alignment on the query sequence; legacy
  void AddToPos1( int a ) { pos1_ += a; }
  /// Advance start position of the alignment on the query sequence
  void AddToStartOnQuery( int a ) { pos1_ += a; }
  /// Advance start position of the alignment on the target sequence; legacy
  void AddToPos2( int a ) { pos2_ += a; }
  /// Advance start position of the alignment on the target sequence
  void AddToStartOnTarget( int a ) { pos2_ += a; }

  nbases_t extent1( ) const { return Pos1( ) - pos1( ); }
  nbases_t extent2( ) const { return Pos2( ) - pos2( ); }
  ho_interval Extent1( ) const { return ho_interval( pos1( ), Pos1( ) ); }
  ho_interval Extent2( ) const { return ho_interval( pos2( ), Pos2( ) ); }

  int Nblocks( ) const { return nblocks_; }

  void SetNblocks( int n ) 
  {    
    if ( n > (int) gaps_.length ) {
      gaps_.resize(n);
      lengths_.resize(n);    
    }
    nblocks_ = n;    
  }

  const avector<int>& Gaps( ) const { return gaps_; }
  const avector<int>& Lengths( ) const { return lengths_; }

  int Gaps( int i ) const { return gaps_(i); }
  int Lengths( int i ) const { return lengths_(i); }
     
  ///return position on 1 corresponding to pos2 on 2, -1 if not possible.
  ///not possible means pos2 is in indel or off end of alignment.
  int PosOn1(int pos2) const;
  ///return position on 2 corresponding to pos1 on 1, -1 if not possible.
  ///not possible means pos2 is in indel or off end of alignment.
  int PosOn2(int pos1) const;

  void SetGap( int i, int val ) 
  {
    AssertLt( i, nblocks_ );
    gaps_(i) = val;    
  }

  void SetLength( int i, int val ) 
  {
    AssertLt( i, nblocks_ );
    lengths_(i) = val;    
  }

  void AddToGap( int i, int addend ) { gaps_(i) += addend; }
  void AddToLength( int i, int addend ) { lengths_(i) += addend; }

  void ReverseThis( int b1_len, int b2_len );

  void Compactify( int len1, int len2 );

  void Kill( ) { SetNblocks(0); }

  void Flip( );

  /// Return new alignment corresponding to basevector 1 being trimmed.
  /// This is meant to be parallel to the call SetToSubOf(b1,startOn1,len);
  align TrimmedTo1(int startOn1, int len) const;

  /// Adjust the align object into sync with the 454-cycles of TACG.
  void Sync_to_TACG( const basevector & seq1,
                     const basevector & seq2,
                     Bool  isRC = False );

  // The following Read and Write functions are NOT correlated with each other.
  ///Not correlated with Write
  void Read( istream& in, int& errors, int& id1, int& id2, Bool& rc );

  ///Not correlated with Read
  void Write( ostream& out, int id1, int id2, Bool rc, int errors );

  friend Bool operator==( const align& a1, const align& a2 )
  {
    if ( a1.pos1_ != a2.pos1_ ) return False;
    if ( a1.pos2_ != a2.pos2_ ) return False;
    if ( a1.nblocks_ != a2.nblocks_ ) return False;
    for ( int i = 0; i < a1.nblocks_; i++ )
      if ( a1.gaps_(i) != a2.gaps_(i) 
           || a1.lengths_(i) != a2.lengths_(i) ) return False;
    return True; 
  }

  friend Bool operator<( const align& a1, const align& a2 )
  {    if ( a1.pos1_ < a2.pos1_ ) return True;
       if ( a1.pos1_ > a2.pos1_ ) return False;
       if ( a1.pos2_ < a2.pos2_ ) return True;
       if ( a1.pos2_ > a2.pos2_ ) return False;
       if ( a1.nblocks_ < a2.nblocks_ ) return True;
       if ( a1.nblocks_ > a2.nblocks_ ) return False;
       if ( a1.gaps_ < a2.gaps_ ) return True;
       if ( a1.gaps_ > a2.gaps_ ) return False;
       if ( a1.lengths_ < a2.lengths_ ) return True;
       return False;    }

  ///Return a pair with insertion and deletion information.
  /// first = gaps on 1, second = gaps on 2.
  std::pair<int,int> Gap1Gap2() const;

  Bool GapFree( ) const { return Gap1Gap2( ) == make_pair( 0, 0 ); } 

  /// Return the total number of errors.
  /// Calls MutationsGap1Gap2.
  int Errors( const basevector& rd1, const basevector& rd2 ) const;

  // Perfect is the same as Errors == 0, but faster.

  Bool Perfect( const basevector& rd1, const basevector& rd2 ) const;

  ///Return a three element vector with error information.
  /// elem 0 = substitutions, elem 1= gaps on 1, elem2 = gaps on 2.
  vector<int> MutationsGap1Gap2( const basevector& rd1, 
                                 const basevector& rd2 ) const;

  ///Return a three element vector with error information.
  /// elem 0 = substitutions, elem 1= gaps on 1, elem2 = gaps on 2.
  /// The errors will only be counted if they are between from1 and to1 on rd1
  /// and also between from2 and to2 on rd2.
  vector<int> MutationsGap1Gap2( const basevector& rd1, 
                                 int from1, int to1, const basevector& rd2, 
                                 int from2, int to2 ) const;

  int Mutations( const basevector& rd1, const basevector& rd2,
                 const qualvector& q1, int min_score ) const;
  void PrintMutations( const basevector& rd1, const basevector& rd2, ostream& log) const;
  int Indels( const basevector& rd1, const basevector& rd2,
              const qualvector& q1, int min_score ) const;
  int MatchingBases( const basevector& rd1, const basevector& rd2 );

  // PerfectIntervals1, PerfectIntervals2: return vector containing intervals of 
  // perfect match on first or second sequence

  void PerfectIntervals1( const basevector& rd1, const basevector& rd2,
                          vec<ho_interval>& perfs ) const;
  void PerfectIntervals2( const basevector& rd1, const basevector& rd2,
                          vec<ho_interval>& perfs ) const;
  void PerfectIntervals2( const fastavector& rd1, const fastavector& rd2,
                          vec<ho_interval>& perfs ) const;

  void CreateFromMutmers( int k, shortvector<mutmer>& m, const basevector& rd1, 
                          const basevector& rd2, int max_errors, float max_badness,
                          int local_max_errors, int end_stretch, int local_max_errors_done, 
                          int& errors_found );

  void CreateFromMutmersAndSW( int k, shortvector<mutmer>& m, 
                               const basevector& rd1, const basevector& rd2, int max_errors, 
                               int end_stretch, int& errors_found, bool affine_penalties );

  void CreateFromMutmersMT( int k, shortvector<mutmer>& m, const basevector& rd1, 
                            const basevector& rd2, int max_errors, float max_badness,
                            int local_max_errors, int end_stretch, int local_max_errors_done, 
                            int& errors_found );

 private:
  int pos1_, pos2_; // TODO: potentially dangerous truncation of indices
  int nblocks_;
  avector<int> gaps_, lengths_;

};  // class align

ostream & operator<<(ostream & os, const align & a);

Bool Proper( const align& a, int len1, int len2 );

void RequireProper( const align& a, int id1, int id2, Bool rc2, 
                    const vecbasevector& EE, int test_no, Bool fatal = True );

void RequireProper( const align& a, int id1, int id2, Bool rc2, 
                    const vec<int>& EE_length, int test_no, Bool fatal = True );

// Compute the number of errors in an alignment.

template<class BASEVEC1, class BASEVEC2>
int ActualErrors( const BASEVEC1& rd1, const BASEVEC2& rd2, 
                  const align& a, int mismatch_penalty = 1, int gap_penalty = 2 );

template<class BASEVEC1, class BASEVEC2>
int ActualErrors( Bool rc, const BASEVEC1& rd1, const BASEVEC2& rd2, 
                  const align& a, int mismatch_penalty = 1, int gap_penalty = 2 );

int ActualErrors( const basevector& rd1, const basevector& rd2, 
                  const align& a, int mismatch_penalty = 1, int gap_penalty = 2 );
int ActualErrors( const fastavector& rd1, const basevector& rd2, 
                  const align& a, int mismatch_penalty = 1, int gap_penalty = 2 );
int ActualErrors( const vec<char>& rd1, const basevector& rd2, 
                  const align& a, int mismatch_penalty = 1, int gap_penalty = 2 );
int ActualErrors( const fastavector& rd1, const fastavector& rd2, 
                  const align& a, int mismatch_penalty = 1, int gap_penalty = 2 );

int Bandwidth( align& a );

int CorrelatePositions( const align& a, int x1 );

/// Trim a basevector and an alignment on that basevector in sync.

void Trim1Together(const basevector & b1, const basevector & b2, 
		   const align & a, int startOn1, int len, 
		   basevector & trimmedb1, align & trimmeda); 

/// Class placement_mark: small-size structure to keep only the
/// target (contig) index/start position and query sequence orientation 
/// for the alignment of the query sequence against the target(s).

class placement_mark {

     public:

     placement_mark( ) { }

     /// Constructor. Create placement mark for an alignment of a query
     /// (no information on the query will be kept except the orientation!) 
     /// starting at position 
     /// \c pos on the target contig \c tig. If the alignment is for the 
     /// query itself (for reverse complemented query), use \c fw1=true
     /// (\c fw1=false).
     placement_mark( const int tig, const bool fw1, const unsigned int pos )
          : tig_or_( tig | ( fw1 ? TopBit32 : 0 ) ), pos_(pos) { }

     /// Set placement mark for a query
     /// (no information on the query will be kept except the orientation!) 
     /// aligning at position 
     /// \c pos on the target contig \c tig. If the alignment is for the 
     /// query itself (for reverse complemented query), use \c fw1=true
     /// (\c fw1=false).
     void Set( const int tig, const bool fw1, const unsigned int pos )
     {    tig_or_ = tig | ( fw1 ? TopBit32 : 0 );
          pos_ = pos;    }

     /// Index of the target sequence (usually contig); this is old interface
     int Tig( ) const { return tig_or_ & Bits31; }
     
     /// Index of the target sequence (usually contig) (more verbose
     /// alias for Tig(); 
     /// this method signature is expected to be preserved among
     /// different alignment classes: intended for use with generic interfaces)
     int TargetId( ) const { return tig_or_ & Bits31; }
     
     /// Returns \c true if query itself (fw) was aligned to the target; old interface.
     bool Fw1( ) const { return ( tig_or_ & TopBit32 ) != 0; }

     /// Returns \c true if query itself (fw) was aligned to the target
     /// (more verbose alias for Fw1(); this method signature is expected to 
     /// be preserved among different alignment classes: intended for use 
     /// with generic interfaces)
     Bool IsQueryFW( ) const { return ( tig_or_ & TopBit32 ) != 0; }

     /// Returns \c true if reverse complemented query was aligned to the target; 
     /// old interface.
     bool Rc1( ) const { return (tig_or_ & TopBit32) == 0 ; }

     /// Returns \c true if reverse complemented query was aligned to the target
     /// (more verbose alias for Rc11(); this method signature is expected to 
     /// be preserved among different alignment classes: intended for use 
     /// with generic interfaces)
     Bool IsQueryRC( ) const { return (tig_or_ & TopBit32) == 0 ; }

     /// Position on the target sequence; this is old interface.
     unsigned int Pos( ) const { return pos_; }
     
     /// Position on the target sequence (more verbose alias for Pos();
     /// this method signature is expected to be preserved among
     /// different alignment classes: intended for use with generic interfaces)
     unsigned int StartOnTarget() const { return pos_; }

     friend Bool operator<( const placement_mark& m1, const placement_mark& m2 )
     {    if ( m1.TargetId( ) < m2.TargetId( ) ) return True;
          if ( m1.TargetId( ) > m2.TargetId( ) ) return False;
          if ( m1.StartOnTarget( ) < m2.StartOnTarget( ) ) return True;
          if ( m1.StartOnTarget( ) > m2.StartOnTarget( ) ) return False;
          if ( m1.IsQueryFW( ) < m2.IsQueryFW( ) ) return True;
          return False;    }

     friend Bool operator==( const placement_mark& m1, const placement_mark& m2 )
     {    return m1.tig_or_ == m2.tig_or_ && m1.pos_ == m2.pos_;    }

     friend Bool operator!=( const placement_mark& m1, const placement_mark& m2 )
     {    return m1.tig_or_ != m2.tig_or_ || m1.pos_ != m2.pos_;    }

     friend ostream& operator<<( ostream& out, const placement_mark& m )
     {    return out << m.Tig( ) << " " << m.Pos( ) << " " 
               << ( m.Fw1( ) ? "fw" : "rc" ) << "\n";    }

     private:

     unsigned int tig_or_;
     unsigned int pos_;

};

TRIVIALLY_SERIALIZABLE(placement_mark);
typedef SerfVec<placement_mark> PlacementMarkVec;
typedef MasterVec<PlacementMarkVec> VecPlacementMarkVec;

/**
   Type concept: BasicAlign

   An alignment of query to target.   Defines methods that give
   of the query and the range of the target involved in the alignment.
   Does _not_ record whether the query or its rc aligns to the target.

   >int StartOnQuery() const;
   >int EndOnQuery() const;
   >int StartOnTarget() const;
   >int EndOnTarget() const;
 */

/**
   Type concept: Align

   An alignment of query to target.  Refines BasicAlign.
   Adds methods to tell whether the query itself or its
   rc aligns to the target, and the query and target ids.

   >Bool IsQueryFW( ) const;
   >Bool IsQueryRC( ) const;
   >int QueryId() const;
   >int TargetId() const;
   
 */

// Semantic Type: align_id_t
// Id of an Align (its index in a vector of alignments).
SemanticTypeStd( int, align_id_t );

// Semantic Type: nmuts_t
// Number of substitutions when aligning query to target
typedef int nmuts_t;


#endif
