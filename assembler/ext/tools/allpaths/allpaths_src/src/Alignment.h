///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// This file contains classes for storing alignments between two DNA sequences
// (basevectors), and various related functions.  These are derived from class
// packalign, in PackAlign.h and PackAlign.cc.  See also Alignment.cc.

// The example of an alignment in PackAlign.h is an example of a proper alignment:
// it goes to the end on one of the sequences on the left (the first sequence in the
// example), and goes to the end on one of the sequences on the right (the second
// sequence in the example).  For a proper alignment, we also require that the
// lengths and gaps are all nonzero.  See RequireProper, below.

// Class "alignment_plus" stores an alignment of two sequences, plus numerical
// identifiers for the sequences, plus a flag to indicate if the second sequence
// was reverse complemented, plus a score for the alignment.  This score is
// normally obtained from ScoreAlignment (see ScoreAlignment.cc).

// There is a variant class "augmented_alignment", and another variant class
// "nobbit", not in this file.  Some redundancy could be stripped out.

#ifndef ALIGNMENT
#define ALIGNMENT

#include "Basevector.h"
#include "CoreTools.h"
#include "Fastavector.h"
#include "PackAlign.h"
#include "Qualvector.h"
#include "ShortVector.h"
#include "math/Arith.h"
#include "pairwise_aligners/Mutmer.h"


// Sante --- Mon Oct 15 09:06:24 EDT 2001
//  I had to move BufSize in here, because it is used also in AlignsToLinks.
const int BufSize = 1000000/sizeof(int), MaxRecordSize = 10000/sizeof(int);



// ================================================================================
//
// Class alignment is a structure for storing an alignment of two sequences,
// which consists of the following data:
//
// 1.  A packalign (as described in PackAlign.h).
// 2.  An error count for the alignment.
//
// Note: Usage of the error count is spotty at best.  Some routines modify
// alignments without updating the error count.
//
// ================================================================================

class alignment : public packalign {

     public:

     alignment( ) : packalign( ) { }

     alignment( int pos1, int pos2, int errors,
          const avector<int>& gaps, const avector<int>& lengths )
          : packalign( pos1, pos2, gaps, lengths ) { errors_ = errors; }

     alignment( const packalign& p, int errors )
          : packalign(p), errors_(errors) { }

     void Set( int pos1, int pos2, int errors,
          const avector<int>& gaps, const avector<int>& lengths )
     {    packalign::Set( pos1, pos2, gaps, lengths );
          errors_ = errors;    }

     void Set( int pos1, int pos2, int errors,
          const avector<int>& gaps, const avector<int>& lengths, int nblocks )
     {    packalign::Set( pos1, pos2, gaps, lengths, nblocks );
          errors_ = errors;    }

     void Set( const align& al )
     {    packalign::Set( al.pos1( ), al.pos2( ), al.Gaps( ), al.Lengths( ),
               al.Nblocks( ) );    }

     void Set( const align& al, int errors )
     {    packalign::Set( al.pos1( ), al.pos2( ), al.Gaps( ), al.Lengths( ),
               al.Nblocks( ) );
          errors_ = errors;    }

     void Unpack( int& pos1, int& pos2, int& errors,
          avector<int>& gaps, avector<int>& lengths ) const
     {    packalign::Unpack( pos1, pos2, gaps, lengths );
          errors = errors_;    }

     void SetPosPosErrors( int pos1, int pos2, int errors )
     {    avector<int> gaps, lengths;
          int p1, p2;
          int err;
          Unpack( p1, p2, err, gaps, lengths );
          packalign::Set( pos1, pos2, gaps, lengths );
          errors_ = errors;    }

     int Errors( ) const { return errors_; }

     void SetErrors( int e ) { errors_ = e; }

     // Return a vector with three entries: the number of mutation errors, the
     // number of gaps on the first sequence, the number of gaps on the second
     // sequence.

     vector<int> MutationsGap1Gap2( const basevector& rd1, const basevector& rd2 );
     vector<int> MutationsGap1Gap2( const basevector& rd1,
          int from1, int to1, const basevector& rd2, int from2, int to2 );
     int Mutations( const basevector& rd1, const basevector& rd2,
          const qualvector& q1, int min_score );
     int Indels( const basevector& rd1, const basevector& rd2,
          const qualvector& q1, int min_score );

     void Compactify( int len1, int len2 ); // Remove zero gaps and lengths.

     void Write( ostream& out, int id1, int id2, Bool rc );
     void Read( istream& in, int& id1, int& id2, Bool& rc );

     void Write( ostream& out ) const;
     void Read( istream& in );

     private:

     int errors_;

};

template<class BASEVEC> int ActualErrors(const BASEVEC& rd1, 
     const BASEVEC& rd2, const alignment& a, int mismatch_penalty = 1, 
     int gap_penalty = 2)
{    return ActualErrors( rd1, rd2, align(packalign(a)), mismatch_penalty,
          gap_penalty );    }

const int Uninitialized = 1111111111;

class alignment_plus {

     public:

     alignment_plus( int read_id1,
		     int read_id2,
		     int rd1length,
		     int rd2length,
		     Bool if_read2_is_rc,
		     const alignment& a_arg,
		     float score_arg );

     alignment_plus( ) { read_id2_ = Uninitialized; }

     int Id1( ) const { return read_id1_; }
     int Id2( ) const
     {    if ( read_id2_ >= 0 ) return read_id2_;
          else return -read_id2_ - 1;    }
     void SetId1( int id1 ) { read_id1_ = id1; }
     void SetId2( int id2 ) { read_id2_ = id2; }
     Bool Rc2( ) const { return read_id2_ < 0; }

     int pos1( ) const { return a.pos1( ); }
     int pos2( ) const { return a.pos2( ); }
     int Pos1( ) const { return a.Pos1( ); }
     int Pos2( ) const { return a.Pos2( ); }
     int Extent1( ) const { return a.Pos1( ) - a.pos1( ); }
     int Extent2( ) const { return a.Pos2( ) - a.pos2( ); }

     // Note: SetRc2 must be called AFTER SetId2.

     void SetRc2( Bool rc )
     {    Assert( read_id2_ != Uninitialized );
          if (rc) read_id2_ = -read_id2_ - 1;    }

     // Use this method if you want to change Id2() without affecting Rc2().
     void SafeSetId2( int id2 ) {
       if ( read_id2_ < 0 )
         read_id2_ = -id2 - 1;
       else
         read_id2_ = id2;
     }

     float score;
     alignment a;

     // human-readable write, and a version which does not require rd2rc and
     // q2rc to exist in advance:

     void Print( Bool abbreviate, ostream& out, const basevector& rd1,
          const basevector& rd2, const basevector& rd2rc,
          const qualvector& q1, const qualvector& q2,
          const qualvector& q2rc, int begin = 0, Bool one_frame = False,
          Bool highlight_score = False );

     void Print( Bool abbreviate, ostream& out, const basevector& rd1,
          const basevector& rd2, const qualvector& q1,
          const qualvector& q2, int begin = 0, Bool one_frame = False,
          Bool highlight_score = False );

     // Swap converts a given alignment_plus between id1 and id2 into the
     // corresponding alignment_plus between id2 and id1.

     void Swap( int rd1length, int rd2length );

     void SetToSwapOf( const alignment_plus& x, int rd1length, int rd2length );

     void SetToSwapOf( const align& x, int id1, int id2, Bool rc2,
          float s, int rd1length, int rd2length );

     friend Bool operator<( const alignment_plus& p1, const alignment_plus& p2 )
     {    return p1.read_id1_ < p2.read_id1_ ||
          ( p1.read_id1_ == p2.read_id1_ && p1.a.pos1( ) < p2.a.pos1( ) );    }

     friend Bool operator==( const alignment_plus& p1, const alignment_plus& p2 )
     {    return ( p1.read_id1_ == p2.read_id1_ &&
                   p1.read_id2_ == p2.read_id2_ &&
                   p1.score == p2.score &&
                   p1.a == p2.a );    }

     void Write( ostream& out ) const;
     void Read( istream& in );

     private:

     // Note that read_id2_ contains both the second read id and the rc bit.

     int read_id1_, read_id2_;

};

#ifdef __DECCXX_VER
#pragma do_not_instantiate ostream& operator<<(ostream&, const vec<alignment_plus>&)
#pragma do_not_instantiate istream& operator>>(istream&, vec<alignment_plus>&)
#endif

// The following operator is designed to efficiently read in an vector
// of alignment_plus.

istream& operator>>( istream& in, vec<alignment_plus>& vecAligns );

// The following operator is designed to efficiently write out an
// entire (unindexed) vector of alignment_plus.

ostream& operator<<( ostream&, const vec<alignment_plus>& );

void WriteAppend( const String& f, const vec<alignment_plus>& v );

// The following function efficiently writes (or appends to) vectors
// of alignment_plus with some indexing.  The vector MUST be sorted,
// and, if appending, the id1 of the last alignment in the existing
// file can be no greater than the id1 of the first alignment in the
// vector to be written.

void WriteAlignsWithIndex( const String& filename, const vec<alignment_plus> &v );

void WriteAppendWithIndex( const String& filename, const vec<alignment_plus>& v );

// TODO: Potentially dangerous truncation of IDs
class augmented_alignment {

     public:

     alignment a;
     int RC;
     int length1, length2;
     int id1, id2;
     int pos1, pos2;
     float score;
     int Pos1, Pos2;

     augmented_alignment( ) { }
     augmented_alignment( const alignment& a_arg,
			  int RC_arg,
			  int length1_arg,
			  int length2_arg,
			  int id1_arg,
			  int id2_arg,
			  int pos1_arg,
			  int pos2_arg,
			  float score_arg,
			  int Pos1_arg,
			  int Pos2_arg )
       : a(a_arg),
	 RC(RC_arg),
	 length1(length1_arg),
	 length2(length2_arg),
	 id1(id1_arg),
	 id2(id2_arg),
	 pos1(pos1_arg),
	 pos2(pos2_arg),
	 score(score_arg),
	 Pos1(Pos1_arg),
	 Pos2(Pos2_arg)
      { }
};

const int OffTheEnd = -1;
const int AtGap = -2;

// CorrelatePositions: Given an alignment, and a position on the first read,
// determine the corresponding position on the second read.

int CorrelatePositions( const alignment& a, int x1 );

int ErrorsAt( const alignment& a, int x1, const basevector& rd1,
     const basevector& rd2 );

Float DepthOfCoverage( const vec<int>& reads, int contig_length,
     const vec<int>& read_lengths, const vec<alignment_plus>& all_aligns,
     const vec<int>& all_aligns_index );

Bool RequireProper( const alignment_plus& ap, const vecbasevector& EE,
     int test_no, Bool fatal = True );

Bool RequireProper( const alignment_plus& ap, const vec<int>& EE_length,
     int test_no, Bool fatal );

// Let ap1 be an alignment between reads id1 and id2.  Let ap2 be an alignment
// between reads id2 and id3.  Return the offset for an alignment between id1 and
// id3, obtained by "composing" the two given alignments.  The returned value is not
// meaningful unless the implicit overlap between id1 and id3 is positive.

int TransitiveOffset( const alignment_plus& ap1, const alignment_plus& ap2,
     int rd1length, int rd2length, int rd3length );

// Return an appropriate value for bandwidth, assuming that you have an alignment
// between two reads, but want to align them again using a banded Smith-Waterman.

int Bandwidth( alignment& a );

// Change a given alignment a, trimming it so that pos1 is incremented by n,
// or slightly more (in case the trimming would terminate in a gap).  The number
// of errors is not adjusted.

void TrimAlignmentFront( alignment& a, int n );
void TrimAlignmentFront( align& a, int n );

// Change a given align a, trimming it so that Pos1 is decremented by n,
// or slightly more (in case the trimming would terminate in a gap).

void TrimAlignmentBack( align& a, int n );

int MaxPerfectMatch( Bool rd1_is_rc, const align& a, const basevector& rd1,
     const basevector& rd2 );

#endif
