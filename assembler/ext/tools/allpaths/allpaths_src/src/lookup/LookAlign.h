///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef LOOK_ALIGN_H
#define LOOK_ALIGN_H

#include "math/Arith.h"
#include "Basevector.h"
#include "CoreTools.h"
#include "Fastavector.h"
#include "PackAlign.h"
#include "Qualvector.h"
#include "ReadPairing.h"
#include "SemanticTypes.h"
#include "SeqInterval.h"
#include "math/Functions.h"
#include "math/HoInterval.h"
#include "system/file/FileReader.h"

#include <set>

class look_align_plus;
class alignment_plus;

/** Alignment structure with information about query and target ids/lengths.

  \class look_align

  Alignment structure used in QueryLookupTable and associated code.

  The alignment object "a" is between the query sequence (as in pos1()
  and Pos1()) and the target sequence (as in pos2() and Pos2()).
  That is, the segment of the query sequence between positions pos1() and Pos1()
  is aligned to the segment of the target sequence between positions
  pos2() and Pos2(). Both [pos1,Pos1) and [pos2,Pos2) are half-open intervals.

  Models type concept Align.

  There are several methods to read and write look_aligns and the subclass
  look_align_plus:

  1. Class methods

  1.1. Binary methods:

  - look_align::BinaryWrite. Deprecated, not portable across architectures.
  - look_align::BinaryRead. Deprecated, not portable across architectures.
  - look_align::BinaryWritePortable. Files are identical for i386 and ia64.
  - look_align::BinaryReadPortable.  Files are identical for i386 and ia64.

  - look_align_plus::BinaryWrite. Files are identical for i386 and ia64.
  - look_align_plus::BinaryRead.  Files are identical for i386 and ia64.


  1.2. Text methods

  PrintXXX: many different output forms, see individual documentation.

  2. Vector methods

  2.1. Binary methods

  - LoadLookAlignBinary/WriteLookAlignBinary: read/write a binary vec<look_align>
  - LoadLookAlignPlusBinary/WriteLookAlignPlus Binary: read/write a binary
    vec<look_align_plus>, preserving the mutations_by_block information.


    2.2 Text methods

    - LoadLookAlignPlus: read a vec<look_align_plus> from a text file of
    alignments written out with PrintParseable(). The file is allowed to
    have additional information in it.


  See also: <Alignvector.h>

*/
class look_align {

 public:

    // TODO: potentially dangerous truncation of index
  int query_id, target_id;
  unsigned int query_length, target_length;
  /// This is only ever used in QueryLookupTable and can
  ///be ignored everywhere else.
  int nhits;
  nmuts_t mutations;
  int indels;
  align a;
  Bool rc1;

  look_align( )
    : query_id(0),
      target_id(0),
      query_length(0),
      target_length(0),
      nhits(0),
      mutations(0),
      indels(0),
      a(),
      rc1(false) { }

  look_align( int query_id_arg,
	      int target_id_arg,
	      unsigned int query_length_arg,
	      unsigned int target_length_arg,
	      Bool rc1_arg,
	      const align& a_arg,
	      int nhits_arg,
	      nmuts_t mutations_arg,
	      int indels_arg )
    : query_id(query_id_arg),
      target_id(target_id_arg),
      query_length(query_length_arg),
      target_length(target_length_arg),
      nhits(nhits_arg),
      mutations(mutations_arg),
      indels(indels_arg),
      a(a_arg),
      rc1(rc1_arg)
  { }

  /// Reset the align, mutations and indels parameters, keep
  /// the original target and query info.
  void ResetFromAlign(const align & a, const basevector &b1,
		      const basevector &b2);

  // Copy some member functions for align a.

  /// start position of the alignment on the query sequence; cryptic legacy interface
  int pos1( ) const { return a.pos1( ); }
  /// start position of the alignment on the query sequence; more verbose alias
  /// for pos1(). This method is expected to be present by generic templatized intefaces.
  int StartOnQuery( ) const { return a.StartOnQuery( ); }

  /// Sets alignment start position on the query sequence
  void SetStartOnQuery(int i) { a.SetStartOnQuery(i); }
  /// Shifts alignment start position on the query sequence
  void AddToStartOnQuery(int i) { a.AddToStartOnQuery(i); }
  /// start position of the alignment on the target sequence; cryptic legacy interface
  int pos2( ) const { return a.pos2( ); }
  ///  start position of the alignment on the target sequence; more verbose alias
  /// for pos2(). This method is expected to be present by generic templatized intefaces.
  int StartOnTarget() const { return a.StartOnTarget(); }

  /// sets alignment start position on the target
  void SetStartOnTarget(int i)  { a.SetStartOnTarget(i); }
  /// shifts alignment start position on the target
  void AddToStartOnTarget(int i)  { a.AddToStartOnTarget(i); }
  /// end position of the alignment on the query sequence; cryptic legacy interface
  int Pos1( ) const { return a.Pos1( ); }
  /// end position of the alignment on the query sequence; more verbose alias
  /// for Pos1(). This method is expected to be present by generic templatized intefaces.
  int EndOnQuery( ) const { return a.EndOnQuery( ); }
  /// end position of the alignment on the target sequence; cryptic legacy interface
  int Pos2( ) const { return a.Pos2( ); }
  /// end position of the alignment on the target sequence; more verbose alias
  /// for Pos2(). This method is expected to be present by generic templatized intefaces.
  int EndOnTarget( ) const { return a.EndOnTarget( ); }

  nbases_t extent1( ) const { return a.extent1( ); }
  nbases_t extent2( ) const { return a.extent2( ); }
  ho_interval Extent1( ) const { return a.Extent1( ); }
  ho_interval Extent2( ) const { return a.Extent2( ); }

  Bool Rc1( ) const { return rc1; }
  Bool Fw1( ) const { return !rc1; }

  Bool IsQueryRC( ) const { return Rc1(); }
  Bool IsQueryFW( ) const { return Fw1(); }

  int QueryId() const { return query_id; }
  int TargetId() const { return target_id; }
  void SetTargetId(int t_id) { target_id = t_id; }

  /// Returns length of the query, alignment of which is described by this object
  unsigned int QueryLength() const { return query_length; }

  /// Returns length of the target against which this alignment is made
  unsigned int TargetLength() const { return target_length; }

  /// Returns true if alignment is proper.
  Bool IsProper( ) const { return Proper( this->a, (int)this->query_length, (int)this->target_length ); }

  off_t BinaryWrite( int fd ) const
  {
    off_t bytes_written = 0;

    // First write all the "int" or "unsigned int" members of this class,
    // together with the three int members of class align.

    bytes_written += SafeWrite( fd, &query_id,
				sizeof(query_id) + sizeof(target_id)
				+ sizeof(query_length) + sizeof(target_length) + sizeof(nhits)
				+ sizeof(mutations) + sizeof(indels) + 3 * sizeof(int) );

    // Write the number of blocks in the alignment, followed by the gaps
    // and lengths elements of the alignments.

    int nblocks = a.Nblocks( );
    bytes_written += SafeWrite( fd, &nblocks, sizeof(nblocks) );
    bytes_written += SafeWrite( fd, &(a.Gaps( )(0)), nblocks * sizeof(int) );
    bytes_written += SafeWrite( fd, &(a.Lengths( )(0)),
				nblocks * sizeof(int) );

    // Finally, write rc1.

    bytes_written += SafeWrite( fd, &rc1, sizeof(rc1) );
    return bytes_written;    }


  ///This method produces the same files under i386 and ia64.
  ///BinaryRead and BinaryWrite operate differently on i386 and ia64.
  off_t BinaryWritePortable( int fd ) const
  {
    off_t bytes_written = 0;

    // First write all the "int" or "unsigned int" members of this class,

    bytes_written += SafeWrite( fd, &query_id, sizeof(query_id));
    bytes_written += SafeWrite( fd, &target_id, sizeof(target_id));
    bytes_written += SafeWrite( fd, &query_length, sizeof(query_length));
    bytes_written += SafeWrite( fd, &target_length, sizeof(target_length));
    bytes_written += SafeWrite( fd, &nhits, sizeof(nhits));
    bytes_written += SafeWrite( fd, &mutations, sizeof(mutations));
    bytes_written += SafeWrite( fd, &indels, sizeof(indels));
    bytes_written += a.BinaryWrite(fd);

    // Finally, write rc1.
    bytes_written += SafeWrite( fd, &rc1, sizeof(rc1) );
    return bytes_written;
  }

  ///This method produces the same files under i386 and ia64.
  ///BinaryRead and BinaryWrite operate differently on i386 and ia64.
  off_t BinaryReadPortable( int fd )
  {
    off_t bytes_read = 0;

    // First read all the "int" or "unsigned int" members of this class,

    bytes_read += read( fd, &query_id, sizeof(query_id));
    bytes_read += read( fd, &target_id, sizeof(target_id));
    bytes_read += read( fd, &query_length, sizeof(query_length));
    bytes_read += read( fd, &target_length, sizeof(target_length));
    bytes_read += read( fd, &nhits, sizeof(nhits));
    bytes_read += read( fd, &mutations, sizeof(mutations));
    bytes_read += read( fd, &indels, sizeof(indels));
    bytes_read += a.BinaryRead(fd);
    // Finally, read rc1.
    bytes_read += read( fd, &rc1, sizeof(rc1) );
    return bytes_read;
  }

  void BinaryRead( int fd )
  {
    read( fd, &query_id, sizeof(query_id) + sizeof(target_id)
          + sizeof(query_length) + sizeof(target_length) + sizeof(nhits)
          + sizeof(mutations) + sizeof(indels) + 3 * sizeof(int) );
    int nblocks;
    read( fd, &nblocks, sizeof(nblocks) );
    a.SetNblocks(nblocks);
    vec<int> gaps, lengths;
    gaps.resize(nblocks), lengths.resize(nblocks);
    read( fd, &gaps[0], nblocks * sizeof(int) );
    read( fd, &lengths[0], nblocks * sizeof(int) );
    for ( int i = 0; i < nblocks; i++ )
    {
      a.SetGap( i, gaps[i] );
      a.SetLength( i, lengths[i] );
    }
    read( fd, &rc1, sizeof(rc1) );
  }

  Bool FullLength( ) const
  {    return a.StartOnQuery( ) == 0 && a.EndOnQuery( ) == (int) query_length;    }

  Bool BeginOfTarget( ) const
  {    return a.StartOnTarget( ) == 0;    }

  Bool EndOfTarget( ) const
  {    return a.EndOnTarget( ) == (int) target_length;    }

  int Errors( ) const { return mutations + indels; }

  float ErrorRate( ) const
  {    return float( mutations + indels ) / float( a.Pos1( ) - a.pos1( ) );    }

  float MutationRate( ) const
  {    return float(mutations) / float( a.EndOnQuery( ) - a.StartOnQuery( ) );    }

  int GaplessLength( ) const
  {    return a.EndOnQuery() - a.StartOnQuery() + a.EndOnTarget() - a.StartOnTarget() - indels; }

  /// Return a look_align corresponding to a trimmed basevector b1.
  /// This is meant to be parallel to the call SetToSubOf(b1,startOn1,len);
  look_align TrimmedTo1(int startOn1, int len,
			const basevector & b1, const basevector & b2) const;

  /// Return a look_align_plus corresponding to a trimmed basevector b1.
  /// This is meant to be parallel to the call SetToSubOf(b1,startOn1,len);
  look_align_plus TrimmedTo1Plus(int startOn1, int len,
				 const basevector & b1,
				 const basevector & b2) const;

  ///Transform this look_align into a seq_interval.
  ///The interval is given on the reference (sequence 2).
  ///The id of the seq_interval is the query_id, just in case it's useful.
  operator seq_interval() const {
    return seq_interval(query_id, target_id, pos2(), Pos2());
  }

  friend Bool operator==( const look_align& q1, const look_align& q2 ) {
    if ( q1.query_id == q2.query_id &&
	 q1.target_id == q2.target_id &&
	 q1.query_length == q2.query_length &&
	 q1.target_length == q2.target_length &&
	 q1.nhits == q2.nhits &&
	 q1.mutations == q2.mutations &&
	 q1.indels == q2.indels &&
	 q1.a == q2.a &&
	 q1.rc1 == q2.rc1 ) return True;
    return False;
  }

  friend Bool operator!=( const look_align& q1, const look_align& q2 )
  {    return !( q1 == q2 );    }

  friend Bool operator<( const look_align& q1, const look_align& q2 ) {
    if ( q1.query_id < q2.query_id ) return True;
    if ( q1.query_id > q2.query_id ) return False;
    if ( q1.a.pos1( ) < q2.a.pos1( ) ) return True;
    if ( q1.a.pos1( ) > q2.a.pos1( ) ) return False;
    if ( q1.a.Pos1( ) > q2.a.Pos1( ) ) return True;
    if ( q1.a.Pos1( ) < q2.a.Pos1( ) ) return False;
    if ( q1.mutations < q2.mutations ) return True;
    if ( q1.mutations > q2.mutations ) return False;
    if ( q1.indels < q2.indels ) return True;
    if ( q1.indels > q2.indels ) return False;
    return False;
  }

  void PrintParseable( ostream& out, const basevector& query,
		       const qualvector& query_qual, const basevector& target,
		       int start_on_target, const vec<String>& seq_names,
		       const vec<String>& target_names ) const;

  /// Prints alignment in parseable format (qltout); if <query> and <target> are not specified,
  /// total number of mutations will be saved in the first block with all the remaining blocks
  /// (if any) recorded as having no mutations. If <query> and <target> are specified,
  /// then the correct numbers of mutations in every continuous block of gapped alignment (with indels)
  /// will be written (there is no difference for gapless aligns since they have only one block!).
  void PrintParseable( ostream &out,
		       const basevector *query=0,
		       const basevector *target=0,
		       bool endl=true) const;

  /// Delegates to ReadParseable(istream &)
  void ReadParseable( const String& in );
  Bool ReadParseableOrFail( const String& in );

  /// Read it in from a one-line description.
  void ReadParseable( istream & ins );

  void PrintParseableBrief( ostream& out, const basevector& query,
			    const qualvector& query_qual, const basevector& target,
			    int start_on_target, const vec<String>& seq_names,
			    const vec<String>& target_names ) const;

  ///Print a human-readable summary on one line.
  ///Note that only out and the XXX_names parameters are actually used!
  void PrintReadableBrief( ostream& out, const basevector& query,
			   const qualvector& query_qual, const basevector& target,
			   int start_on_target, const vec<String>& seq_names,
			   const vec<String>& target_names ) const;

  ///Print a human-readable readable summary on one line.
  void PrintReadableBrief( ostream& out, const String & query_name,
			   const String & target_name ) const;

  ///Print a human-readable readable summary with ids instead of names.
  void PrintReadableBrief( ostream& out ) const;

  Float RmrPercent( const basevector& query,
		    const qualvector& query_qual, const basevector& target,
		    int start_on_target ) const;

  Float MMPercent( const basevector& query,
		   const qualvector& query_qual, const basevector& target,
		   int start_on_target ) const;

  void PrintRmr( ostream& out, const basevector& query,
		 const qualvector& query_qual, const basevector& target,
		 int start_on_target ) const;

  void PrintMM( ostream& out, const basevector& query,
		const qualvector& query_qual, const basevector& target,
		int start_on_target ) const;

  void PrintRmrByBlock( ostream& out, const basevector& query,
			const qualvector& query_qual, const basevector& target,
			int start_on_target, const vec<String>& seq_names,
			const vec<String>& target_names ) const;

  int CountNqs( const basevector& query, const qualvector& query_qual,
		const basevector& target, int start_on_target ) const;

  int CountNqsWithReturns( const basevector& query, const qualvector& query_qual,
			   const basevector& target, int start_on_target,
			   const vec<String>& seq_names,
			   const vec<String>& target_names, int &qualitySnpBase,
			   int  &targetSnpBase );

  void PrintNqs( ostream& out, const basevector& query,
		 const qualvector& query_qual, const basevector& target,
		 int start_on_target, const vec<String>& seq_names,
		 const vec<String>& target_names ) const;

  Float QualScore( const basevector& query,
		   const qualvector& query_qual, const basevector& target,
		   int start_on_target, const qualvector& qual_g ) const;

  void PrintQualScore( ostream& out, const basevector& query,
		       const qualvector& query_qual, const basevector& target,
		       int start_on_target, const qualvector& qual_g ) const;

  // PrintVisual: the option reverse_display causes the alignment to be shown with
  // the read forward even if the alignment is rc.  But you have to supply the
  // rc of the target too.

  void PrintVisual( ostream& out, const basevector& query,
		    const qualvector& query_qual, const basevector& target,
		    int start_on_target, Bool abbr = True ) const;
  void PrintVisual( ostream& out, const basevector& query,
		    const qualvector& query_qual, const fastavector& target,
		    int start_on_target, Bool abbr = True ) const;
  void PrintVisual( ostream& out, const fastavector& query,
		    const qualvector& query_qual, const basevector& target,
		    int start_on_target, Bool abbr = True ) const;
  void PrintVisual( ostream& out, const fastavector& query,
		    const qualvector& query_qual, const fastavector& target,
		    int start_on_target, Bool abbr = True ) const;

  void PrintVisual( ostream& out, const basevector& query,
		    const qualvector& query_qual, const qualvector& target_qual,
                    const basevector& target,
		    int start_on_target, Bool abbr = True ) const;
  void PrintVisual( ostream& out, const basevector& query,
		    const qualvector& query_qual, const qualvector& target_qual,
                    const fastavector& target,
		    int start_on_target, Bool abbr = True ) const;
  void PrintVisual( ostream& out, const fastavector& query,
		    const qualvector& query_qual, const qualvector& target_qual,
                    const basevector& target,
		    int start_on_target, Bool abbr = True ) const;
  void PrintVisual( ostream& out, const fastavector& query,
		    const qualvector& query_qual, const qualvector& target_qual,
                    const fastavector& target,
		    int start_on_target, Bool abbr = True ) const;

  void PrintVisual( ostream& out, const basevector& query,
		    const basevector& target, Bool abbr = True ) const;
  void PrintVisual( ostream& out, const fastavector& query,
		    const basevector& target, Bool abbr = True ) const;

  void PrintVisual( ostream& out, const basevector& query,
		    const qualvector& query_qual, const basevector& target,
                    const basevector& target_rc, int start_on_target,
                    const Bool abbr = True,
                    const Bool reverse_display = False ) const;

  void PrintVisual( ostream& out, const basevector& query,
		    const basevector& target, const basevector& target_rc,
                    const Bool abbr = True,
                    const Bool reverse_display = False ) const;

};  // class look_align


inline void BinaryWrite( int fd, const look_align& a ) { a.BinaryWritePortable( fd ); }
inline void BinaryRead( int fd, look_align& a ) { a.BinaryReadPortable( fd ); }
inline void BinaryWrite( int fd, const vec<look_align>& v ) { BinaryWriteComplex( fd, v ); }
inline void BinaryRead( int fd, vec<look_align>& v ) { BinaryReadComplex( fd, v ); }

/// Load look_aligns from a binary file.
///  q_ids: if not null, load only hits with given query_id (must be sorted)
///  t_ids: if not null, load only hits with given target_id (must be sorted)
void LoadLookAlignBinary( const String &file_name,
			  vec<look_align> &hits,
			  const set<int> *q_ids = 0,
			  const set<int> *t_ids = 0 );

/// Write look_aligns to a binary file.
///  q_ids: if not null, load only hits with given query_id
///  t_ids: if not null, load only hits with given target_id
void WriteLookAlignBinary( const String &file_name,
			   const vec<look_align> &hits,
			   const set<int> *q_ids = 0,
			   const set<int> *t_ids = 0 );

/// Return true if the pair of look_aligns best, secondBest represent
/// ambiguous alignments of a read.  Two tests are supported: by
/// difference or ratio of #errors.
inline bool IsAmbiguous(const look_align &best, const look_align &secondBest,
			int ERR_DIFF, double ERR_RATE_MULT = 4.0, bool COMPARE_BY_RATE = false)
{
  if (COMPARE_BY_RATE)
    return (secondBest.ErrorRate()
	    <= ERR_RATE_MULT * max(best.ErrorRate(), float(1.0/best.query_length)));
  return (secondBest.Errors() <= best.Errors() + ERR_DIFF);
}

/// Builds index into a vec of look aligns.
void BuildLookAlignsIndex( const vec<look_align>& aligns,
			   vec< vec<align_id_t> >& aligns_index, int nqueries );


/// Class look_align_plus is designed to mirror the output of PrintParseable.

class look_align_plus : public look_align {

 public:

  vec<int> mutations_by_block;

  look_align_plus( )
    : look_align(), mutations_by_block() { }

  look_align_plus(const look_align & orig)
    : look_align(orig), mutations_by_block(orig.a.Nblocks()) {}

  void ReadParseable( const String& in );

  void WriteParseable( ostream& out ) const;

  ///Write a complete look_align_plus portably, matches BinaryRead().
  off_t BinaryWrite(int fd) const;

  ///Read a complete look_align_plus portably, matches BinaryWrite().
  off_t BinaryRead(int fd);
};

/// Load look_aligns from a file, and generate an index to them.

void LoadLookAligns( const String& file_name, vec<look_align>& aligns,
                     vec< vec<align_id_t> >& aligns_index, unsigned long nqueries );

/// Same, but skip the index.

void LoadLookAligns( const String& file_name, vec<look_align>& aligns,
     const Bool ignore_bads = False );

/// Write out the look_aligns to a file in parseable text form.
void WriteLookAligns( const String& file_name, const vec<look_align>& aligns );


/// Load look_align_pluses from file.
///  q_ids: if not null, load only hits with given query_id (must be sorted)
///  t_ids: if not null, load only hits with given target_id (must be sorted)
void LoadLookAlignPlus( const String &file_name,
			vec<look_align_plus> &hits,
			const vec<int> *q_ids = 0,
			const vec<int> *t_ids = 0 );

/// Load look_align_pluses from a binary file.
///  q_ids: if not null, load only hits with given query_id (must be sorted)
///  t_ids: if not null, load only hits with given target_id (must be sorted)
void LoadLookAlignPlusBinary( const String &file_name,
                              vec<look_align_plus> &hits,
                              const set<int> *q_ids = NULL,
                              const set<int> *t_ids = NULL );

/// Write look_align_pluses to a binary file.
///  q_ids: if not null, load only hits with given query_id
///  t_ids: if not null, load only hits with given target_id
void WriteLookAlignPlusBinary( const String &file_name,
                               const vec<look_align_plus> &hits,
                               const set<int> *q_ids = NULL,
                               const set<int> *t_ids = NULL );


/// Function: LoadLookAlignPlusBinary
///
/// Load look_align_pluses from an open binary file.
///  q_ids: if not null, load only hits with given query_id (must be sorted)
///  t_ids: if not null, load only hits with given target_id (must be sorted)
void LoadLookAlignPlusBinary( int fd,
                              vec<look_align_plus> &hits,
                              const set<int> *q_ids = NULL,
                              const set<int> *t_ids = NULL );

/// Function: WriteLookAlignPlusBinary
///
/// Write look_align_plusses to an open binary file.
///
/// Parameters:
///  fd - a file descriptor, open for writing
///  q_ids - if not null, load only hits with given query_id
///  t_ids - if not null, load only hits with given target_id
void WriteLookAlignPlusBinary( int fd,
			       const vec<look_align_plus> &hits,
			       const set<int> *q_ids = NULL,
			       const set<int> *t_ids = NULL );

/**
   Class: GaplessAlign

   Represents a simple gapless alignment of a query to a target.

   A vector of such alignments and be read and written to disk faster
   then a vec< look_align >.
 */
class GaplessAlign {
 public:
  int query_id, target_id;
  int startOnQuery, endOnQuery,
    startOnTarget, endOnTarget;
  Bool rc1;

 GaplessAlign():
  query_id( -1 ), target_id( -1 ),
    startOnQuery( -1 ), endOnQuery( -1 ),
    startOnTarget( -1 ), endOnTarget( -1 ),
    rc1( False ) { }

  GaplessAlign( int _query_id, int _target_id,
		int _startOnQuery, int _endOnQuery,
		int _startOnTarget, int _endOnTarget,
		Bool _rc1 ):
  query_id( _query_id ), target_id( _target_id ),
    startOnQuery( _startOnQuery ), endOnQuery( _endOnQuery ),
    startOnTarget( _startOnTarget ), endOnTarget( _endOnTarget ),
    rc1( _rc1 ) { }

 GaplessAlign( const look_align& la ):
  query_id( la.query_id ), target_id( la.target_id ),
    startOnQuery( la.StartOnQuery() ), endOnQuery( la.EndOnQuery() ),
    startOnTarget( la.StartOnTarget() ), endOnTarget( la.EndOnTarget() ),
    rc1( la.Rc1() ) {
   ForceAssert( la.indels == 0 );
 }

  int QueryId() const { return query_id; }
  int TargetId() const { return target_id; }
  int StartOnQuery() const { return startOnQuery; }
  int EndOnQuery() const { return endOnQuery; }
  int StartOnTarget() const { return startOnTarget; }
  int EndOnTarget() const { return endOnTarget; }
  Bool Rc1() const { return rc1; }
  Bool Fw1() const { return !rc1; }

  Bool IsQueryRC( ) const { return Rc1(); }
  Bool IsQueryFW( ) const { return Fw1(); }

};  // class GaplessAlign


class genome_pos {

 public:

  int pos1, Pos1;
  int t; // target contig
  int pos2, Pos2;
  int rc1;
  float mismatch_percent;

  genome_pos( ) : pos1(-1), Pos1(-1), t(-1), pos2(-1), Pos2(-1), rc1(-1),
		  mismatch_percent(-1) { }
  genome_pos( int pos1_arg,
	      int Pos1_arg,
	      int t_arg,
	      int pos2_arg,
	      int Pos2_arg,
	      int rc1_arg,
	      float mismatch_percent_arg )
    : pos1(pos1_arg),
      Pos1(Pos1_arg),
      t(t_arg),
      pos2(pos2_arg),
      Pos2(Pos2_arg),
      rc1(rc1_arg),
      mismatch_percent(mismatch_percent_arg)
  { }

  friend istream& operator>>( istream& in, genome_pos& x )
  {
    return in >> x.pos1 >> x.Pos1 >> x.t >> x.pos2 >> x.Pos2 >> x.rc1
              >> x.mismatch_percent;
  }

  friend ostream& operator<<( ostream& out, const genome_pos& x )
  {
    return out << x.pos1 << " " << x.Pos1 << " " << x.t << " "
               << x.pos2 << " " << x.Pos2 << " " << x.rc1 << " "
               << x.mismatch_percent
               << "\n";
  }

  friend void LoadUniqueLocs( const String& run_dir, vec<genome_pos>& up,
			      const String& basename = "reads" )
  {
    int N = MastervecFileObjectCount( run_dir + "/reads.fastb" );
    up.resize(N);
    FileReader fr( (run_dir + "/" + basename + ".uniquely_placed").c_str() );
    fr.read( &up[0], up.size( ) * sizeof(genome_pos) );
  }

};

class genome_pos_plus {

 public:

  int pos1, Pos1;
  int t; // target contig
  int pos2, Pos2;
  int length2; // target length
  int rc1;
  float mismatch_percent;

  genome_pos_plus( ) : pos1(-1), Pos1(-1), t(-1), pos2(-1), Pos2(-1),
		       length2(-1), rc1(-1), mismatch_percent(-1) { }
  genome_pos_plus( int pos1_arg,
		   int Pos1_arg,
		   int t_arg,
		   int pos2_arg,
		   int Pos2_arg,
		   int length2_arg,
		   int rc1_arg,
		   float mismatch_percent_arg )
    : pos1(pos1_arg),
      Pos1(Pos1_arg),
      t(t_arg),
      pos2(pos2_arg),
      Pos2(Pos2_arg),
      length2(length2_arg),
      rc1(rc1_arg),
      mismatch_percent(mismatch_percent_arg)
  { }

  friend istream& operator>>( istream& in, genome_pos_plus& x )
  {
    return in >> x.pos1 >> x.Pos1 >> x.t >> x.pos2 >> x.Pos2
              >> x.length2 >> x.rc1 >> x.mismatch_percent;
  }

  friend ostream& operator<<( ostream& out, const genome_pos_plus& x )
  {
    return out << x.pos1 << " " << x.Pos1 << " " << x.t << " "
               << x.pos2 << " " << x.Pos2 << " " << x.length2 << " "
               << x.rc1 << " " << x.mismatch_percent << "\n";
  }

  friend void LoadUniqueLocs( const String& run_dir,
			      const String& placement_file, vec<genome_pos_plus>& up )
  {
    int N = MastervecFileObjectCount( run_dir + "/reads.fastb" );
    up.resize(N);
    FileReader fr( (run_dir + "/" + placement_file).c_str() );
    fr.read( &up[0], up.size( ) * sizeof(genome_pos_plus) );
  }

};

/// Sort by: begin_on_target (begin defined by start of aligned portion).
struct order_lookalign_Begin
  : public binary_function<const look_align&, const look_align&, bool>
{
 public:
  bool operator() ( const look_align &left, const look_align &right ) const {
    return ( left.a.pos2( ) < right.a. pos2( ) );
  }
};

/// Sort by: target_id - begin_on_target (begin of aligned portion).
struct order_lookalign_TargetBegin
  : public binary_function<const look_align&, const look_align&, bool>
{
 public:
  bool operator() ( const look_align &left, const look_align &right ) const {
    if ( left.target_id == right.target_id )
      return ( left.a.pos2( ) < right.a. pos2( ) );
    return ( left.target_id < right.target_id );
  }
};

/// Sort by: target_id - end_on_target (end of aligned portion).
struct order_lookalign_TargetEnd
  : public binary_function<const look_align&, const look_align&, bool>
{
 public:
  bool operator() ( const look_align &left, const look_align &right ) const {
    if ( left.target_id == right.target_id )
      return ( left.a.Pos2( ) < right.a. Pos2( ) );
    return ( left.target_id < right.target_id );
  }
};

/// Sort by: target_id - actual_begin_on_target (the difference between
///  ActualTargetBegin and TargetBegin is that the former looks at the
///  beginning of the read, the latter looks at the beginning of the
///  aligned part).
struct order_lookalign_ActualTargetBegin
  : public binary_function<const look_align&, const look_align&, bool>
{
 public:
  bool operator() ( const look_align &left, const look_align &right ) const {
    if ( left.target_id == right.target_id )
      return ( left.a.pos2( ) - left.a.pos1( ) < right.a.pos2( ) - right.a.pos1( ) );
    return ( left.target_id < right.target_id );
  }
};

/// Sort by: target_id - if the same, sort by query_id
struct order_lookalign_by_target_id
  : public binary_function<const look_align&, const look_align&, bool>
{
 public:
  bool operator() ( const look_align &left, const look_align &right ) const {
    if ( left.target_id == right.target_id )
      return ( left.query_id < right.query_id );
    return ( left.target_id < right.target_id );
  }
};


/// Sort by: target_id - if the same, sort by query_id, if same, sort by start,end on target
struct order_lookalign_TargetQueryStartEnd
  : public binary_function<const look_align&, const look_align&, bool>
{
 public:
  bool operator() ( const look_align &left, const look_align &right ) const {
    if ( left.target_id < right.target_id )
      return true;
    if ( left.target_id > right.target_id )
      return false;
    if ( left.query_id < right.query_id )
      return true;
    if ( left.query_id > right.query_id )
      return false;
    if ( left.a.Pos1() < right.a.Pos1() )
      return true;
    if ( left.a.Pos1() > right.a.Pos1() )
      return false;
    if ( left.a.Pos2() < right.a.Pos2() )
      return true;
    if ( left.a.Pos2() > right.a.Pos2() )
      return false;
    if ( left.mutations < right.mutations )
      return true;
    if ( left.mutations > right.mutations )
      return false;
    if ( left.indels < right.indels )
      return true;

    return false;
  }
};

struct equal_lookalign_TargetQueryStartEnd
  : public binary_function<const look_align&, const look_align&, bool>
{
 public:
  bool operator() ( const look_align &left, const look_align &right ) const {
    if ( left.target_id == right.target_id &&
	 left.query_id == right.query_id &&
	 left.a.Pos2() == right.a.Pos2() &&
	 left.a.Pos1() == right.a.Pos1() &&
	 left.mutations == right.mutations &&
	 left.indels == right.indels )
      return true;

    return false;
  }
};



/**
 * LookAlignOffset
 *
 * Let say that hit is the alignment of a chimp read on human. This
 * incarnation of LookAlignOffset returns the position of the start of
 * the read on human, defined by anchoring the read by using the best
 * aligning block. The definition of "best" takes into account both
 * the length of the aligning block, and the number of mutations in
 * the block.
 */
int LookAlignOffset( const look_align_plus &hit );

/**
 * LookAlignOffset
 *
 * By example: say that two chimp reads align an human chromosome as
 * per hit1 and hit2 . LookAlignOffset fills offset with the offset
 * between the two reads determined by the best shared aligning
 * block. It returns false if they do not share a (decent) common
 * block.
 */
bool LookAlignOffset( const look_align_plus &hit1,
		      const look_align_plus &hit2,
		      int &offset );

/**
 * LookAlignOffsetOverlap
 *
 * Returns the overlap between two reads implied by the LookAlignOffset above.
 */
int LookAlignOffsetOverlap( const look_align_plus &hit1,
			    const look_align_plus &hit2 );

/// Sort by: target_id - offset.
struct order_lookalign_TargetLookAlignOffset
  : public binary_function<const look_align&,
			   const look_align&, bool>
{
public:
  bool operator() ( const look_align &left,
		    const look_align &right ) const {
    if ( left.target_id == right.target_id )
      return ( LookAlignOffset( left ) < LookAlignOffset( right ) );
    return ( left.target_id < right.target_id );
  }
};

/// Sort by: query_id
struct order_lookalign_Query
  : public binary_function<const look_align&, const look_align&, bool>
{
 public:
  bool operator() ( const look_align &left,
		    const look_align &right ) const {
    return ( left.query_id < right.query_id );
  }
};

/// Sort by: query_id, number of indels, number of mutations.
struct order_lookalign_QueryIndelsMutations :
  public binary_function<const look_align&, const look_align&, bool>
{
public:
  bool operator( ) ( const look_align &left,
		     const look_align &right ) const {
    if ( left.query_id == right.query_id ) {
      if ( left.indels == right.indels ) {
	return left.mutations < right.mutations; }
      return ( left.indels < right.indels ); }
    return ( left.query_id < right.query_id ); }
};

/// Sort by: query_id, number of errors (as from method Errors( )).
struct order_lookalign_QueryErrors :
  public binary_function<const look_align&, const look_align&, bool>
{
public:
  bool operator() ( const look_align &left,
		    const look_align &right ) const {
    if ( left.query_id == right.query_id )
      return ( left.Errors( ) < right.Errors( ) );
    return ( left.query_id < right.query_id );
  }
};

/// Sort by: number of errors (as from method ErrorRate( )).
struct order_lookalign_ErrorRate :
  public binary_function<const look_align&, const look_align&, bool>
{
public:
  bool operator() ( const look_align &left,
		    const look_align &right ) const {
    return ( left.ErrorRate( ) < right.ErrorRate( ) );
  }
};


/// Sort by: number of mutations (as from method MutationRate( )).
struct order_lookalign_MutationRate :
  public binary_function<const look_align&, const look_align&, bool>
{
public:
  bool operator() ( const look_align &left,
		    const look_align &right ) const {
    return ( left.MutationRate( ) < right.MutationRate( ) );
  }
};

/// Orders look aligns by 1) query id, then 2) complete position/structure of the target region(s)
/// the query aligns to: first by target id, then start on target, then end on target, then number
/// of gaps, then starts/ends of all blocks, and finally by number of mismatches.
struct order_lookalign_QueryIdFullTargetMutations : public binary_function<const look_align &, const look_align &, bool> {
  bool operator()(const look_align & a, const look_align & b) const {
    if ( a.QueryId() < b.QueryId() ) return true;
    if ( a.QueryId() > b.QueryId() ) return false;
    if ( a.TargetId() < b.TargetId() ) return true;
    if ( a.TargetId() > b.TargetId() ) return false;
    if ( a.StartOnTarget() < b.StartOnTarget() ) return true;
    if ( a.StartOnTarget() > b.StartOnTarget() ) return false;
    // for alignments with indels EndOnTraget requires some cycles to do the math; do it only once:
    longlong d = a.EndOnTarget() - b.EndOnTarget();
    if ( d < 0 ) return true;
    if ( d > 0 ) return false;
    // ok, these are alignments for the same read (query id) to the same contig
    // (target id) and their spans (Start,End) are exactly
    // the same; now let's try to order by number and positions of individual
    // continuous alignment blocks (if indels are present) and finally by the number
    // of mismatches:
    if ( a.indels < b.indels ) return true;
    if ( a.indels > b.indels ) return false;
    if ( a.a.Lengths(0) < b.a.Lengths(0) ) return true;
    if ( a.a.Lengths(0) > b.a.Lengths(0) ) return false;
    for ( unsigned int i = 1 ; i < static_cast<unsigned int>(a.a.Nblocks()) ; i++ ) {
      if ( a.a.Gaps(i) < b.a.Gaps(i) ) return true;
      if ( a.a.Gaps(i) > b.a.Gaps(i) ) return false;
      if ( a.a.Lengths(i) < b.a.Lengths(i) ) return true;
      if ( a.a.Lengths(i) > b.a.Lengths(i) ) return false;
    }
    if ( a.mutations < b.mutations ) return true;
    if ( a.mutations > b.mutations ) return false;
    return false;
  }
};

/// counterpart to order_lookalign_QueryIdFullTargetMutations
struct equal_lookalign_QueryIdFullTargetMutations : public binary_function<const look_align &, const look_align &, bool> {
  bool operator()(const look_align & a, const look_align & b) const {
    if ( a.QueryId() != b.QueryId() ) return false;
    if ( a.TargetId() != b.TargetId() ) return false;
    if ( a.mutations != b.mutations ) return false;
    if ( a.indels != b.indels ) return false;
    if ( a.StartOnTarget() != b.StartOnTarget() ) return false;
    if ( a.a.Lengths(0) != b.a.Lengths(0) ) return false;
    for ( unsigned int i = 1 ; i < static_cast<unsigned int>(a.a.Nblocks()) ; i++ ) {
      if ( a.a.Gaps(i) != b.a.Gaps(i) ) return false;
      if ( a.a.Lengths(i) != b.a.Lengths(i) ) return false;
    }
    return true;
  }
};



/// Boolean for look_align errors (as from method Errors( )).
struct gt_lookalign_Errors :
  public unary_function<const look_align&, bool>
{
private:
  int errors;
public:
  gt_lookalign_Errors(int errors) : errors(errors) { }
  bool operator() ( const look_align &l) const { return ( l.Errors( ) > errors ); }
};

// return true if the beginning of a query (left side for fw orientation)
// aligns to the target
bool IsLeft( const look_align &align );

// return true if the end of a query (right side for fw orientation)
// aligns to the target
bool IsRight( const look_align &align );

// return true if neither end of the query
// aligns to the target
bool IsCenter( const look_align &align );

// return the begin, end bases of the query which do not align to the target
pair<int,int> hangingEnd( const look_align &align );

// return the maximum hang amount (largest between head and tail)
int MaxHang( const look_align &align );

// return size of largest gap
int MaxGap( const look_align &align );

// a pair is logical iff hits belong to same target and have different orient
bool IsLogicalPair( const look_align &hit1,
		    const look_align &hit2 );

// an insert is valid if and only if it is both logical and not stretched
bool IsValidPair( const look_align &hit1,
		  const look_align &hit2,
		  const read_pairing &pair,
		  double max_stretch,
		  double max_mult );

// return observed separation (warning! It will assert if IsLogicalPair( )
//  fails on pair)
int ObservedSeparation( const look_align &hit1,
			const look_align &hit2 );

// return stretch (warning! It will assert if IsLogicalPair( ) fails on pair)
//  as_multiplier: if true return ( observed_separation / given_separation )
float Stretch ( const look_align &hit1,
		const look_align &hit2,
		const read_pairing &pair,
		bool as_multiplier = false );

///Trim a basevector and a look_align on that basevector in sync.
void Trim1Together(const basevector & b1, const basevector & b2,
		   const look_align & la, int startOn1, int len,
		   basevector & trimmedb1, look_align & trimmedla);

///Trim a basevector and a look_align_plus on that basevector in sync.
void Trim1Together(const basevector & b1, const basevector & b2,
		   const look_align & la, int startOn1, int len,
		   basevector & trimmedb1, look_align_plus & trimmedla);

///Trim the first alignment so it does not overlap with the second.
///Make this cleverer later to divide according to quality.
///If verbose is true, then show start and end points.
/// If la1 starts after la2, we swap the two look_aligns.
/// If la1 and la2 start at the same place, we zero the shorter
/// of the two alignments.
///Preconditions:
/// 1. la1 and la2 refer to the same query.
/// \todo there is a bug that appears to remove the last base in
/// the alignment when trimming a reversed complement alignment
/// (that is, the first base in the fw vector
void RemoveOverlap(look_align & la1, look_align & la2,
		   bool verbose = false);

///Create a look_align from an alignment_plus.
look_align FromAlignmentPlus(const alignment_plus & alp, int query_length,
			     int target_length);

/// GetBestAligns: pick out "the" best alignment of each read to the reference.
void GetBestAligns(
     const vecbasevector& bases,             // the reads
     const vec<look_align>& aligns,          // alignments of reads to reference
     const vec< vec<int> >& aligns_index,    // index by read of alignments
     vec<int>& best                          // index of best alignment (returned)
         );

/// Return the number of lines beginning with QUERY.
uint CountLookAligns( const String & fname );


/**
   Type: look_align_plus_vec

   A vector of alignments (<look_align_pluses>).  Can represent, for example,
   all alignments of a read to the reference.
*/
typedef vec< look_align_plus > look_align_plus_vec;

/**
   Type: vec_look_align_plus_vec

   A vector of vectors of alignments (<look_align_plus>es).  Can represent, for example,
   for each read, the alignments of that read to the reference.
*/
typedef vec< look_align_plus_vec > vec_look_align_plus_vec;

/**
   Type: read_aligns_plus_t

   For each read, the alignments of that read to the reference,
   as <look_align_pluses>.
*/
typedef vec_look_align_plus_vec read_aligns_plus_t;


/**
   Function: WriteVecLookAlignPlusVec

   Write a vector of vectors of <look_align_pluses> to a binary file.
*/
void WriteVecLookAlignPlusVec( const String& file_name,
			       const vec_look_align_plus_vec& vecLookAlignPlusVec );

/**
   Function: LoadVecLookAlignPlusVec

   Load a vector of vectors of <look_align_pluses> from a binary file.
*/
void LoadVecLookAlignPlusVec( const String& file_name,
			      vec_look_align_plus_vec& vecLookAlignPlusVec );








ofstream & operator << (ofstream & os, const look_align & la);

void SetWritePrettyLookAligns(ofstream & os);






#endif
