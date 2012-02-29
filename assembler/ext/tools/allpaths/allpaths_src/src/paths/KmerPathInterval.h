// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

#ifndef PATHS_KMERPATHINTERVAL_H
#define PATHS_KMERPATHINTERVAL_H

#include "CoreTools.h"
#include "math/Functions.h"
#include "CommonSemanticTypes.h"

// Portability note: endianness
// The implementations here would need to be changed for big endian architectures.

#ifndef Little_Endian
     #error KmerPathInterval.h is designed for little endian architectures.
#endif


/*
   Type: kmer_id_t 
   
   A number identifying a particular kmer.
   
   Define constants and functions involved in kmer numbering.  A valid numbering
   of kmers must satisfy the following:
    - Palindromic kmers are numbered (arbitrarily) starting with <first_palindrome>,
   and there are at most max_palindromes of them.
    - Non-palindromic kmers must have numbers in [0, first_palindrome - 1).
    - If x is the number of a non-palindromic kmer, then the number of its
   reverse complement is flip_kmer(x).

   *NOTE*: The numbering of kmers depends not on just their sequence but also
   on the total set of kmers.  Kmer numbers are assigned by <ReadsToPaths>,
   based on a specific set of kmers.  Kmer numbers assigned based on different
   sets of kmers have no relation to each other.

   Computing with kmer numbers, rather than actual kmer sequences, is useful
   for several reasons:

      - kmer numbers are more compact
      - for many purposes you only care whether two kmers are adjacent or not --
        not what their actual sequences are
      - you can represent very long sequences compactly, if the kmer numbers
      are assigned in such a way that long stretches of the sequence have
      kmer numbers in a contiguous interval; then the whole interval can be
      represented by its start & stop kmer numbers, regardless of the
      interval's length.
*/
SemanticTypeStd( longlong, kmer_id_t );

const kmer_id_t NULL_KMER_ID = -1;

const kmer_id_t max_palindromes = 1000 * 1000 * 1000; // arbitrary
const kmer_id_t max_in_five_bytes = ( kmer_id_t(1)<<40 ) - 1;    // 2^40 - 1

// Constant: first_palindrome
// Splits the space of <kmer numbers> into two halves: kmer numbers below this
// denote non-palindromic kmers (those not equal to their own reverse complement),
// and kmer numbers starting with this denote palindromic kmers.
const kmer_id_t first_palindrome = max_in_five_bytes - max_palindromes + 1;
const kmer_id_t halfway_kmer = first_palindrome/2 - 2;
inline Bool is_palindrome( kmer_id_t x ) { return x >= first_palindrome; }

// Function: flip_kmer
// Given a <kmer number> of a _non-palindromic_ kmer, give the kmer number of its reverse complement.
// Do not call this routine directly -- call <reverse_kmer()> instead, which handles both palindromic
// and non-palindromic kmers!
inline kmer_id_t flip_kmer( kmer_id_t x ) 
{    return 2*halfway_kmer - x + 1;    }

// Function: reverse_kmer
// Given a <kmer number>, give the kmer number of its reverse complement
// (which is the kmer itself, if the kmer is a palindrome).
inline kmer_id_t reverse_kmer( kmer_id_t x ) 
{    if ( is_palindrome(x) ) return x;
     return flip_kmer(x);    }

class tagged_rpint;
class big_tagged_rpint;
class new_tagged_rpint;

/*
   Class: KmerPathInterval

   A KmerPathInterval is an object which represents either of the two types 
   of segments which make up a <kmer path>: a range of <kmer numbers> or a gap.
   They are both just intervals of nonnegative integers.
   For ALLPATHS, gaps are not used and a KmerPathInterval always represents
   a range of kmer numbers.
  
   The data is stored in 8 bytes: 5 bytes for interval start,
   one bit for the interval-vs-gap flag, and the rest for the
   interval length.  For now the length is capped at 65536, since
   there is only room for a 2-byte length in the current tagged_rpint.
   But if we need longer intervals later, this class can go up to 2^23
   with the trivial change noted below.  (But other classes will die.)

   See also: <path interval database>

   Notes:
   
     - a gap is measured in K-mers, not in bases.  This is needed
   to ensure that the gap length represented is never negative!!
   But to interact with the user, we should always represent gaps in
   sequence space, not K-mer space.

     - gaps in kmer paths are not used for assembly from short reads

     - gaps in kmer paths are not to be confused with <gapped kmers>.

     - a <kmer path> is a sequence of <kmer path intervals>, but the
       interval boundaries are not reflective of any underlying genome
       features -- they a only reflective of our particular choice
       of kmer numbering.

   To do: check that "gaps in kmer paths are not used for assembly from short reads"
   is true.
*/
class KmerPathInterval {
 public:
  KmerPathInterval( ) { }

  KmerPathInterval( const kmer_id_t start, const kmer_id_t stop, const Bool is_gap = False ) 
  { Set( start, stop, is_gap ); }

  // Restriction to enforce compatibility with tagged_rpint and big_tagged_rpint.
  static const int maxLength = 0xffff;

 private:
  // We have to artificially limit the maximum difference (i.e. what
  // we actually store) in order to maintain the limit on Length().
  // The -1 is necessary because of the +1 in the Length() function.
  static const int maxDiff = maxLength - 1;

 public:
  // large integer constants specified by bit shifting to keep i686 happy
  void Set( const kmer_id_t start, const kmer_id_t stop, const Bool is_gap = False ) {
    AssertGe( start, 0 );
    AssertLe( start, ((kmer_id_t)1<<40)-1 );  /* 40 bits, 0xffffffffff */
    AssertLe( start, stop );
    // Changing this bound to 0x7fffff would allow intervals up to 8388607
    AssertLe( stop - start, maxDiff );     /* 16 bits, 65535 */
    // pack it in:
    data_ =  (start)  |  ( (stop-start) << 40 )  
                      |  ( (ulonglong)(is_gap==True) << 63 );
  }

  kmer_id_t Start( )  const { return(start_( )); }
  kmer_id_t Stop( )  const { return(start_( ) + diff_( )); }
  int Length( ) const { return((int)(diff_( ) + 1)); }
  // Same functions with different names, for thinking about gaps:
  kmer_id_t Minimum( )  const { return(start_( )); }
  kmer_id_t Maximum( )  const { return(start_( ) + diff_( )); }
  int Stretch( ) const { return((int)(diff_())); } // Note that this is different than Length()!

  Bool isGap( ) const { return ( (data_ & ((kmer_id_t)1<<63)) != 0 ); }
  Bool isSeq( ) const { return ( (data_ & ((kmer_id_t)1<<63)) == 0 ); }

  Bool Contains( const kmer_id_t kmer ) const 
    { return( Start() <= kmer && kmer <= Stop() ); }
  Bool Overlaps( const KmerPathInterval& other ) const
    { return( max(Start(),other.Start()) <= min(Stop(),other.Stop()) ); }

  template<class TAG> 
  void AppendToDatabase( vec<TAG>& segs, int i, int j ) const;

  inline friend Bool operator==( const KmerPathInterval& i1, const KmerPathInterval& i2 )
    { return i1.data_ == i2.data_; }
  inline friend Bool operator!=( const KmerPathInterval& i1, const KmerPathInterval& i2 )
    { return i1.data_ != i2.data_; }
  inline friend Bool operator<( const KmerPathInterval& i1, const KmerPathInterval& i2 )
    { return i1.data_ < i2.data_; }

  // Pretty-printed output:
  friend ostream& operator<<(ostream& out, const KmerPathInterval& rpi) {
    if( rpi.isSeq() )
      out << "[" << setw(12) << rpi.Start() << "-" << setw(12) << rpi.Stop() << "]";
    else
      out << "      (gap " 
          << setw(3) << rpi.Minimum() 
          << " - " 
          << std::left << setw(3) << rpi.Maximum() 
          << std::right << ")      ";
    return out;
  }

  // A class that is output as as many spaces as a KmerPathInterval takes up.
  struct Blank {
    friend ostream& operator<<(ostream& out, const Blank& s ) {
      return out << "                           ";
    }
  };

  ulonglong GetHash() const { return data_; }

  struct OrderByStart : public binary_function<KmerPathInterval,KmerPathInterval,bool> {
    bool operator() ( const KmerPathInterval& i1, const KmerPathInterval& i2 ) const
    {    
      return i1.Start( ) < i2.Start( );
    }
  };

 private:
  ulonglong data_;
  inline kmer_id_t start_() const { return ( data_ & (((kmer_id_t)1<<40)-1) ); }
  inline kmer_id_t diff_()  const { return ( (data_ >> 40) & 0x7fffff ); }
  // This & 0x7fffff just strips off the gap bit; it does not enforce
  // the current Length < 65536 requirement.
};

TRIVIALLY_SERIALIZABLE(KmerPathInterval);

inline Bool cmp_start( const KmerPathInterval& i1, const KmerPathInterval& i2 )
{    return i1.Start( ) < i2.Start( );    }

/**
   Class: path interval database

   An indexed vector of <tagged path intervals>.  Supports quickly finding
   all intervals containing a given kmer, or intersecting a given <interval>.

   Represented as vec<TAG>, where TAG models <tagged path interval>.
   To index a path interval database for subsequent fast searching, you call
   Prepare() on a vec<TAG>.
   
   Depending on the tag type, different kinds of path interval databases
   can be created.  You can take all path intervals in read paths, and
   with each interval record the read it's from.  This can be used to implicitly
   represent the adjacency relation on kmers: when you need to know the successors
   or predecessors of a kmer, you find all read path intervals containing this kmer,
   and look to the right/left of the kmer in the read.
   Or, a path interval database can be used to represent all path intervals in unipaths.
   Then, each path interval is disjoint from all other path intervals (because each
   kmer occurs in exactly one unipath), and with each path interval we can associate
   the unipath id of the unipath containing the interval.
*/


// Semantic Type: path_interval_id_t
// Index of a KmerPathInterval in a <path interval database>.
// Returned by Contains() methods.
SemanticTypeStd( longlong, path_interval_id_t );



/**
   Type Concept: tagged path interval

   A tagged_rpint ("tagged path interval") represents a <path interval>
   taken from a <read path>, together with information on its origin (which read
   path and where in it), and a lookback entry, which facilitates searching of
   a vec<tagged_rpint>.  A tagged_rpint is packed so that it occupies 16 bytes.
  
   The origin of a read path interval is +i if it came from paths[i],
   and is -i-1 if it came from paths_rc[i].
  
   If v is a vec<tagged_rpint> which has been <Prepare()>'d and index is a k-mer 
   address, then Contains returns {i} such that
   v[i].Start( ) <= index <= v[i].Stop( ).
  
   Only the "sequence" types of read path intervals are stored in this object,
   so there is no extra bit to distinguish the "gap" types of intervals.

   Modeled by <tagged_rpint>, <big_tagged_rpint>, <new_tagged_rpint>.
*/

/**
   Class: tagged_rpint
  
   A tagged_rpint ("tagged read path interval") represents a read path interval,
   taken from a <read path>, together with information on its origin (which read
   path and where in it), and a lookback entry, which facilitates searching of
   a vec<tagged_rpint>.  A tagged_rpint is packed so that it occupies 16 bytes.
  
   The origin of a read path interval is +i if it came from paths[i],
   and is -i-1 if it came from paths_rc[i].
  
   If v is a vec<tagged_rpint> which has been <Prepare()>'d and index is a k-mer 
   address, then Contains returns {i} such that
   v[i].Start( ) <= index <= v[i].Stop( ).
  
   Only the "sequence" types of read path intervals are stored in this object,
   so there is no extra bit to distinguish the "gap" types of intervals.

   Models <tagged read path interval>.
*/
// TODO: potentially dangerous truncation of index by PathId, ReadId
class tagged_rpint {

     public: /* PUBLIC METHODS */

     tagged_rpint( ) : data1_(0), data2_(0) { }

     tagged_rpint( kmer_id_t start, unsigned short length, int path_id, 
          unsigned short path_pos )
     {    data1_ = (start << 24);
          data2s_[3] = length;
          data2s_[2] = path_pos;
          data2i_[0] = path_id;    }

     void Set( kmer_id_t start, unsigned short length, int path_id, 
          unsigned short path_pos )
     {    data1_ = (start << 24);
          data2s_[3] = length;
          data2s_[2] = path_pos;
          data2i_[0] = path_id;    }

     kmer_id_t Start( ) const
     {    return data1_ >> 24;    }

     unsigned short Length( ) const
     {    return data2s_[3];    }

     kmer_id_t Stop( ) const
     {    return Start( ) + Length( ) - 1;    }

     int PathId( ) const
     {    return data2i_[0];    }

     Bool Fw( ) const { return PathId( ) >= 0; }
     Bool Rc( ) const { return PathId( ) < 0; }

     read_id_t ReadId( ) const
     {    int id = PathId( );
          return ( id >= 0 ? id : -id-1 );    }

     // Method: PathPos
     // Return the position of this KmerPathInterval within
     // its path (i.e. its index in the path's list of KmerPathIntervals).
     unsigned short PathPos( ) const
     {    return data2s_[2];    }

     Bool Overlaps( const KmerPathInterval& other ) const
     { return( max(Start(),other.Start()) <= min(Stop(),other.Stop()) ); }
  
     Bool Overlaps( const tagged_rpint& other ) const
     { return( max(Start(),other.Start()) <= min(Stop(),other.Stop()) ); }
  
     int Lookback( ) const
     {    return data1i_[0] & (16777216 - 1);    }

     friend Bool operator==( const tagged_rpint& s1, const tagged_rpint& s2 )
     {    return ( s1.Start( ) == s2.Start( ) &&
                   s1.Length( ) == s2.Length( ) &&
                   s1.PathId( ) == s2.PathId( ) &&
                   s1.PathPos( ) == s2.PathPos( ) );    }

     private: /* PRIVATE METHODS */

     void SetLookback( int lookback )
     {    data1i_[0] = ( data1i_[0] & ((512 - 1) << 24) ) ^ lookback;    }

     friend Bool operator<( const tagged_rpint& s1, const tagged_rpint& s2 )
     {    return s1.Start( ) < s2.Start( );    }

     public: /* PUBLIC VECTOR METHODS */

     template<class TAG> friend void Prepare( vec<TAG>& segs );

     // Note for Contains: if cap is set to a positive value, then no more than
     // that many values will be put in answer.

     template<class TAG>
     friend void Contains( const vec<TAG>& segs, kmer_id_t index,
          vec<longlong>& answer, bool append = false, int cap = -1 );

     template<class TAG>
     friend void Contains( const vec<TAG>& segs, KmerPathInterval rpi,
          vec<longlong>& answer, bool append = false, int cap = -1 );

     template<class TAG>
     friend longlong Instance( const vec<TAG>& segs, kmer_id_t k );


     private: /* THE DATA */

     union {
          ulonglong data1_;
          int data1i_[2];
          unsigned short data1s_[4];
     };

     union {
          ulonglong data2_;
          int data2i_[2];
          unsigned short data2s_[4];
     };

     // data1_:
     // start of k-mer segment [5 bytes]
     // lookback [3 bytes]
     //
     // data2_:
     // length of k-mer segment [2 bytes]
     // position in read path (which rpint) [2 bytes]
     // read path id [4 bytes] -- explicitly signed:
     //   paths[i] is PathID i; paths_rc[i] is PathID -i-1.
public:
  static const unsigned int LOOKBACK_MAX = 0xffffff; // 3 bytes
  static const unsigned int POSITION_MAX = 0xffff;   // 2 bytes
  static const unsigned int LENGTH_MAX   = 0xffff;   // 2 bytes

  static const Bool IS_BIG = False;
};

/**
   Class: big_tagged_rpint
   
   A big_tagged_rpint is the same as a tagged_rpint, but allows more space for
   the path_pos.  There is room to spare, so other data components could be
   expanded as well.

   For short reads, you might think that much space for read position or
   interval length will never be needed.  But even for short reads, we sometimes represent
   <genome parts> as <genome paths>, and genome parts may be arbitrarily long.
*/
// TODO: potentially dangerous truncation of index by id and pos args
class big_tagged_rpint {

     public: /* PUBLIC METHODS */

     big_tagged_rpint( ) : data1_(0), data2_(0), data3_(0) { }

     big_tagged_rpint( kmer_id_t start, unsigned short length, int path_id, 
          int path_pos )
     {    data1_ = (start << 24);
          data2s_[3] = length;
          data3i_[0] = path_pos;
          data2i_[0] = path_id;    }

     void Set( kmer_id_t start, unsigned short length, int path_id, int path_pos )
     {    data1_ = (start << 24);
          data2s_[3] = length;
          data3i_[0] = path_pos;
          data2i_[0] = path_id;    }

     kmer_id_t Start( ) const
     {    return data1_ >> 24;    }

     unsigned short Length( ) const
     {    return data2s_[3];    }

     kmer_id_t Stop( ) const
     {    return Start( ) + Length( ) - 1;    }

     int PathId( ) const
     {    return data2i_[0];    }

     Bool Fw( ) const { return PathId( ) >= 0; }
     Bool Rc( ) const { return PathId( ) < 0; }

     int ReadId( ) const
     {    int id = PathId( );
          return ( id >= 0 ? id : -id-1 );    }

     int PathPos( ) const
     {    return data3i_[0];    }

     Bool Overlaps( const KmerPathInterval& other ) const
     { return( max(Start(),other.Start()) <= min(Stop(),other.Stop()) ); }
  
     Bool Overlaps( const big_tagged_rpint& other ) const
     { return( max(Start(),other.Start()) <= min(Stop(),other.Stop()) ); }
  
     int Lookback( ) const
     {    return data3i_[1];    }

     friend Bool operator==( const big_tagged_rpint& s1, const big_tagged_rpint& s2 )
     {    return ( s1.Start( ) == s2.Start( ) &&
                   s1.Length( ) == s2.Length( ) &&
                   s1.PathId( ) == s2.PathId( ) &&
                   s1.PathPos( ) == s2.PathPos( ) );    }

     private: /* PRIVATE METHODS */

     void SetLookback( int lookback )
     {    data3i_[1] = lookback;    }

     friend Bool operator<( const big_tagged_rpint& s1, const big_tagged_rpint& s2 )
     {    return s1.Start( ) < s2.Start( );    }

     public: /* PUBLIC VECTOR METHODS */

     template<class TAG> friend void Prepare( vec<TAG>& segs );

     // Note for Contains: if cap is set to a positive value, then no more than
     // that many values will be put in answer.

     template<class TAG>
     friend void Contains( const vec<TAG>& segs, kmer_id_t index,
          vec<longlong>& answer, bool append = false, int cap = -1 );

     template<class TAG>
     friend void Contains( const vec<TAG>& segs, KmerPathInterval rpi,
          vec<longlong>& answer, bool append = false, int cap = -1 );

     template<class TAG>
     friend longlong Instance( const vec<TAG>& segs, kmer_id_t k );


     private: /* THE DATA */

     union {
          ulonglong data1_;
          int data1i_[2];
          unsigned short data1s_[4];
     };

     union {
          ulonglong data2_;
          int data2i_[2];
          unsigned short data2s_[4];
     };

     union {
          ulonglong data3_;
          int data3i_[2];
          unsigned short data3s_[4];
     };

     
     // data1_:
     // start of k-mer segment [5 bytes]
     // unused [3 bytes]
     //
     // data2_:
     // length of k-mer segment [2 bytes]
     // unused [2 bytes]
     // read path id [4 bytes] -- explicitly signed:
     // (paths[i] is PathID i; paths_rc[i] is PathID -i-1)
     //
     // data3_:
     // lookback [4 bytes]
     // position in read path [4 bytes]
public:
  static const unsigned int LOOKBACK_MAX = 0xffffffff; // 4 bytes
  static const unsigned int POSITION_MAX = 0xffffffff; // 4 bytes
  static const unsigned int LENGTH_MAX   = 0xffff;     // 2 bytes
  
  static const Bool IS_BIG = True;
};

// A new_tagged_rpint is the same as a big_tagged_rpint, but allows more space for
// the length.  Presumably we should change new_tagged_rpint to big_tagged_rpint,
// and get rid of the old definition.

class new_tagged_rpint {

     public: /* PUBLIC METHODS */

     new_tagged_rpint( ) : data1_(0), data2_(0), data3_(0) { }

     new_tagged_rpint( kmer_id_t start, unsigned int length, int path_id, 
          int path_pos )
     {    data1_ = (start << 24);
          data2i_[1] = length;
          data3i_[0] = path_pos;
          data2i_[0] = path_id;    }

     void Set( kmer_id_t start, unsigned int length, int path_id, int path_pos )
     {    data1_ = (start << 24);
          data2i_[1] = length;
          data3i_[0] = path_pos;
          data2i_[0] = path_id;    }

     kmer_id_t Start( ) const
     {    return data1_ >> 24;    }

     unsigned int Length( ) const
     {    return data2i_[1];    }

     kmer_id_t Stop( ) const
     {    return Start( ) + Length( ) - 1;    }

     int PathId( ) const
     {    return data2i_[0];    }

     Bool Fw( ) const { return PathId( ) >= 0; }
     Bool Rc( ) const { return PathId( ) < 0; }

     int ReadId( ) const
     {    int id = PathId( );
          return ( id >= 0 ? id : -id-1 );    }

     int PathPos( ) const
     {    return data3i_[0];    }

     Bool Overlaps( const KmerPathInterval& other ) const
     { return( max(Start(),other.Start()) <= min(Stop(),other.Stop()) ); }
  
     Bool Overlaps( const new_tagged_rpint& other ) const
     { return( max(Start(),other.Start()) <= min(Stop(),other.Stop()) ); }
  
     int Lookback( ) const
     {    return data3i_[1];    }


     private: /* PRIVATE METHODS */

     void SetLookback( int lookback )
     {    data3i_[1] = lookback;    }

     friend Bool operator<( const new_tagged_rpint& s1, const new_tagged_rpint& s2 )
     {    return s1.Start( ) < s2.Start( );    }

     friend Bool operator==( const new_tagged_rpint& s1, const new_tagged_rpint& s2 )
     {    return ( s1.Start( ) == s2.Start( ) &&
                   s1.Length( ) == s2.Length( ) &&
                   s1.PathId( ) == s2.PathId( ) &&
                   s1.PathPos( ) == s2.PathPos( ) );    }

     public: /* PUBLIC VECTOR METHODS */

     template<class TAG> friend void Prepare( vec<TAG>& segs );

     // Note for Contains: if cap is set to a positive value, then no more than
     // that many values will be put in answer.

     template<class TAG>
     friend void Contains( const vec<TAG>& segs, kmer_id_t index,
          vec<longlong>& answer, bool append = false, int cap = -1 );

     template<class TAG>
     friend void Contains( const vec<TAG>& segs, KmerPathInterval rpi,
          vec<longlong>& answer, bool append = false, int cap = -1 );

     template<class TAG>
     friend kmer_id_t Instance( const vec<TAG>& segs, kmer_id_t k );

     private: /* THE DATA */

     union {
          ulonglong data1_;
          int data1i_[2];
          unsigned short data1s_[4];
     };

     union {
          ulonglong data2_;
          int data2i_[2];
          unsigned short data2s_[4];
     };

     union {
          ulonglong data3_;
          int data3i_[2];
          unsigned short data3s_[4];
     };

     // data1_:
     // start of k-mer segment [5 bytes]
     // unused [3 bytes]
     //
     // data2_:
     // length of k-mer segment [4 bytes]
     // read path id [4 bytes] -- explicitly signed:
     // (paths[i] is PathID i; paths_rc[i] is PathID -i-1)
     //
     // data3_:
     // lookback [4 bytes]
     // position in read path [4 bytes]
public:
  static const unsigned int LOOKBACK_MAX = 0xffffffff; // 4 bytes
  static const unsigned int POSITION_MAX = 0xffffffff; // 4 bytes
  static const unsigned int LENGTH_MAX   = 0xffffffff; // 4 bytes
};


#endif

// Synonyms: Various synonyms
//    kmer id - See <kmer_id_t>
//    kmer number - See <kmer_id_t>
//    kmer index - See <kmer_id_t>
//    k-mer index - See <kmer_id_t>
//    path interval - See <KmerPathInterval>
//    interval - See <KmerPathInterval>
//    kmer range - See <KmerPathInterval>



