///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef KMERRECORD
#define KMERRECORD

#include "Basevector.h"
#include "BasevectorTools.h"
#include "CoreTools.h"
#include "kmers/KmerShape.h"
#include "CommonSemanticTypes.h"
#include "math/Functions.h"
#include "math/HoInterval.h"

/**
   File: KmerRecord.h

   Data structures for representing kmers and kmer occurrences,
   together with some associated information (for example the
   frequency of the kmer or the location (read and position) of the
   kmer occurrence).
*/


// Table for converting bytes of data from basevectors into actual bases.
// Basically, it swaps the first two bits in the byte with the last two bits,
// and also swaps bits 3/4 with bits 5/6.
//
// This is used in comparing basevectors alphanumerically:
// For example, if byte X represents the bases "ACAT" and Y represents "TGGA",
// then X is stored as '11-00-01-00' = 0xc4, while Y = '00-10-10-11' = 0x2b.
// But bv_table[X] = 0x13 and bv_table[Y] = 0xe4, so bv_table[X] < bv_table[Y].
//
// Idea stolen from basevector::RCtable.
static const unsigned char bv_table[256] = {
  0x00, 0x40, 0x80, 0xc0, 0x10, 0x50, 0x90, 0xd0,
  0x20, 0x60, 0xa0, 0xe0, 0x30, 0x70, 0xb0, 0xf0,
  0x04, 0x44, 0x84, 0xc4, 0x14, 0x54, 0x94, 0xd4,
  0x24, 0x64, 0xa4, 0xe4, 0x34, 0x74, 0xb4, 0xf4,
  0x08, 0x48, 0x88, 0xc8, 0x18, 0x58, 0x98, 0xd8,
  0x28, 0x68, 0xa8, 0xe8, 0x38, 0x78, 0xb8, 0xf8,
  0x0c, 0x4c, 0x8c, 0xcc, 0x1c, 0x5c, 0x9c, 0xdc,
  0x2c, 0x6c, 0xac, 0xec, 0x3c, 0x7c, 0xbc, 0xfc,
  0x01, 0x41, 0x81, 0xc1, 0x11, 0x51, 0x91, 0xd1,
  0x21, 0x61, 0xa1, 0xe1, 0x31, 0x71, 0xb1, 0xf1,
  0x05, 0x45, 0x85, 0xc5, 0x15, 0x55, 0x95, 0xd5,
  0x25, 0x65, 0xa5, 0xe5, 0x35, 0x75, 0xb5, 0xf5,
  0x09, 0x49, 0x89, 0xc9, 0x19, 0x59, 0x99, 0xd9,
  0x29, 0x69, 0xa9, 0xe9, 0x39, 0x79, 0xb9, 0xf9,
  0x0d, 0x4d, 0x8d, 0xcd, 0x1d, 0x5d, 0x9d, 0xdd,
  0x2d, 0x6d, 0xad, 0xed, 0x3d, 0x7d, 0xbd, 0xfd,
  0x02, 0x42, 0x82, 0xc2, 0x12, 0x52, 0x92, 0xd2,
  0x22, 0x62, 0xa2, 0xe2, 0x32, 0x72, 0xb2, 0xf2,
  0x06, 0x46, 0x86, 0xc6, 0x16, 0x56, 0x96, 0xd6,
  0x26, 0x66, 0xa6, 0xe6, 0x36, 0x76, 0xb6, 0xf6,
  0x0a, 0x4a, 0x8a, 0xca, 0x1a, 0x5a, 0x9a, 0xda,
  0x2a, 0x6a, 0xaa, 0xea, 0x3a, 0x7a, 0xba, 0xfa,
  0x0e, 0x4e, 0x8e, 0xce, 0x1e, 0x5e, 0x9e, 0xde,
  0x2e, 0x6e, 0xae, 0xee, 0x3e, 0x7e, 0xbe, 0xfe,
  0x03, 0x43, 0x83, 0xc3, 0x13, 0x53, 0x93, 0xd3,
  0x23, 0x63, 0xa3, 0xe3, 0x33, 0x73, 0xb3, 0xf3,
  0x07, 0x47, 0x87, 0xc7, 0x17, 0x57, 0x97, 0xd7,
  0x27, 0x67, 0xa7, 0xe7, 0x37, 0x77, 0xb7, 0xf7,
  0x0b, 0x4b, 0x8b, 0xcb, 0x1b, 0x5b, 0x9b, 0xdb,
  0x2b, 0x6b, 0xab, 0xeb, 0x3b, 0x7b, 0xbb, 0xfb,
  0x0f, 0x4f, 0x8f, 0xcf, 0x1f, 0x5f, 0x9f, 0xdf,
  0x2f, 0x6f, 0xaf, 0xef, 0x3f, 0x7f, 0xbf, 0xff
};


// Convert a four-byte unsigned int into a different unsigned int, by flipping each
// of the bytes independently according to bv_table.
//
// Idea stolen from basevector::ReverseComplement(unsigned int).
inline unsigned int bv_table_convert( unsigned int word4 ) {
  unsigned int ret = 0;
  ret = bv_table[ word4 & 0xFF ]; // 1st byte
  ret <<= 8;
  word4 >>= 8;
  ret |= bv_table[ word4 & 0xFF ]; // 2nd byte
  ret <<= 8;
  word4 >>= 8;
  ret |= bv_table[ word4 & 0xFF ]; // 3rd byte
  ret <<= 8;
  word4 >>= 8;
  ret |= bv_table[ word4 & 0xFF ]; // 4th byte
  return ret;
}

template < int Kmod4 >
struct ZeroOutUnusedBits {
  static void ZeroUnusedBits( unsigned char &x ) { STATIC_ASSERT_M( 0, UnknownKmod4 ); }
};

template <>
struct ZeroOutUnusedBits<0> {
  static void ZeroUnusedBits( unsigned char &x ) {  }
};

template <>
struct ZeroOutUnusedBits<1> {
  static void ZeroUnusedBits( unsigned char &x ) {
    x &= 3;
  }
};

template <>
struct ZeroOutUnusedBits<2> {
  static void ZeroUnusedBits( unsigned char &x ) {
    x &= 0xf;
  }
};

template <>
struct ZeroOutUnusedBits<3> {
  static void ZeroUnusedBits( unsigned char &x ) {
    x &= 0x3f;
  }
};

/// Class: byte_pac
///
/// A byte_pac is a structure with A+B bytes of data, stored using ints, having
/// a < operator which facilitates sorting using only the first A bytes of
/// data.  This is needed in SortKmers() so that we can group together
/// SortKmerOutputRecords that record different occurrences of the same kmer,
/// by sorting a vec of such records on kmer sequence only (so that
/// records recording different occurrences of the same kmer appear
/// equal to the sorting routine).
/// Slightly more than A bytes will be used for sorting if it is more efficient.
/// operator== compares using all data, while EqualKmers compares using
/// only the kmer data (the first A bytes).

template< int A, int B >
class byte_pac {
public:
  static const int A_VAL = A;
  static const int B_VAL = B;

  static const int NUM_INTS = (A+B+3)/4;
  static const int NUM_SHORTS = NUM_INTS * 2;
  static const int NUM_BYTES = NUM_INTS * 4;
  union
  {
    unsigned int ints[ NUM_INTS ];
    unsigned short shorts[ NUM_SHORTS ];
    unsigned char bytes[ NUM_BYTES ];
  };

  byte_pac() {
    memset( bytes, 0, NUM_BYTES );
  }

  unsigned char * begin() { return bytes; }
  const unsigned char * begin() const { return bytes; }
  unsigned char * end() { return bytes + size(); }
  const unsigned char * end() const { return bytes + size(); }
  unsigned int size() const { return  NUM_INTS * 4; }

  ///Compare only the bytes that contain kmer data.

  bool EqualKmers(const byte_pac & b) const {
    return 0==memcmp(bytes,b.bytes,A);
  }

  int CmpKmers(const byte_pac & b) const {
    if ( *this < b ) return -1;
    if ( b < *this ) return 1;
    return 0;
  }

  ///Compare only the bytes that contain data.
  friend bool operator==(const byte_pac & l, const byte_pac & r) {
    return 0==memcmp(l.bytes,r.bytes,A+B);
  }

};  // class byte_pac

template < int A, int B>
void SetKmer( byte_pac<A,B>& bpac, const basevector& b, nbases_t K ) {
    int end = K & ~15;
    unsigned int* dst = bpac.ints;
    for ( int idx = 0; idx < end; idx += 16 )
        *dst++ = b.extractKmer(idx,16);
    if ( end < K ) *dst = b.extractKmer(end,K-end);
}

template < int A, int B, nbases_t K>
void SetKmer( byte_pac<A,B>& bpac, const unsigned int *b ) {
  memset( bpac.bytes, 0, A );
    for ( int j = 0; j < (K+15)/16; j++ )
      bpac.ints[j] = b[j];
    ZeroOutUnusedBits< K % 4 >::ZeroUnusedBits( bpac.bytes[ A-1 ] );
}

// K=4
template< int B >
bool operator< ( const byte_pac<1,B>& b1, const byte_pac<1,B>& b2 )
{
  return b1.bytes[0] < b2.bytes[0];
}

// K=8
template< int B >
bool operator< ( const byte_pac<2,B>& b1, const byte_pac<2,B>& b2 )
{
  if ( b1.shorts[0] < b2.shorts[0] ) return true;
  return false;
}

// K=12
template< int B >
bool operator< ( const byte_pac<3,B>& b1, const byte_pac<3,B>& b2 )
{
  if ( b1.shorts[0] < b2.shorts[0] ) return true;
  if ( b1.shorts[0] > b2.shorts[0] ) return false;
  if ( b1.bytes[2] < b2.bytes[2] ) return true;
  return false;
}

// K=16
template< int B >
bool operator< ( const byte_pac<4,B>& b1, const byte_pac<4,B>& b2 )
{
  if ( b1.ints[0] < b2.ints[0] ) return true;
  return false;
}

// K=20
template< int B >
bool operator< ( const byte_pac<5,B>& b1, const byte_pac<5,B>& b2 )
{
  if ( b1.ints[0] < b2.ints[0] ) return true;
  if ( b1.ints[0] > b2.ints[0] ) return false;
  if ( b1.bytes[4] < b2.bytes[4] ) return true;
  return false;
}

// K=24
template< int B >
bool operator< ( const byte_pac<6,B>& b1, const byte_pac<6,B>& b2 )
{
  if ( b1.ints[0] < b2.ints[0] ) return true;
  if ( b1.ints[0] > b2.ints[0] ) return false;
  if ( b1.shorts[2] < b2.shorts[2] ) return true;
  return false;
}

// Warning: the following definition will only work correctly if K is divisible by 4!

#define DEFINE_BYTE_PAC_LESS(K)                                               \
     template< int B >                                                        \
     bool operator< ( const byte_pac<K/4,B>& b1, const byte_pac<K/4,B>& b2 )  \
     {    return memcmp( &b1.ints[0], &b2.ints[0], K/4 ) < 0;    }

DEFINE_BYTE_PAC_LESS(28)
DEFINE_BYTE_PAC_LESS(32)
DEFINE_BYTE_PAC_LESS(36)
DEFINE_BYTE_PAC_LESS(40)
DEFINE_BYTE_PAC_LESS(44)
DEFINE_BYTE_PAC_LESS(48)
DEFINE_BYTE_PAC_LESS(52)
DEFINE_BYTE_PAC_LESS(60)
DEFINE_BYTE_PAC_LESS(64)
DEFINE_BYTE_PAC_LESS(68)
DEFINE_BYTE_PAC_LESS(80)
DEFINE_BYTE_PAC_LESS(88)
DEFINE_BYTE_PAC_LESS(96)
DEFINE_BYTE_PAC_LESS(100)
DEFINE_BYTE_PAC_LESS(128)
DEFINE_BYTE_PAC_LESS(144)
DEFINE_BYTE_PAC_LESS(192)
DEFINE_BYTE_PAC_LESS(200)
DEFINE_BYTE_PAC_LESS(320)
DEFINE_BYTE_PAC_LESS(368)
DEFINE_BYTE_PAC_LESS(400)
DEFINE_BYTE_PAC_LESS(500)
DEFINE_BYTE_PAC_LESS(544)
DEFINE_BYTE_PAC_LESS(640)
DEFINE_BYTE_PAC_LESS(720)
DEFINE_BYTE_PAC_LESS(1000)
DEFINE_BYTE_PAC_LESS(1200)
DEFINE_BYTE_PAC_LESS(1600)
DEFINE_BYTE_PAC_LESS(2000)
DEFINE_BYTE_PAC_LESS(4000)
DEFINE_BYTE_PAC_LESS(10000)

/// Class: kmer_record
///
/// Stores a kmer together with its origin and position.
///
/// Class kmer_record comes in small and big versions, depending on the value
/// of the template parameter I (1 or 2).
///
/// A kmer_record consists of the following
///   - k bases, stored 4 per byte. k must be a multiple of 4.
///   - a 4-byte readID
///   - a position, negated in the case where the reverse complement of
///   the k-mer occurs at the given position.  The position occupies 2 bytes
///   when I = 1 and 4 bytes when I = 2.
///
/// Models <SortKmersOutputRecord>.  See also <kmer>, which stores just
/// the kmer without the origin and position.
template<int K, int I = 1> class kmer_record {
private:

  // byte_pac< (K+3)/4, 5 + 3*I > data_;
  typedef byte_pac< (K+3)/4, 4 + 2*I > data_t;
  data_t data_;

public:

  /// Constructor is empty and does no work after compilation.
  /// Here only for the static assert: if we use kmers of size that
  /// is not a multiple of 4, our comparison and equality operators
  /// can mess up. So let's make sure we don't.
  kmer_record() { /*STATIC_ASSERT_M(0==K%4, K_not_multiple_4);*/ }

  /// Size of the bases contained in this kmer_record.
  /// For use when creating a basevector from our data.
  static const int BASES_SIZE=K;

  // number of 2-byte words to represent a position
  static const int POSITION_WORDS=I;

  const unsigned char* Bytes( ) const
  {    return data_.bytes;    }


  const unsigned char * BytesEnd() const {
    return data_.end();
  }

  int BytesSize() const {
    return data_.size();
  }

  unsigned int* Ints( )
  {    return data_.ints;    }

  void GetBasevector( basevector& kmer ) const {
      kmer.assignBaseBits(K,data_.begin());
  }

  // TODO: potentially dangerous truncation of readID
  // unfortunately, the whole scheme has to be reworked to accommodate 64-bit ids
  read_id_t GetId( ) const
  {
    if ( I == 1 )
      {
#ifdef Little_Endian
	return data_.bytes[(K+3)/4] |
	  (short(data_.bytes[(K+3)/4 + 1]) << 8) |
	  (int(data_.bytes[(K+3)/4 + 2]) << 16) |
	  (int(data_.bytes[(K+3)/4 + 3]) << 24);
#endif
#ifdef Big_Endian
	return data_.bytes[(K+3)/4 + 3] |
	  (short(data_.bytes[(K+3)/4 + 2]) << 8) |
	  (int(data_.bytes[(K+3)/4 + 1]) << 16) |
	  (int(data_.bytes[(K+3)/4]) << 24);
#endif
      }
    else if ( I == 2 ) return data_.ints[ (K+15)/16 ];
    else return 1; /* unreachable statement */    }

  read_pos_t GetPos( ) const
  {
    if ( I == 1 )
      {
#ifdef Little_Endian
	return (int) (short) (data_.bytes[(K+3)/4 + 4]
			      | (short(data_.bytes[(K+3)/4 + 5]) << 8));
#endif
#ifdef Big_Endian
	return (int) (short) (data_.bytes[(K+3)/4 + 5]
			      | (short(data_.bytes[(K+3)/4 + 4]) << 8));
#endif
      }

    else if ( I == 2 ) return data_.ints[ (K+15)/16 + 1 ];
    else return 1; /* unreachable statement */
  }

  /// Returns position of this kmer_record on the sequence (i.e. read) it was
  /// drawn from as genomic interval with id, i.e. sequence_id:[start,stop);
  /// the coordinates returned are always on forward strand on the sequence regardless
  /// of whether this kmer or its reverse complement are found there.
  HoIntervalWithId GetInterval() const {
    read_pos_t pos = GetPos(); // unpack position
    // if the kmer is forward-oriented, just use its coords to generate
    // the interval, otherwise don't forget to invert position to make it positive!
    if ( pos < 0 ) {
      pos = (-pos); // invert
      --pos;
    }
    return HoIntervalWithId( pos, pos + K, GetId() );
  }

  // TODO: potentially dangerous truncation of index by read_id and read_pos
  void Set( const basevector& b, int read_id, int read_pos ) {
    Assert( I == 1 || I == 2 );
    Assert( read_id >= 0 );
    SetKmer( data_, b, K );

    if ( I == 1 ) {
#ifdef Little_Endian
      for ( int j = 0; j < 4; j++ )
	data_.bytes[ (K+3)/4 + j ] =
	  ((unsigned char*) (&read_id))[j];
#endif
#ifdef Big_Endian
      for ( int j = 0; j < 4; j++ )
	data_.bytes[ (K+3)/4 + j ] =
	  ((unsigned char*) (&read_id))[j+1];
#endif
      short rp = (short) read_pos;
      for ( int j = 0; j < 2; j++ ) {
	data_.bytes[ (K+3)/4 + 4 + j ] = ((unsigned char*) (&rp))[j];
      }
    }
    else if ( I == 2 ) {
      data_.ints[ (K+15)/16 ] = read_id;
      data_.ints[ (K+15)/16 + 1 ] = read_pos;
    }
  }

  // Set a kmer_record from an array of uints: -dnave 2001/10/30
  void Set( const unsigned int* b, int read_id, int read_pos ) {
    Assert( I == 1 || I == 2 );
    Assert( read_id >= 0 );
    SetKmer< data_t::A_VAL, data_t::B_VAL, K > ( data_, b );

    if ( I == 1 ) {
#ifdef Little_Endian
      for ( int j = 0; j < 4; j++ )
	data_.bytes[ (K+3)/4 + j ] =
	  ((unsigned char*) (&read_id))[j];
#endif
#ifdef Big_Endian
      for ( int j = 0; j < 4; j++ )
	data_.bytes[ (K+3)/4 + j ] =
	  ((unsigned char*) (&read_id))[j+1];
#endif
      short rp = (short) read_pos;
      for ( int j = 0; j < 2; j++ ) {
	data_.bytes[ (K+3)/4 + 4 + j ] = ((unsigned char*) (&rp))[j];
      }
    }
    else if ( I == 2 ) {
      data_.ints[ (K+15)/16 ] = read_id;
      data_.ints[ (K+15)/16 + 1 ] = read_pos;
    }
  }

  ///Only works for I=2
  void SetId( int read_id ) {
    AssertEq(I, 2);
    data_.ints[ (K+15)/16 ] = read_id;
  }

  ///Only works for I=2
  void SetPos( int read_pos ) {
    AssertEq(I, 2);
    data_.ints[ (K+15)/16 + 1 ] = read_pos;
  }

  ///We mark bad by setting the id to -1
  void MarkAsBad() { SetId(-1); };

  /// Check whether id is -1.
  bool IsBad() const { return -1 == GetId(); }

  ///Set all data from an unsigned char * buffer (used when reading from file).
  void Set( const unsigned char * rawdata ) {
    copy(rawdata, rawdata + BytesSize(), data_.begin());
  }

  ///True if kmers are equal, even if GetPos() and GetId() are different.

  bool EqualKmers( const kmer_record& k2 ) const {
      return data_.EqualKmers(k2.data_);
  }

  int CmpKmers( const kmer_record& k2 ) const {
      return data_.CmpKmers(k2.data_);
  }

  ///True if kmers are rc, even if GetPos() and GetId() are different.
  bool ReverseKmers( const kmer_record& k2 ) const {
      kmer_record rck2;
      rck2.Set(bvec(K,k2.Bytes()).ReverseComplement(), 0, 0);
      return EqualKmers(rck2);
  }

  bool EqualOrReverseKmers( const kmer_record& k2 ) const {
    return (EqualKmers(k2) || ReverseKmers(k2));
  }

  ///Pick the lower of myself and my reverse complement.
  ///If I reverse myself, set pos negative as a marker, and refer it to the
  /// back of the basevector which becomes the front as we reverse.
  void Canonicalize(int length) {
    bvec bv(K,Bytes());
    if ( bv.getCanonicalForm() == bvec::rc_form )
    {
        kmer_record rc;
        // Set the position counting from the back, and add -1 to distinguish
        // 0 forward from 0 rc
        rc.Set(bv,GetId(),-(length - (GetPos() + K))  - 1);
        *this = rc;
    }
  }

  /// Looks at the position to see if it is < 0
  bool IsReversed() const { return GetPos() < 0; }

  /// Looks at the position to see if it is < 0 ( better name )
  Bool IsRc() const { return GetPos() < 0; }
  Bool IsFw() const { return !IsRc(); }

  /// Return a positive position, even if the kmer_record is reversed.
  /// For reversed kmer, position is from back of basevector!
  ///
  /// BEWARE!  This method is only valid for kmer_records produced by
  /// SimpleSortKmers.  The equivalent for kmer_records produced by
  /// SortKmers is abs(GetPos())-1, which *always* gives a position in
  /// the forward version of the basevector.
  int TruePos() const {
    int pos = GetPos();
    return pos >= 0 ? pos : -pos-1;
  }

  read_pos_t GetReadPos() const {
    return abs( GetPos() ) - 1;
  }

  friend bool operator==( const kmer_record& k1, const kmer_record& k2 ) {
    return k1.data_ == k2.data_;
  }

  friend bool operator!=( const kmer_record& k1, const kmer_record& k2 ) {
    return !(k1 == k2);
  }

  friend bool operator<( const kmer_record& k1, const kmer_record& k2 )
  {    return k1.data_ < k2.data_;    }

  friend bool operator>( const kmer_record& k1, const kmer_record& k2 )
  {    return k2 < k1;    }

  // this is actually a 'lt' not a 'cmp'
  static Bool id_cmp( const kmer_record& k1, const kmer_record& k2 )
  {    return k1.GetId( ) < k2.GetId( );    }

  // this is actually a 'lt' not a 'cmp'
  static Bool id_cmp_pos( const kmer_record& k1, const kmer_record& k2 )
  {
    read_id_t id1 = k1.GetId( ), id2 = k2.GetId( );
    if ( id1 < id2 ) return True;
    if ( id1 > id2 ) return False;
    return Abs( k1.GetPos( ) ) < Abs( k2.GetPos( ) );
  }

  // this is actually a 'lt' not a 'cmp'
  static Bool cmp_pos( const kmer_record& k1, const kmer_record& k2 )
  {
    return k1.GetPos( ) < k2.GetPos( );
  }

  // Compares these kmer_records' basevectors as basevectors, i.e., uses the
  // same sort order as Basevector::operator<.  This is NOT the same as the
  // "CmpKmers" sort order used elsewhere in kmer_record and in byte_pac.
  // this is actually a 'lt' not a 'cmp'
  static Bool cmp_basevector( const kmer_record& k1, const kmer_record& k2 )
  {
    // Look at the integers representing the basevector, one by one.
    // Use the bv_table to get the proper alphanumeric ordering.
    for ( int k = 0; k < (K+15)/16; k++ )
      if ( k1.data_.ints[k] != k2.data_.ints[k] )
	return bv_table_convert( k1.data_.ints[k] ) < bv_table_convert( k2.data_.ints[k] );
    return false;
  }

  // less-than for basevector, id, and pos
  static Bool lt_basevector_id_pos( const kmer_record& k1, const kmer_record& k2 )
  {
    if (cmp_basevector(k1, k2)) return true;
    if (cmp_basevector(k2, k1)) return false;
    return id_cmp_pos(k1, k2);
  }

  void ToString(ostream & out) const {
    out << GetId() << " " << GetPos() << " ";
    bvec(K,data_.begin()).Print(out);
  }

  String ToString() const {
      return bvec(K,data_.begin()).ToString();
  }

  void FromString( istream & in ) {
    int id;
    int pos;
    String str_bases;
    in >> id >> pos >> str_bases;
    if ( !in ) return;
    bvec bases( str_bases );
    this->Set( bases, id, pos );
  }

  friend ostream & operator<<(ostream & out, const kmer_record& krec) {
    krec.ToString(out);
    return out;
  }

  friend istream & operator>> ( istream &in, kmer_record &krec ) {
    krec.FromString( in );
    return in;
  }

};  // class kmer_record

//For sorting by id and pos, keeping negative positions negative.
template<class KmerRecord>
struct LessByIdAndPos:
  public binary_function<KmerRecord, KmerRecord, bool> {
  bool operator()(const KmerRecord & k1, const KmerRecord & k2) const {
    read_id_t id1 = k1.GetId( ), id2 = k2.GetId( );
    if ( id1 < id2 ) return true;
    if ( id1 > id2 ) return false;
    return k1.GetPos( )  < k2.GetPos( ) ;
  }
};

//For unique sorting by id and pos to be used in conjunction with less-than operator above.
template<class KmerRecord>
struct EqualByIdAndPos:
  public binary_function<KmerRecord, KmerRecord, bool> {
  bool operator()(const KmerRecord & k1, const KmerRecord & k2) const {
    return  ( k1.GetId() == k2.GetId() ) && ( k1.GetPos() ==  k2.GetPos( ));
  }
};

///Return true if kmers are equal, even if GetPos() and GetId() are different.
template<class KmerRecord>
struct CompareForwardKmers:
  public binary_function<KmerRecord, KmerRecord, bool> {
  bool operator()(const KmerRecord & k1, const KmerRecord & k2) const {
    return k1.EqualKmers(k2);
  }
};

///Return true if kmers are equal in forward or rc directions.
template<class KmerRecord>
struct CompareForwardReverseKmers:
  public binary_function<KmerRecord, KmerRecord, bool> {
  bool operator()(const KmerRecord & k1, const KmerRecord & k2) const {
    return k1.EqualOrReverseKmers(k2);
  }
};



/**
   Class: kmer

   A kmer holds just a kmer.  It allocates just enough space to hold
   exactly K bases.

   Note that there is also a logical type <kmer_t>, which is just a type
   for <basevectors> that happen to be used to hold a single kmer.

   Models <SortKmersOutputRecord>.  See also <kmer_record>, which stores,
   in addition to the kmer sequence, the origin (<read id>) and position
   of the kmer occurrence.
*/
template<int K> class kmer {

  typedef byte_pac< (K+3)/4, 0 > data_t;

     public:
  /// Size of the bases contained in this kmer_record.
  /// For use when creating a basevector from our data.
  static const int BASES_SIZE=K;

     kmer( ) { }

     kmer( const basevector& b )
     {
       Set( b );
     }

     const unsigned char* Bytes( ) const
     {    return data_.bytes;    }

     const unsigned int* Ints( ) const
     {    return data_.ints;    }

     void GetBasevector( basevector& kmer ) const {
         kmer.assignBaseBits(K,data_.begin());
     }

     friend bool operator<( const kmer& k1, const kmer& k2 )
     {    return k1.data_ < k2.data_;    }

     // In the following Set function, read_id and read_pos
     // arguments are not used.  They must still be declared
     // so that class kmer is a valid model of the SortKmersOutputRecord
     // concept.
     void Set( const basevector& b, int read_id, int read_pos )
     {
       SetKmer( data_, b, K );
     }

  // Extract the K-mer from source that begins at start.

  void SetToSubOf( const basevector& source, const size_type start );

  void Set( const basevector& b ) { Set( b, -1, -1 ); }

  void ReverseComplement( ) {
    Set( bvec(K,data_.begin()).ReverseComplement(), -1, -1 );
  }

  friend bool operator==(const kmer & l, const kmer& r) {
    return l.data_ == r.data_;
  }
  friend bool operator!=(const kmer & l, const kmer& r) {
    return !( l == r );
  }

  //   friend bool operator<( const kmer& k1, const kmer& k2 )
  // {    return k1.data_ < k2.data_;    }

  friend bool operator>( const kmer& k1, const kmer& k2 )
  {    return k2 < k1;    }
  friend bool operator>=( const kmer& k1, const kmer& k2 )
  {    return !( k1 < k2 ); }
  friend bool operator<=( const kmer& k1, const kmer& k2 )
  {    return !( k1 > k2 );   }


  String ToString() const {
      return bvec(K,data_.begin()).ToString();
  }

     private:

     data_t data_;

};  // class kmer

template <int K>
struct Serializability< kmer<K> > : public TriviallySerializable {};

template < nbases_t K > inline
ostream& operator<< ( ostream& s, const kmer< K >& k ) {
  basevector b;
  k.GetBasevector( b );
  b.Print( s );
  return s;
}


inline
void CanonicalizeKmer( basevector& b ) {
  b.Canonicalize();
}


// Class: kmer_with_count
//
// A kmer_with_count holds a kmer and a multiplicity.  It consists of:
// * k bases, stored 4 per byte (with the last byte zero-filled if need be);
// * a 2-byte count.
//
// The "count" associated with a kmer is sometimes used to represent not counts
// but other values we may want to associate with the kmer -- for example,
// a 0/1 value indicating whether the kmer is <strong>.
// <KmerShortMap> uses vectors of this class to represent
// vectors of (kmer, short int) pairs.
class kmer_with_count_base
{
public:
    static int const max_count = 65535;
};
template<int K> class kmer_with_count : public kmer_with_count_base {
     private:
  //           (kmer )  (  optional pad   )   (count)
     typedef byte_pac< (K+3)/4, (K+7)/8*2 - (K+3)/4 +    2   > data_t;
  data_t data_;


     public:

     kmer_with_count( ) { }

     kmer_with_count( const basevector& b, unsigned short count ) {
       Set(b, count);
     }

     void GetBasevector( basevector& kmer ) const {
         kmer.assignBaseBits(K,data_.begin());
     }

     void Set(const basevector &b, unsigned short count) {
       SetKmer( data_, b, K );
       for ( int j = (K+3)/4; j < (K+7)/8*2; j++ )
	 data_.bytes[j] = 0;
       data_.shorts[ (K+7)/8 ] = count;
     }

     unsigned short Count( ) const { return data_.shorts[ (K+7)/8 ]; }

     const unsigned int* Ints( ) const
     {    return data_.ints;    }

     const unsigned short* Shorts( ) const
     {    return data_.shorts;    }

     friend Bool operator<( const kmer_with_count& k1, const kmer_with_count& k2 )
     {    if ( k1.data_ < k2.data_ ) return True;
          if ( k2.data_ < k1.data_ ) return False;
          if ( k1.Count( ) < k2.Count( ) ) return True;
          return False;    }

     friend Bool operator==( const kmer_with_count& k1, const kmer_with_count& k2 )
     {    return ( k1.data_ == k2.data_ );    }

     friend Bool eq_kmer( const kmer_with_count& k1, const kmer_with_count& k2 )
     {    return ( 0 == memcmp( k1.data_.bytes, k2.data_.bytes, (K+3)/4 ) );    }

     friend Bool lt_kmer( const kmer_with_count& k1, const kmer_with_count& k2 )
     {    return ( k1.data_ < k2.data_ );     }

  String ToString() const {
      return bvec(K,data_.begin()).ToString();
  }

};  // class kmer_with_count

template <int K>
struct Serializability< kmer_with_count<K> > : public TriviallySerializable {};

template < nbases_t K > inline
ostream& operator<< ( ostream& s, const kmer_with_count< K >& k ) {
  basevector b;
  k.GetBasevector( b );
  b.Print( s );
  return s;
}

#define KmerRecordType(KSHAPE,I) kmer_record_ ## KSHAPE ## _ ## I
#define CreateKmerRecordType(KSHAPE,dummy) \
    typedef kmer_record<KSHAPE::KSIZE,1> KmerRecordType(KSHAPE,1) ;  \
    typedef kmer_record<KSHAPE::KSIZE,2> KmerRecordType(KSHAPE,2)

FOR_ALL_KSHAPES(CreateKmerRecordType,unused);

#define INSTANTIATE_KMER_RECORD_FOR_K(K, dummy) \
  BINARY2_DEF( kmer_with_count<K> ); \
  BINARY3_DEF( kmer_with_count<K> )


#endif
