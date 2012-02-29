///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef READ_LOC_H
#define READ_LOC_H

#include "CoreTools.h"
#include "PackAlign.h"
#include "Superb.h"
#include "math/Functions.h"
#include "polymorphism/DumbCall.h"
#include "system/file/FileReader.h"

// A read_loc is yet another read location data structure.  It takes 40 bytes in
// memory and has a 29 byte packed representation on disk:
//
//   bytes+bits      type          content
// disk    memory
//
// 4+1     8         unsigned      read id
// 0+0     4         unsigned      contig id (4+0 but no storage required on disk)
// 4+0     4         signed        pos on contig
// 0+1     1         unsigned      orientation
// 2+0     2         unsigned      read length
// 1+0     1         unsigned      bandwidth for alignment
// 4+1     8         unsigned      partner read id
// 4+0     4         unsigned      partner contig id
// 4+0     4         signed        partner pos on contig
// 0+1     1         unsigned      partner orientation
// 2+0     2         unsigned      partner read length
// 1+0     1         unsigned      partner bandwidth for alignment
// 0+4     1         unsigned      read class
// 1+0     1         unsigned      lib id in read class
// 0+1     1         unsigned      partner placed?
//
// where read class is
//
// 0 = frag reads
// 1 = jump reads
// 2 = long jump reads
//
// Contig id: no storage is required because the locs are sorted by contig and
// a separate index is maintained.
//
// Bandwidth for alignment: this field was added after various pieces of code
// were written.  It facilitates the recovery of alignments by running 
// SmithWatBandedA using the negative of "pos on contig" as the offset and the 
// given bandwidth.  We make the convention that the bandwidth should be as small
// as possible, and subject to that constraint, the "pos on contig" should be as
// small as possible.

const int8_t    pow_2_4  = 16;
const int32_t   pow_2_24 = 16777216;
const int64_t   pow_2_33 = 8589934592;

// ==============================================================================
// CLASSES PackedStuff, UnpackedStuff
// ==============================================================================

// Helper classes to pack/unpack core data structures of ReadLoc.

class PackedStuff;

class UnpackedStuff
{
public:
    UnpackedStuff( ) { }

    UnpackedStuff( uint64_t readID, uint64_t partnerReadID, 
         uint32_t partnerContigID, int32_t contigPosn, int32_t partnerContigPosn, 
         bool isRC, bool partnerIsRC, uint8_t readClass, uint8_t libID, 
         uint16_t readLen, uint16_t partnerReadLen, uint8_t bandwidth,
         uint8_t partnerBandwidth, bool partnerPlaced )
    : mReadID(readID), mPartnerReadID(partnerReadID), 
      mPartnerContigID(partnerContigID), mContigPosn(contigPosn),
      mPartnerContigPosn(partnerContigPosn), mIsRC(isRC), mPartnerIsRC(partnerIsRC), 
      mReadClass(readClass), mLibID(libID), mReadLen(readLen), 
      mPartnerReadLen(partnerReadLen), mBandwidth(bandwidth),
      mPartnerBandwidth(partnerBandwidth), mPartnerPlaced(partnerPlaced) {}

    explicit UnpackedStuff( PackedStuff const& );

    // compiler-supplied copying and destructor are OK

    uint64_t getReadID() const { return mReadID; }
    int32_t getContigPosn() const { return mContigPosn; }
    bool isRC() const { return mIsRC; }
    bool partnerIsRC() const { return mPartnerIsRC; }
    uint16_t getReadLen() const { return mReadLen; }
    uint8_t getBandwidth() const { return mBandwidth; }
    uint64_t getPartnerReadID() const { return mPartnerReadID; }
    uint32_t getPartnerContigID() const { return mPartnerContigID; }
    int32_t getPartnerContigPosn() const { return mPartnerContigPosn; }
    bool PartnerisRC() const { return mPartnerIsRC; }
    uint16_t getPartnerReadLen() const { return mPartnerReadLen; }
    uint8_t getPartnerBandwidth() const { return mPartnerBandwidth; }
    uint8_t getReadClass() const { return mReadClass; }
    uint8_t getLibID() const { return mLibID; }
    bool partnerPlaced() const { return mPartnerPlaced; }

private:
    uint64_t mReadID;
    uint64_t mPartnerReadID;
    int64_t  mPartnerContigID;
    int32_t  mContigPosn;
    int32_t  mPartnerContigPosn;
    bool mIsRC;
    bool mPartnerIsRC;
    uint8_t mReadClass;
    uint8_t mLibID;
    uint16_t mReadLen;
    uint16_t mPartnerReadLen;
    uint8_t mBandwidth;
    uint8_t mPartnerBandwidth;
    bool mPartnerPlaced;

};

class PackedStuff
{
public:
    PackedStuff( ) { }

    PackedStuff( uint64_t readID, uint64_t partnerReadID, 
         uint32_t partnerContigID, int32_t contigPosn, int32_t partnerContigPosn, 
         bool isRC, bool partnerIsRC, uint8_t readClass, uint8_t libID, 
         uint16_t readLen, uint16_t partnerReadLen, uint8_t bandwidth,
         uint8_t partnerBandwidth, bool partnerPlaced )
    : mReadID(readID), mPartnerReadID(partnerReadID), 
      mPartnerContigID(partnerContigID), mContigPosn(contigPosn),
      mPartnerContigPosn(partnerContigPosn), mIsRC(isRC), mPartnerIsRC(partnerIsRC), 
      mReadClass(readClass), mLibID(libID), mReadLen(readLen), 
      mPartnerReadLen(partnerReadLen), mBandwidth(bandwidth),
      mPartnerBandwidth(partnerBandwidth), mPartnerPlaced(partnerPlaced)
    { ForceAssertLt(readID,1ul<<33); 
      ForceAssertLt(partnerReadID,1ul<<33); 
      ForceAssertLt(partnerContigID,1ul<<32); 
      ForceAssertLt(contigPosn,1l<<32); 
      ForceAssertLt(partnerContigPosn,1l<<32); 
      ForceAssertLt(readClass,1ul<<4);
      ForceAssertLt(libID,1ul<<8); 
      ForceAssertLt(readLen,1ul<<16);
      ForceAssertLt(partnerReadLen,1ul<<16);
      ForceAssertLt(bandwidth,1ul<<8);
      ForceAssertLt(partnerBandwidth,1ul<<8); }

    explicit PackedStuff( UnpackedStuff const& );

    // compiler-supplied copying and destructor are OK

    uint64_t getReadID() const { return mReadID; }
    int32_t getContigPosn() const { return mContigPosn; }
    bool isRC() const { return mIsRC; }
    bool partnerIsRC() const { return mPartnerIsRC; }
    uint16_t getReadLen() const { return mReadLen; }
    uint8_t getBandwidth() const { return mBandwidth; }
    uint64_t getPartnerReadID() const { return mPartnerReadID; }
    uint32_t getPartnerContigID() const { return mPartnerContigID; }
    int32_t getPartnerContigPosn() const { return mPartnerContigPosn; }
    bool PartnerisRC() const { return mPartnerIsRC; }
    uint16_t getPartnerReadLen() const { return mPartnerReadLen; }
    uint8_t getPartnerBandwidth() const { return mPartnerBandwidth; }
    uint8_t getReadClass() const { return mReadClass; }
    uint8_t getLibID() const { return mLibID; }
    bool partnerPlaced() const { return mPartnerPlaced; }

private:
    uint64_t mReadID : 33;
    uint64_t mPartnerReadID : 33;
    int64_t mPartnerContigID : 32;
    int32_t mContigPosn : 32;
    int32_t mPartnerContigPosn : 32;
    uint8_t mIsRC : 1;
    uint8_t mPartnerIsRC : 1;
    uint8_t mReadClass : 4;
    uint8_t mLibID : 8;
    uint16_t mReadLen : 16;
    uint16_t mPartnerReadLen : 16;
    uint8_t mBandwidth : 8;
    uint8_t mPartnerBandwidth : 8;
    uint8_t mPartnerPlaced : 1;

};

inline UnpackedStuff::UnpackedStuff( PackedStuff const& ps )
: mReadID(ps.getReadID()), mPartnerReadID(ps.getPartnerReadID()), 
      mPartnerContigID(ps.getPartnerContigID()), mContigPosn(ps.getContigPosn()),
      mPartnerContigPosn(ps.getPartnerContigPosn()), mIsRC(ps.isRC()), 
      mPartnerIsRC(ps.partnerIsRC()), mReadClass(ps.getReadClass()), 
      mLibID(ps.getLibID()), mReadLen(ps.getReadLen()), 
      mPartnerReadLen(ps.getPartnerReadLen()), mBandwidth(ps.getBandwidth()), 
      mPartnerBandwidth(ps.getPartnerBandwidth()), 
      mPartnerPlaced(ps.partnerPlaced())
{}

inline PackedStuff::PackedStuff( UnpackedStuff const& ups )
: mReadID(ups.getReadID()), mPartnerReadID(ups.getPartnerReadID()), 
      mPartnerContigID(ups.getPartnerContigID()), mContigPosn(ups.getContigPosn()),
      mPartnerContigPosn(ups.getPartnerContigPosn()), mIsRC(ups.isRC()), 
      mPartnerIsRC(ups.partnerIsRC()), mReadClass(ups.getReadClass()), 
      mLibID(ups.getLibID()), mReadLen(ups.getReadLen()), 
      mPartnerReadLen(ups.getPartnerReadLen()), mBandwidth(ups.getBandwidth()), 
      mPartnerBandwidth(ups.getPartnerBandwidth()), 
      mPartnerPlaced(ups.partnerPlaced())
    { ForceAssertLt(ups.getReadID(),1ul<<33); 
      ForceAssertLt(ups.getPartnerReadID(),1ul<<33); 
      ForceAssertLt(ups.getPartnerContigID(),1ul<<32); 
      ForceAssertLt(ups.getContigPosn(),1l<<32); 
      ForceAssertLt(ups.getPartnerContigPosn(),1l<<32); 
      ForceAssertLt(ups.getReadClass(),1ul<<4);
      ForceAssertLt(ups.getLibID(),1ul<<8); 
      ForceAssertLt(ups.getReadLen(),1ul<<16);
      ForceAssertLt(ups.getPartnerReadLen(),1ul<<16);
      ForceAssertLt(ups.getBandwidth(),1ul<<8);
      ForceAssertLt(ups.getPartnerBandwidth(),1ul<<8); }

// ==============================================================================
// CLASS ReadLoc
// ==============================================================================

class read_loc {

public:

  read_loc( ) { memset(this,0,sizeof(*this)); }
  
  read_loc( const uint64_t read_id, const uint64_t partner_read_id, 
         const uint32_t contig_id,
         const uint32_t partner_contig_id, const int32_t pos_on_contig, 
         const int32_t partner_pos_on_contig, const bool fw_on_contig, 
         const bool partner_fw_on_contig, const uint8_t read_class, 
         const uint8_t library_id, const uint16_t read_length, 
         const uint16_t partner_read_length, const uint8_t bandwidth,
         const uint8_t partner_bandwidth, const bool partner_placed ) :
    read_id_(read_id), partner_read_id_(partner_read_id),
    contig_id_(contig_id), partner_contig_id_(partner_contig_id), 
    pos_on_contig_(pos_on_contig), partner_pos_on_contig_(partner_pos_on_contig), 
    read_length_(read_length), partner_read_length_(partner_read_length),
    library_id_(library_id), bandwidth_(bandwidth), 
    partner_bandwidth_(partner_bandwidth), fw_on_contig_(fw_on_contig),
    partner_fw_on_contig_(partner_fw_on_contig), read_class_(read_class), 
    partner_placed_(partner_placed)
  { AssertLt(read_id,1ul<<33); 
    AssertLt(partner_read_id,1ul<<33); 
    AssertLt(partner_contig_id,1ul<<32); 
    AssertLt(pos_on_contig,1l<<32); 
    AssertLt(partner_pos_on_contig,1l<<32); 
    AssertLt(read_class,1ul<<4);
    AssertLt(library_id,1ul<<8); 
    AssertLt(read_length,1ul<<16);
    AssertLt(partner_read_length,1ul<<16);
    AssertLt(bandwidth,1ul<<8);
    AssertLt(partner_bandwidth,1ul<<8); }

  read_loc( const uint64_t read_id, const uint32_t contig_id,
         const int32_t pos_on_contig, const bool fw_on_contig, 
         const uint8_t read_class, const uint8_t library_id, 
         const uint16_t read_length ) : read_id_(read_id), partner_read_id_(0),
    contig_id_(contig_id), partner_contig_id_(0), pos_on_contig_(pos_on_contig),
    partner_pos_on_contig_(0), read_length_(read_length), partner_read_length_(0),
    library_id_(library_id), bandwidth_(0), partner_bandwidth_(0), 
    fw_on_contig_(fw_on_contig), partner_fw_on_contig_(false), 
    read_class_(read_class), partner_placed_(false)
  { AssertLt(read_id,1ul<<33); 
    AssertLt(pos_on_contig,1l<<32); 
    AssertLt(read_class,1ul<<4);
    AssertLt(library_id,1ul<<8); 
    AssertLt(read_length,1ul<<16); }

  read_loc( const align& a, const uint64_t read_id, const uint32_t contig_id,
         const bool fw_on_contig, const uint8_t read_class, 
         uint8_t library_id, const uint16_t read_length );

  read_loc( const PackedStuff &packed, const uint32_t contig_id )
  {
    contig_id_ = contig_id;
    UnpackedStuff unpacked( packed );

    read_id_ = unpacked.getReadID();
    pos_on_contig_ = unpacked.getContigPosn();
    fw_on_contig_ = ! unpacked.isRC();
    read_length_ = unpacked.getReadLen();
    bandwidth_ = unpacked.getBandwidth();
    partner_read_id_ = unpacked.getPartnerReadID();
    partner_contig_id_ = unpacked.getPartnerContigID();
    partner_pos_on_contig_ = unpacked.getPartnerContigPosn();
    partner_fw_on_contig_ = ! unpacked.PartnerisRC();
    partner_read_length_ = unpacked.getPartnerReadLen();
    partner_bandwidth_ = unpacked.getPartnerBandwidth();
    read_class_ = unpacked.getReadClass();
    library_id_ = unpacked.getLibID();
    partner_placed_ = unpacked.partnerPlaced();
  }

  uint64_t ReadId( ) const { return read_id_; }
  bool PartnerPlaced( ) const { return partner_placed_; }
  void SetPartnerPlaced( ) { partner_placed_ = true; }
  uint64_t PartnerReadId( ) const { return partner_read_id_; }
  void SetReadId( const uint64_t read_id) { read_id_ = read_id; }
  void SetPartnerReadId( const uint64_t partner_read_id ) 
  { partner_read_id_ = partner_read_id; }
  uint32_t ContigId( ) const { return contig_id_; }
  uint32_t PartnerContigId( ) const { return partner_contig_id_; }
  void SetPartnerContigId( const uint32_t partner_contig_id ) 
  { partner_contig_id_ = partner_contig_id; }
  int32_t Start( ) const { return pos_on_contig_; }
  int32_t PartnerStart( ) const { return partner_pos_on_contig_; }
  void SetPartnerStart( const uint32_t partner_pos_on_contig )
  { partner_pos_on_contig_ = partner_pos_on_contig; }
  int32_t Stop( ) const { return pos_on_contig_ + read_length_; }
  int32_t PartnerStop( ) const 
  { return partner_pos_on_contig_ + partner_read_length_; }
  bool Fw( ) const { return fw_on_contig_; }
  bool Rc( ) const { return !fw_on_contig_; }
  bool PartnerFw( ) const { return partner_fw_on_contig_; }
  void SetPartnerFw( const bool partner_fw_on_contig )
  { partner_fw_on_contig_ = partner_fw_on_contig; }
  bool PartnerRc( ) const { return !partner_fw_on_contig_; }
  uint8_t LibraryId( ) const { return library_id_; }
  uint16_t ReadLength( ) const { return read_length_; }
  uint16_t PartnerReadLength( ) const { return partner_read_length_; }
  void SetPartnerReadLength( const uint16_t partner_read_length )
  { partner_read_length_ = partner_read_length; }
  uint8_t Bandwidth( ) const { return bandwidth_; }
  uint8_t PartnerBandwidth( ) const { return partner_bandwidth_; }
  void SetPartnerBandwidth( const uint8_t partner_bandwidth )
  { partner_bandwidth_ = partner_bandwidth; }
  int Sep( ) const;
  int Dev( ) const;

  // Read class stuff.

  uint8_t ReadClass( ) const { return read_class_; }
  String ReadClassName( ) const;
  bool Frag( ) const { return ReadClass( ) == 0; }
  bool Jump( ) const { return ReadClass( ) == 1; }
  bool LongJump( ) const { return ReadClass( ) == 2; }

  friend Bool SamePartnerContigId( const read_loc& rl1, const read_loc& rl2 )
  {    return rl1.PartnerContigId( ) == rl2.PartnerContigId( );    }
  friend Bool SameStart( const read_loc& rl1, const read_loc& rl2 )
  {    return rl1.Start( ) == rl2.Start( );    }
  friend Bool SamePartnerStart( const read_loc& rl1, const read_loc& rl2 )
  {    return rl1.PartnerStart( ) == rl2.PartnerStart( );    }
  friend Bool SameFw( const read_loc& rl1, const read_loc& rl2 )
  {    return rl1.Fw( ) == rl2.Fw( );    }
  friend Bool SamePartnerFw( const read_loc& rl1, const read_loc& rl2 )
  {    return rl1.PartnerFw( ) == rl2.PartnerFw( );    }
  friend Bool SameLibraryId( const read_loc& rl1, const read_loc& rl2 )
  {    return rl1.LibraryId( ) == rl2.LibraryId( );    }

  uint64_t PseudoPairId( ) const { return Min( ReadId( ), PartnerReadId( ) ); }

  // GetAlign: recover the alignment corresponding to a read location.  The 
  // basevector b should be reversed before calling if Rc( ).

  void GetAlign( align& a, const basevector& b, const basevector& t,
     const Bool partner = False ) const;

  void PrintVisualLoc( const Bool visual_abbr, ostream& out, 
       const String& read_name, const basevector& b, const basevector& t,
       const qualvector* q = 0, const Bool partner = False ) const;

  String AsCsvString( ) const {                     // OLD!
    return ToString( read_id_ ) + ","
      + ToString( contig_id_ ) + ","
      + ToString( pos_on_contig_ ) + ","
      + ( fw_on_contig_ ? "+" : "-" ) + ","
      + ToString( (int)read_class_ ) + ","
      + ToString( (uint)library_id_ ) + ","
      + ToString( (uint)read_length_ )
      + ToString( (uint)bandwidth_ );
  }

  void LoadSepDev( const String& run_dir );

  friend void LoadLibNames( const String& run_dir,
			    vec< vec<String> >& lib_names );

  friend ostream& operator<<( ostream& out, const read_loc& rl );
  
  PackedStuff PackWithoutContigId( ) const
  {
    return PackedStuff( read_id_, partner_read_id_,
         partner_contig_id_, pos_on_contig_, partner_pos_on_contig_, 
         !fw_on_contig_, !partner_fw_on_contig_, read_class_, library_id_, 
         read_length_, partner_read_length_, bandwidth_,
         partner_bandwidth_, partner_placed_ );
  }

  friend bool operator==( const read_loc& rl1, const read_loc& rl2 )  // OLD!!
  {    return rl1.ReadId( ) == rl2.ReadId( )
            && rl1.ContigId( ) == rl2.ContigId( ) && rl1.Start( ) == rl2.Start( )
            && rl1.Stop( ) == rl2.Stop( ) && rl1.Fw( ) == rl2.Fw( )
            && rl1.ReadClass( ) == rl2.ReadClass( )
            && rl1.LibraryId( ) == rl2.LibraryId( )
            && rl1.ReadLength( ) == rl2.ReadLength( )
            && rl1.Bandwidth( ) == rl2.Bandwidth( );    }
  friend bool operator!=( const read_loc& rl1, const read_loc& rl2 )
  {    return !( rl1 == rl2 );    }
  
  friend bool operator<( const read_loc& rl1, const read_loc& rl2 )
  {
    if ( rl1.ContigId( ) < rl2.ContigId( ) ) return True;
    if ( rl1.ContigId( ) > rl2.ContigId( ) ) return False;
    if ( rl1.Start( ) < rl2.Start( ) ) return True;
    if ( rl1.Start( ) > rl2.Start( ) ) return False;
    if ( rl1.Fw( ) < rl2.Fw( ) ) return True;
    if ( rl1.Fw( ) > rl2.Fw( ) ) return False;
    return rl1.ReadId( ) < rl2.ReadId( );
  }

  // Catywampus: given two reads from a pair that land on the same contig, return
  // True if they are facing the wrong direction.

  Bool Catywampus( ) const
  {    Assert( PartnerPlaced( ) );
       AssertEq( ContigId( ), PartnerContigId( ) );
       return Fw( ) == PartnerFw( );    }

  // Backwards: given two reads from a pair that land on the same contig, return
  // True if they are facing the wrong direction.

  Bool Backwards( ) const
  {    Assert( PartnerPlaced( ) );
       AssertEq( ContigId( ), PartnerContigId( ) );
       if ( Fw( ) && PartnerRc( ) && Start( ) > PartnerStart( ) ) return True;
       if ( PartnerFw( ) && Rc( ) && PartnerStart( ) > Start( ) ) return True;
       return False;    }

  // Degenerate: give two reads from a pair that land on the same contig, return
  // True if they're a jump and overlap, but are not Catywampus or Backwards.
  // Such pairs presumably arise by circularizing very small fragments.

  Bool Degenerate( ) const
  {    if ( Catywampus( ) || Backwards( ) ) return False;
       if ( ReadClass( ) == 0 ) return False;
       if ( Fw( ) && Stop( ) > PartnerStart( ) ) return True;
       if ( PartnerFw( ) && PartnerStop( ) > Start( ) ) return True;
       return False;    }

  int ActualSep( ) const
  {    Assert( PartnerPlaced( ) );
       AssertEq( ContigId( ), PartnerContigId( ) );
       if ( Fw( ) && PartnerRc( ) ) return PartnerStart( ) - Stop( );
       if ( PartnerFw( ) && Rc( ) ) return Start( ) - PartnerStop( );
       return 0;    }

  // SepDelta: given two reads from a pair that land on the same contig, return
  // the deviance of their separation from the expected value.  Return zero if 
  // reads are in the same orientation.

  int SepDelta( ) const
  {    Assert( PartnerPlaced( ) );
       AssertEq( ContigId( ), PartnerContigId( ) );
       if ( Fw( ) && PartnerRc( ) ) return PartnerStart( ) - Stop( ) - Sep( );
       if ( PartnerFw( ) && Rc( ) ) return Start( ) - PartnerStop( ) - Sep( );
       return 0;    }

  // PrintWithPairInfo takes as input a sorted vec of read_locs for one contig.

  friend void PrintWithPairInfo( ostream& out, const vec<read_loc>& rl,
      const vec<superb>& scaffolds, const vec<int>& to_super, 
      const vec<int>& to_super_pos,
      const vec<int>& to_super_posr, const String& run_dir, const String& sub_dir, 
      const String& assembly, const Bool SHOW_ALIGNS, 
      const Bool SHOW_PARTNER_ALIGNS, const Bool SHOW_QUALS_IN_ALIGNS );

  // PrintSAMAligns takes as input a sorted vec of read_locs.

  friend void PrintSAMAligns( ostream& sam, const vec<read_loc>& rl,
      const vec<int>& to_super, const vec<int>& to_super_pos,
      const vec< vec<int> >& to_super_offset, const String& run_dir,
      const String& sub_dir, const String& assembly );

  // Generate a pileup from the read locations for a single contig.  There are two
  // versions, depending on whether the reads are loaded into memory.

  friend void Pileup( const basevector& tig, const vec<read_loc>& rl,
      const vecbasevector& frag_reads, const vecbasevector& jump_reads,
      const vecbasevector& long_jump_reads, vec<dumbcall>& calls );
  friend void Pileup( const basevector& tig, const vec<read_loc>& rl,
      const String& run_dir, vec<dumbcall>& calls, 
      const Bool get_frags, const Bool get_jumps, const Bool get_longs );
  friend void Pileup( const basevector& tig, const vec<read_loc>& rl,
      const String& run_dir, vec<dumbcall>& calls );

  friend void AddToPileup( const read_loc& rl, const basevector& b, 
     const basevector& tig, vec<dumbcall>& calls );

private:

  uint64_t read_id_;
  uint64_t partner_read_id_;
  uint32_t contig_id_;
  uint32_t partner_contig_id_;
  int32_t pos_on_contig_;
  int32_t partner_pos_on_contig_;
  uint16_t read_length_;
  uint16_t partner_read_length_;
  uint8_t library_id_;
  uint8_t bandwidth_;
  uint8_t partner_bandwidth_ ;
  uint8_t fw_on_contig_ : 1;
  uint8_t partner_fw_on_contig_ : 1;
  uint8_t read_class_ : 4;
  uint8_t partner_placed_ : 1;

  static vec< vec<int> > lib_sep_, lib_dev_;

};

// ==============================================================================
// INTERFACE FOR READING AND WRITING
// ==============================================================================

// The i/o interface consists of a function WriteReadLocs, which generates a file
// from readlocs in memory, and a class read_locs_on_disk, which loads readlocs
// for a given contig.

// WriteReadLocs.  Write a *sorted* vector of read_locs to HEAD.readlocs,
// along with an index, included as part of the file.

void WriteReadLocs( const String& HEAD, const int32_t ntigs, 
     const vec<read_loc>& locs );

class read_locs_on_disk {

     public:

     read_locs_on_disk( const String& HEAD, const String& run_dir );

     void LoadContig( const int tig, vec<read_loc>& locs );
     void WriteContig(const String& file,  const int tig, const vec<read_loc>& locs );
     void LoadOne( const int64_t id, read_loc& rl );
     int64_t getLocsSize() { return locs_size_;}

     private:

     vec<int64_t> index_;
     int64_t locs_size_;
     FileReader fr_;
};

#endif
