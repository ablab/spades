///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef READ_LOCATION
#define READ_LOCATION

#include <fstream>
#include <functional>

#include "CoreTools.h"
#include "CommonSemanticTypes.h"
#include "math/Functions.h"
#include "math/HoInterval.h"

// enum orientation {ForwardOr, ReverseOr};  --  too inefficient -- it uses 4 bytes
typedef char orientation;
#define ForwardOr 0
#define ReverseOr 1

// Class: read_location
//
// Specifies the locations of some sequences (called here 
// "reads") relative to some other sequences (called here "contigs"), which were 
// obtained by some sort of merger process from the reads.  Can also represent
// the location of <simulated> reads on a reference.
//
// If the on_supers_ flag is set, then the positions of the contigs relative to
// ambient supercontigs is also tracked.
//
// Notes (somewhat outdated):
//
// 1. A read can have more than one read_location associated to it, if it crosses
//    more than one contig.
//
// 2. A read need not lie entirely within a contig, because we may not be confident 
//    about the tail of the read.
//
// 3. At present, by "read", we mean a trimmed read.
//
// 4. The bases which a read matches on its contig are always accessible via
//         for ( int i = StartOnContig( ); i <= StopOnContig( ); i++ ).
//    Note that in this sense a read always goes forward on its contig, whether its 
//    orientation is forward or reverse relative to the contig.
//
// The current implementation is not efficient, because LengthOfRead, 
// LengthOfContig, and LengthOfSuperContig are stored for each read_location object.

// NOTE: The read_location class cannot handle datasets with more than 2^31
// reads.  It is not designed for use with the RunAllPathsLG pipeline.
class read_location {

     public:

     read_location( ) { memset(this,0,sizeof(*this)); }

     read_location( int read_id, unsigned int length_of_read, int contig, 
          int start_on_contig, orientation orientation_on_contig, 
          int length_of_contig, short super_contig,
          int start_on_supercontig, int length_of_supercontig ) :
          read_id_(read_id), length_of_read_(length_of_read), contig_(contig),
          start_on_contig_(start_on_contig), 
          orientation_on_contig_(orientation_on_contig),
          length_of_contig_(length_of_contig), on_supers_(True), 
          super_contig_(super_contig), start_on_supercontig_(start_on_supercontig),
          length_of_supercontig_(length_of_supercontig) { }

     read_location( int read_id, unsigned int length_of_read, int contig, 
          int start_on_contig, orientation orientation_on_contig, 
          int length_of_contig ) :
          read_id_(read_id), length_of_read_(length_of_read), contig_(contig),
          start_on_contig_(start_on_contig), 
          orientation_on_contig_(orientation_on_contig),
          length_of_contig_(length_of_contig), on_supers_(False), 
          super_contig_(contig), start_on_supercontig_(start_on_contig),
          length_of_supercontig_(length_of_contig) { }

     read_id_t ReadId( ) const { return read_id_; }
     unsigned int LengthOfRead( ) const { return length_of_read_; }

     int Contig( ) const { return contig_; }

     int StartOnContig( ) const { return start_on_contig_; }
     int StopOnContig( ) const { return start_on_contig_ + length_of_read_ - 1; }
     int Start( ) const { return start_on_contig_; }
     int RcStart( ) const { return LengthOfContig( ) - Start( ) - LengthOfRead( ); }
     int Stop( ) const { return start_on_contig_ + length_of_read_ - 1; }
     ho_interval Interval( ) const 
     {    return ho_interval( StartOnContig( ), StopOnContig( ) );    }
     ho_interval SafeInterval( ) const
     {    return ho_interval( Max( 0, StartOnContig( ) ), 
               Min( LengthOfContig( ), StopOnContig( ) ) );    }

     orientation OrientationOnContig( ) const { return orientation_on_contig_; }
     Bool Fw( ) const { return orientation_on_contig_ == ForwardOr; }
     Bool Rc( ) const { return orientation_on_contig_ == ReverseOr; }
     orient_t Orient() const { return Fw() ? ORIENT_FW : ORIENT_RC; }
     int LengthOfContig( ) const { return length_of_contig_; }

     Bool OnSupers( ) const { return on_supers_; }
     short SuperContig( ) const { return super_contig_; }
     int StartOnSuperContig( ) const { return start_on_supercontig_; }
     int StopOnSuperContig( ) const 
     { return start_on_supercontig_ + length_of_read_ - 1; }
     int LengthOfSuperContig( ) const { return length_of_supercontig_; }

     void SetReadId( int r ) { read_id_ = r; }
     void SetLengthOfRead( int l ) { length_of_read_ = l; }
     void SetContig( int c ) { contig_ = c; }
     void SetLengthOfContig( int l ) { length_of_contig_ = l; }
     void SetLengthOfSuperContig( int l ) { length_of_supercontig_ = l; }
     void SetOrientationOnContig( orientation o ) { orientation_on_contig_ = o; }
     void SetStartOnContig( int s ) { start_on_contig_ = s; }
     void AddToStart( int a ) { start_on_contig_ += a; }
     void SetStartOnSuperContig( int s ) { start_on_supercontig_ = s; }
     void SetOnSupers( Bool v ) { on_supers_ = v; }

     void Kill() { contig_ = -1; }
     bool IsDead() const { return ( contig_ == -1 ); }

     void Reverse( )
     {    SetStartOnContig( LengthOfContig( ) - StopOnContig( ) );
          SetOrientationOnContig(
               OrientationOnContig( ) == ForwardOr ? ReverseOr : ForwardOr );    }

     Bool InBounds( ) const
     {    return StartOnContig( ) >= 0
               && StopOnContig( ) <= LengthOfContig( );    }

     void ForceInBounds( ostream *out_ptr = &cout );

     friend bool operator%(const read_location& oldRead, 
          const read_location& newRead);

     friend bool operator==(const read_location& oldRead, 
          const read_location& newRead);

     friend bool operator!=(const read_location& oldRead,
		     const read_location& newRead) {
       return !operator==(oldRead, newRead);    }

     friend Bool operator<( const read_location& r1, const read_location& r2 )
     {    return r1.Contig( ) < r2.Contig( ) ||
          (r1.Contig( ) == r2.Contig( ) 
               && r1.StartOnContig( ) < r2.StartOnContig( ));    }

     friend ostream& operator<<(ostream &o, const read_location &r);

     private:

     // TODO: potentially dangerous truncation of index by read_id and contig_
     // and other members
     read_id_t read_id_;
     nbases_t length_of_read_;

     int contig_;
     int start_on_contig_;
     orientation orientation_on_contig_;
     int length_of_contig_;

     Bool on_supers_;
     short super_contig_;
     int start_on_supercontig_;
     int length_of_supercontig_;

};  // class read_location

istream& operator>>( istream& s, vec<read_location>& v );

// Function: WriteLocs
//
// Given read locations v, generate three files, locs_file 
// (contains the read locations), locs_file + "_index" (contains the ordered
// read location index at the start of each contig), and locs_file + "_indexr".
// These indices permit random access to locs_file via ReadLocs and ReadsToLocations
// (below).
//
// Upon entry, num_contigs should be the number of contigs, and nreads should be the 
// number of reads.  If either is negative, the corresponding index is not created.
void WriteLocs( const String& locs_file , const vec<read_location>& v, 
     const int num_contigs = -1, int nreads = -1, Bool write_locs = True );

// Function: ReadLocs
//
// Suppose <WriteLocs()> has been run on locs_file, and that contig_ids 
// is a list of contigs, not necessarily sorted.  Then ReadLocs returns a vector v 
// consisting of all the read locations in the given list of contigs.  The order
// within v is first by contig (as given in contig_ids), and then by the order of
// read locations within the contig.
void ReadLocs( const String& locs_file , const vec<int>& contig_ids, 
     vec<read_location>& v );

// Function: NumLocs
//
// Suppose <WriteLocs()> has been run on locs_file.  NumLocs
// resizes num_locs_per_contig to the number of contigs and fills it
// with a map of contig id to number of read locations in that contig.
void NumLocs( const String& locs_file, vec<int> &num_locs_per_contig );

// Function: ReadsToLocations
//
// Suppose WriteLocs has been run on locs_file, and ids is a list
// of read ids.  For each id, return all the read locations for that read.
void ReadsToLocations( const String& locs_file, const vec<int>& ids, 
     vec< vec<read_location> >& locs );

// This functor is useful for searching algorithms like equal_range.

struct order_read_locations_by_contig 
  : public binary_function<read_location,read_location,bool>
{
  bool operator() ( const read_location& lhs, const read_location& rhs) const
    { return ( lhs.Contig() < rhs.Contig() ); }
};

struct order_read_locations_by_readid
  : public binary_function<read_location,read_location,bool>
{
  bool operator() ( const read_location& lhs, const read_location& rhs) const
    { return ( lhs.ReadId() < rhs.ReadId() ); }
};

void AnnotateReadLocations( vec<read_location>& data, String run_dir, Bool orig,
     String human_out, const vec<Bool>& contigs_to_use, Bool gzip = False, 
     Bool untrimmed = False, Bool show_divider = False );

void AnnotateReadLocations2( const String& run_dir, const String& sub_dir,
     const vec<int>& tig_ids, ostream& out );

// Class: read_location_short
//
// Compactly represents the location and orientation of a read on a contig
// (for example on a unipath).
// A read_location_short stores only a read_id, contig_id, start_on_contig, and
// orientation.  It occupies 12 bytes, whereas a read_location occupies 36 bytes.

// NOTE: The read_location_short class is deprecated as of September 2009.
// Instead of using this class, use ReadLocationLG, which can handle larger
// dataset sizes.
// This class continues to exist only in old code from the RunAllPaths pipeline
// that is not part of RunAllPathsLG.

class read_location_short {

     public:

     read_location_short( ) { }

     read_location_short( const read_location& rl )
     {    read_id_ = rl.ReadId( );
          if ( rl.Fw( ) ) contig_rc_ = rl.Contig( );
          else contig_rc_ = -rl.Contig( ) - 1;
          start_ = rl.Start( );    }

     read_location_short( read_id_t read_id, int contig, 
          int start_on_contig, orientation orientation_on_contig )
     {    read_id_ = read_id;
          if ( orientation_on_contig == ForwardOr ) contig_rc_ = contig;
          else contig_rc_ = -contig - 1;
          start_ = start_on_contig;    }

     void SetReadId( read_id_t new_id ) { read_id_ = new_id; }

     read_id_t ReadId( ) const { return read_id_; }
     Bool Rc( ) const { return contig_rc_ < 0; }
     Bool Fw( ) const { return contig_rc_ >= 0; }
     orient_t Orient() const { return Fw() ? ORIENT_FW : ORIENT_RC; }
     int StartOnContig( ) const { return start_; }
     int Start( ) const { return start_; }
     int Contig( ) const
     {    if ( Fw( ) ) return contig_rc_;
          else return -contig_rc_ - 1;    }

     friend Bool operator<( const read_location_short& r1, 
          const read_location_short& r2 )
     {    return r1.Contig( ) < r2.Contig( ) ||
          (r1.Contig( ) == r2.Contig( ) && r1.Start( ) < r2.Start( ) );    }

     private:

     read_id_t read_id_;
     int contig_rc_;
     int start_;

};

inline Bool cmp_contig_read( const read_location_short& r1, 
     const read_location_short& r2 )
{    return r1.Contig( ) < r2.Contig( ) ||
     ( r1.Contig( ) == r2.Contig( ) && r1.ReadId( ) < r2.ReadId( ) );    }

#endif
