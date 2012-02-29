///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef READ_LOCATION_LG
#define READ_LOCATION_LG

#include "CommonSemanticTypes.h" // ORIENT_FW, ORIENT_RC
#include "VecTemplate.h" // BINARY2_DEF

/* ReadLocationLG
 *
 * A ReadLocationLG represents the oriented location of a read on a contig.
 * It typically appears as a member of a vec<ReadLocationLG>, which is a set
 * of all such locations among a set of reads and contigs.  Here the 'reads' are
 * query sequences, while the 'contigs' are target sequences; typically the
 * contigs are either a set of unipaths or a reference genome.
 * 
 * Note that any read and any contig may appear in more than one ReadLocationLG,
 * or in no ReadLocationLG's.
 *
 * The ReadLocationLG class is designed to replace the older read_location and
 * read_location_short classes (ReadLocation.h) throughout the RunAllPathsLG
 * pipeline.  A ReadLocationLG uses much less memory than a read_location.
 * It can handle larger datasets with more than 2^31 reads, unlike read_location
 * and read_location_short.
 *
 *
 * MEMORY USAGE:
 * 16 bytes.
 *
 * NUMERIC LIMITATIONS:
 * Number of reads in the dataset must be less than max(longlong) = 2^63 =~ 1e19
 * Number of contigs in the dataset must be less than max(int) = 2^31
 * All contig lengths must be less than max(int) = 2^31
 *
 *
 * Josh Burton
 * September 2009
 *
 ******************************************************************************/
class ReadLocationLG {
  
public:
  
  // -------- CONSTRUCTORS
  ReadLocationLG( ) { }
  
  ReadLocationLG( longlong read_id, int contig, 
		  int start_on_contig, Bool orientation_on_contig ) {
    read_id_ = read_id;
    if ( orientation_on_contig == ORIENT_FW ) contig_rc_ = contig;
    else contig_rc_ = -contig - 1;
    start_ = start_on_contig;
  }
  
  
  // -------- MODIFIERS
  void SetReadId( longlong new_id ) { read_id_ = new_id; }
  
  
  // -------- QUERY FUNCTIONS
  longlong ReadId( ) const { return read_id_; }
  Bool Rc( ) const { return contig_rc_ < 0; }
  Bool Fw( ) const { return contig_rc_ >= 0; }
  Bool Orient() const { return Fw() ? ORIENT_FW : ORIENT_RC; }
  int StartOnContig( ) const { return start_; }
  int Start( ) const { return start_; }
  int Contig( ) const {
    if ( Fw( ) ) return contig_rc_;
    else return -contig_rc_ - 1;
  }
  
  friend Bool operator<( const ReadLocationLG& r1, 
			 const ReadLocationLG& r2 ) {
    return r1.Contig( ) < r2.Contig( ) ||
      (r1.Contig( ) == r2.Contig( ) && r1.Start( ) < r2.Start( ) );
  }
  
  friend bool operator==( const ReadLocationLG & r1, const ReadLocationLG & r2 ) {
    if ( r1.Contig() != r2.Contig() ) return false;
    if ( r1.Start()  != r2.Start()  ) return false;
    return true;
  }
  
  
  // -------- MEMBER VARIABLES
private:
  
  longlong read_id_;
  int contig_rc_;
  int start_;
  
};

inline Bool cmp_contig_read( const ReadLocationLG& r1, 
			     const ReadLocationLG& r2 )
{
  return r1.Contig( ) < r2.Contig( ) ||
    ( r1.Contig( ) == r2.Contig( ) && r1.ReadId( ) < r2.ReadId( ) );
}



inline bool contig_read_start_lt(const ReadLocationLG & r1, 
                                 const ReadLocationLG & r2)
{
  if (r1.Contig() < r2.Contig()) return true;
  if (r1.Contig() > r2.Contig()) return false;
  if (r1.ReadId() < r2.ReadId()) return true;
  if (r1.ReadId() > r2.ReadId()) return false;
  if (r1.Start()  < r2.Start())  return true;
  if (r1.Start()  > r2.Start())  return false;
  return false;
}

inline bool contig_read_start_eq(const ReadLocationLG & r1, 
                                 const ReadLocationLG & r2)
{
  if (r1.Contig() != r2.Contig()) return false;
  if (r1.ReadId() != r2.ReadId()) return false;
  if (r1.Start()  != r2.Start())  return false;
  return true;
}

  

// Instantiate templatized I/O functions for use with vec<ReadLocationLG>
// (chiefly BinaryRead2 and BinaryWrite2.)
BINARY2_DEF( ReadLocationLG );


#endif
