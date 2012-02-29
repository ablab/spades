// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

#include "paths/ReadFillDatabase.h"

#include <functional>

void ReadFillDatabase::FillsFromFile( const String& filename ) {
  delete mp_RFD;
  
  // Instantiate the appropriate pointed-to RFD implementation
  if( IsRegularFile(filename) )
    mp_RFD = new vectorRFD(filename);
  else
    mp_RFD = new identityRFD();
  
  // Load the left and right trims -- either per-fill trims,
  // or the per-read trims from eg reads.paths.left_trim.k48
  
  String left_trim_file(filename), right_trim_file(filename);
  if( filename.Contains( "_fillrecords." ) ) {
    left_trim_file.ReplaceBy( "_fillrecords.", "_fillrecords.left_trim." );
    right_trim_file.ReplaceBy( "_fillrecords.", "_fillrecords.right_trim." );
  }
  
  if( IsRegularFile(left_trim_file) && IsRegularFile(right_trim_file) ) {
    trim_indexed_by_fill = true;
  }
  else {
    trim_indexed_by_fill = false;
    
    String dirname = Dirname(filename);
    String k48 = FilenameExtension(filename);
    
    left_trim_file = dirname + "/reads.paths.left_trim." + k48;
    right_trim_file = dirname + "/reads.paths.right_trim." + k48;
  }
  
  BinaryRead2( left_trim_file, left_trim );
  BinaryRead2( right_trim_file, right_trim );
}


// The < to binary-search by FirstFilling:
struct ReadFillLess : public binary_function<ReadFillRecord,ReadFillRecord,bool> {
  bool operator() ( const ReadFillRecord& lhs, const ReadFillRecord& rhs ) {
    return( lhs.FirstFilling() < rhs.FirstFilling() );
  }
};

int ReadFillDatabase::vectorRFD::FillToRead( int fill ) const {
  if( FirstFilling(cached_read) <= fill && fill <= LastFilling(cached_read) )
    return cached_read;

  ReadFillRecord dummy( fill, 0,0,0 );
  cached_read = (upper_bound(mp_fills->begin(), mp_fills->end(), dummy, ReadFillLess())
		 - mp_fills->begin()) - 1;
  return cached_read;
}
