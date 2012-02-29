///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


#include "FastaFilestreamPreview.h"

#include <algorithm>

FastaFilestreamPreview::FastaFilestreamPreview(istream& filestream)
  : max_sequence_size_(0), start_offset_(0)
{
  const int bufsize = 81920;
  char buffer[bufsize];
  char *bufptr;

  char breakPoints[256];
  for ( int c = 0; c < 256; ++c )
    breakPoints[ c ] = 0;
  breakPoints[  0  ] = 1;
  breakPoints[ static_cast<unsigned char>( '>' ) ] = 1;

  // read bufsize-1 bytes and put a sentinel after the last byte read
  filestream.read( buffer, bufsize-1 );
  bufptr = buffer + filestream.gcount();
  *bufptr = 0;

  // if nothing was read, return
  if ( bufptr == buffer )
    return;

  bool end_of_file = false;
      
  // what is the offset in the file of the data now in the buffer
  streamoff offset_of_buffer = 0;

  // find the first '>' preceded by nothing at all or a newline
  bufptr = buffer;
  if ( *bufptr == '>' )
    start_offset_ = 0;
  else
  {
    do 
    {
      while ( breakPoints[ static_cast<unsigned char>(*bufptr) ] == 0 )
        ++bufptr;
      if ( *bufptr == '>' && 
           bufptr > buffer &&
           *(bufptr-1) == '\n' )
      {
        start_offset_ = offset_of_buffer + static_cast<int>( bufptr - buffer );
        ++bufptr;
        break;
      }
      if ( *bufptr == 0 )
        // either we hit the end of the buffer or the end of the file
      {
        // try reading in a new buffer
        filestream.unget();
        filestream.read( buffer, bufsize-1 );
        offset_of_buffer += ( bufptr - buffer ) - 1;
        // put a sentinel at the end
        bufptr = buffer + filestream.gcount();
        *bufptr = 0;
        // if data read in is one byte or less, that was the eof
        if ( bufptr - buffer <= 1 )
          end_of_file = true;
        bufptr = buffer;
      }
    }
    while ( ! end_of_file );
  }

  // if we hit the end of the file before we found the first '>', there's no data
  if ( end_of_file )
    return;

  max_sequence_size_ = start_offset_;
  
  streampos size_of_curr_sequence = 0;
  streampos offset_of_last_sequence = start_offset_;

  do
  {
    while ( breakPoints[ static_cast<unsigned char>( *bufptr ) ] == 0 )
      ++bufptr;
    if ( *bufptr == '>' )
    {
      if ( bufptr > buffer &&
           *(bufptr-1) == '\n' )
      {
        streampos offset_of_curr_sequence = offset_of_buffer + static_cast<int>( bufptr - buffer );
        size_of_curr_sequence = offset_of_curr_sequence - offset_of_last_sequence;
        sequence_sizes_.push_back( size_of_curr_sequence );
        max_sequence_size_ = max( max_sequence_size_, size_of_curr_sequence );
        offset_of_last_sequence = offset_of_curr_sequence;
      }
      ++bufptr;
    }
    if ( *bufptr == 0 )
      // either we hit the end of the buffer or the end of the file
    {
      // try reading in a new buffer, including the last character
      filestream.unget();
      filestream.read( buffer, bufsize-1 );
      offset_of_buffer += ( bufptr - buffer ) - 1;
      // insert the sentinel
      bufptr = buffer + filestream.gcount();
      *bufptr = 0;
      // if the data read in is one byte or less, that was the eof
      if ( bufptr - buffer <= 1 )
      {
        end_of_file = true;
        break;
      }
      bufptr = buffer;
    }
  }
  while ( ! end_of_file );

  // we need to make sure it reads to the end of the file,
  // so we insert the offset to the position of the eof
  size_of_curr_sequence = 
    offset_of_buffer
    + static_cast<int>( bufptr - buffer )
    - offset_of_last_sequence + 1;
  max_sequence_size_ = max( max_sequence_size_, size_of_curr_sequence );
  sequence_sizes_.push_back(size_of_curr_sequence);
}

const streampos FastaFilestreamPreview::getMaxSequenceSize()
{
  return max_sequence_size_;
}

vec<streampos>& FastaFilestreamPreview::getSequenceSizes()
{
  return sequence_sizes_;
}

const streampos FastaFilestreamPreview::getStartOffset()
{
  return start_offset_;
}

