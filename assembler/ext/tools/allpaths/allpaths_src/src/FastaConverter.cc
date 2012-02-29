///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


#include "FastaConverter.h"
#include "system/System.h"

#include <strstream>

template<typename sequenceT>
void FastaConverter<sequenceT>::extractNameFromBuffer( char* buffer, String &name )
{ 
  name_parser_->extractNameFromBuffer( buffer, name );
}

template<typename sequenceT>
void FastaConverter<sequenceT>::extractAllFromBuffer( char* buffer, 
                                                      String &name, sequenceT &sequence )
{
  //  get first line of buffer
  char* end_of_first_line = buffer;
  while ( *end_of_first_line != '\n' && 
          *end_of_first_line != 0 )
    ++end_of_first_line;
  
  char* beginning_of_second_line = end_of_first_line+1;
  // ensure we haven't overrun the buffer.
  if ( *end_of_first_line == 0 )
    --beginning_of_second_line;
  
  // terminate the name portion of the buffer
  *end_of_first_line = 0;
  // pass this shortened buffer to extractNameFromBuffer
  extractNameFromBuffer( buffer, name ); 
  
  // pass the rest of the buffer to extractDatumFromBuffer
  if ( ! extractDatumFromBuffer( beginning_of_second_line, sequence ) )
  {
    InputErr( "There was an error processing the quality score for the sequence" << endl
              << name << ".  The file may not be in FASTA format." );
  }
}

template<>
bool
FastaConverter<charvector>::extractDatumFromBuffer( char* buffer, charvector &sequence )
{
  sequence.resize(0);
  while ( *buffer != 0 )
  {
    sequence.push_back( *buffer );
    ++buffer;
  }
  return true;
}


template<>
bool
FastaConverter<CompressedSequence>::extractDatumFromBuffer( char* buffer,
                                                            CompressedSequence &sequence )
{
  static bool initialized = false;
  static int is_valid[256];
  if ( ! initialized )
  {
    for ( int c = 0; c < 256; ++c )
      is_valid[c] = ( isspace( c ) ? 0 : 1 );
    initialized = true;
  }
    
  char *valid_buffer = buffer;
  char *valid_buffer_walk = buffer;

  while ( *buffer != 0 )
  {
    *valid_buffer_walk = *buffer;
    valid_buffer_walk += is_valid[ static_cast<unsigned char>( *buffer++ ) ];
  }
  *valid_buffer_walk = 0;

  sequence = CompressedSequence( valid_buffer );
  
  return true;
}


template<>
bool
FastaConverter<qualvector>::extractDatumFromBuffer( char* buffer,
                                                    qualvector &sequence )
{
  static bool initialized = false;
  static int is_number[256];
  static int is_not_number_nor_null[256];
  static int is_space_or_null[256];
  if ( ! initialized )
  {
    for ( int c = 0; c < 256; ++c )
    {
      is_number[c] = ( isdigit( c ) != 0 );
      is_not_number_nor_null[c] = ! ( c == 0 || isdigit( c ) );
      is_space_or_null[c] = ( c == 0 || isspace( c ) != 0 );
    }
    initialized = true;
  }

  sequence.clear();

  // skip leading whitespace
  while ( is_not_number_nor_null[ static_cast<unsigned char>( *buffer ) ] )
    ++buffer;
  
  if ( *buffer == 0 )
    return true;

  // go through the whole buffer
  int q = 0;
  while ( 1 )
  {
    // if the current char is a digit
    if ( is_number[ static_cast<unsigned char>( *buffer ) ] )
    {
      // then we're extending a number so we
      // add it to the current quality * 10
      q = 10*q + *buffer - '0';

      ++buffer;
    }
    // if it's not a number, we've just read in
    // a quality score, so we save it and reset
    else 
    {
      // if q is too big or the delimiter is not a whitespace char,
      // we inform the user that their file seems to be screwy
      if ( q > 255 || ! is_space_or_null[ static_cast<unsigned char>( *buffer ) ] )
      {
        PRINT2( q, *buffer );
        return false;
      }
      sequence.push_back( static_cast<unsigned char>( q ) );
      q = 0;
      
      // before we continue, we scan up to the next digit
      while ( is_not_number_nor_null[ static_cast<unsigned char>( *buffer ) ] )
        ++buffer;

      if ( *buffer == 0 )
        break;
    }
  }

  return true;
}

#define INSTANTIATE(sequenceT) \
template void FastaConverter<sequenceT>::extractNameFromBuffer( char* buffer, String &name ); \
template void FastaConverter<sequenceT>::extractAllFromBuffer( char* buffer, String &name, \
                                                               sequenceT &sequence );

INSTANTIATE(charvector)
INSTANTIATE(CompressedSequence)
INSTANTIATE(qualvector)
