///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "FastaNameParser.h"
#include "system/Assert.h"

void
FullNameParser::extractNameFromBuffer( char* buffer, String& name ) const
{
  // If the buffer is empty, set name to empty string
  if ( *buffer == 0 ) {
    name.resize(0);
    return;
  }

  AssertEq( *buffer, '>' );

  name = buffer+1;
}

void
FirstWordParser::extractNameFromBuffer( char* buffer, String& name ) const
{
  // If the buffer is empty, clear the name and return.
  if ( *buffer == 0 )
  {
    name.resize(0);
    return;
  }

  AssertEq( *buffer, '>' );

  // The start of the name is the first non-space character.
  char* name_start_ptr = buffer+1;
  while ( *name_start_ptr != 0 && isspace(*name_start_ptr) )
    ++name_start_ptr;

  if ( *name_start_ptr == 0 )
  {
    name.resize(0);
    return;
  }

  // The end of the name is the last non-space character after the start.
  char* space_finder = name_start_ptr+1;
  while ( *space_finder != 0 && !isspace(*space_finder) )
    ++space_finder;
  *space_finder = 0;
  
  name = name_start_ptr;
}

void
LastWordParser::extractNameFromBuffer( char* buffer, String& name ) const
{
  // If the buffer is empty, set name to empty string.
  if ( *buffer == 0 ) {
    name.resize(0);
    return;
  }

  AssertEq( *buffer, '>' );

  ++buffer;
  char* last_word_ptr = buffer;

  while ( *buffer != 0 )
  {
    // zero-out space characters
    while ( *buffer != 0 && isspace(*buffer) ) 
      *buffer++ = 0;
    // if we're not at the end of the buffer, we're at the start of a word
    if ( *buffer != 0 )
      last_word_ptr = buffer;
    // skip non-space characters
    while ( *buffer != 0 && !isspace(*buffer) )
      ++buffer;
  }

  name = last_word_ptr;
}

void
TruncatedLastWordParser::extractNameFromBuffer( char* buffer, String& name ) const
{
  LastWordParser::extractNameFromBuffer( buffer, name );
  
  int len = name.size();
  if ( name.size() > 4 &&
       ( strcmp( name.c_str() + len - 4, ".scf" ) == 0 ||
         strcmp( name.c_str() + len - 4, ".exp" ) == 0 ) )
    name.resize( len - 4 );
}

void
Riken_cDNA_Parser::extractNameFromBuffer( char* buffer, String& name ) const
{
  if ( *buffer == 0 ) {
    name.resize(0);
    return;
  }

  name = String( buffer ).Before ( " " );
  for (int ii=0; ii<2; ii++)
    name = name.After( "|" );
  name = name.Before( "|" );
}

