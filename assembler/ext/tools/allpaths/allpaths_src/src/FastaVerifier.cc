///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


#include "FastaVerifier.h"

#include <cctype>
#include <string.h>

bool FastaVerifier::verifyLine( const char* line )
{
  if ( line[0] == '>' )
  {
    ignore_rest_of_line_ = true;
    return true;
  }

  ignore_rest_of_line_ = false;

  int length = strlen( line );
  
  for ( int i = 0; i < length; ++i )
    if (! this->verifyChar( line[i] ) )
      return false;
  
  return true;
}

bool FastaVerifier::verifyRestOfLine( const char* line )
{
  if ( ignore_rest_of_line_ )
    return true;

  int length = strlen( line );
  
  for ( int i = 0; i < length; ++i )
    if (! this->verifyChar( line[i] ) )
      return false;
  
  return true;
}
 
bool FastaSequenceVerifier::verifyChar( const char c )
{
  if ( c == 'A' ||
       c == 'C' ||
       c == 'G' ||
       c == 'T' ||
       
       c == 'N' ||
       c == 'X' ||
       
       isspace(c) || 
       
       c == 'a' ||
       c == 'c' ||
       c == 'g' ||
       c == 't' ||
       
       c == 'n' ||
       c == 'x' )
    return true;
  
  return false;
}
  
bool FastaQualityVerifier::verifyChar( const char c )
{
  if ( isdigit(c) || isspace(c) )
    return true;
  
  return false;
}
