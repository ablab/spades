///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "String.h"
#include "TokenizeString.h"
#include "Vec.h"

/*
 * Tokenize (with standard separators).
 */
int Tokenize( const String &a_string,
	      vec<String> &tokens )
{
  vec<char> separators;
  separators.push_back( ' ' );
  separators.push_back( '\t' );
  
  return Tokenize( a_string, separators, tokens );
}

/*
 * Tokenize (with one separator).
 */
int Tokenize( const String &a_string,
              const char sep,
	          vec<String> &tokens )
{
  vec<char> separators;
  separators.push_back(sep);
  return Tokenize( a_string, separators, tokens );
}


/*
 * Tokenize (with user defined separators).
 */
int Tokenize( const String &a_string,
	      const vec<char> &separators,
	      vec<String> &tokens )
{
  tokens.clear( );
  
  // Parse string.
  vec< pair<int, int> > token_interval;
  bool in_token = false;
  int token_start = 0;
  for (int ii=0; ii<(int)a_string.size( ); ii++) {
    bool matches_separator = false;
    for (int jj=0; jj<(int)separators.size( ); jj++)
      if ( separators[jj] ==a_string[ii] ) {
	matches_separator = true;
	break;
      }

    if ( matches_separator ) {
      if ( in_token ) {
	token_interval.push_back( pair<int, int>( token_start, ii ) );
	in_token = false;
      }
    }
    else {
      if ( !in_token ) {
	token_start = ii;
	in_token = true;
      }
    }
  }
  if ( in_token )
    token_interval.push_back( pair<int, int>( token_start, (int)a_string.size() ) );

  // Create tokens.
  tokens.reserve( token_interval.size( ) );
  for (int ii=0; ii<(int)token_interval.size( ); ii++) {
    int token_length = token_interval[ii].second - token_interval[ii].first;
    tokens.push_back( a_string.substr( token_interval[ii].first, token_length ) );
  }

  // Return number of tokens.
  return (int)token_interval.size( );
}



/*
 * Tokenize (with user defined separators).
 */
int TokenizeStrictly( const String &a_string,
		      const vec<char> &separators,
		      vec<String> &tokens )
{
  tokens.clear( );
  
  // Parse string.
  vec< pair<int, int> > token_interval;
  int token_start = 0;
  for (int ii=0; ii<(int)a_string.size( ); ii++) {
    bool matches_separator = false;
    for (int jj=0; jj<(int)separators.size( ); jj++)
      if ( separators[jj] == a_string[ii] ) {
	matches_separator = true;
	break;
      }

    if ( matches_separator ) {
      token_interval.push_back( pair<int, int>( token_start, ii ) );
      token_start = ii+1;
    }
  }
  token_interval.push_back( pair<int, int>( token_start, (int)a_string.size() ) );

  // Create tokens.
  tokens.reserve( token_interval.size( ) );
  for (int ii=0; ii<(int)token_interval.size( ); ii++) {
    int token_length = token_interval[ii].second - token_interval[ii].first;
    tokens.push_back( a_string.substr( token_interval[ii].first, token_length ) );
  }

  // Return number of tokens.
  return (int)token_interval.size( );
}

void RegexMatch(String regex, String &str, vec<String> &matches)
{
    char* command = new char[str.size() * 2];
    memset(command, 0x01, (str.size() + regex.size()) * 10);
    int number_of_parens = str.Freq("(");

    sprintf(command, "echo \"%s\" | perl -ne 'my @matches = %s; print join(\",\", @matches) . \"\\n\";'", str.c_str(), regex.c_str()); 
    String command_str(command);
    String result = LineOfOutput(command_str, false);
    
    vector<char> a_comma(1); a_comma[0] = ',';
    Tokenize(result, a_comma, matches);

    delete[] command;
}
