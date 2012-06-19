#ifndef UTIL
#define UTIL
#include <stdlib.h>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#define MAX_STRING_LENGTH_UTIL 256

using namespace std;

class Utils {
public:
  static string trimSpace( const string& s ){
    if( s.length() == 0 )
      return s;
    int b = s.find_first_not_of( " \t\r\n" );
    int e = s.find_last_not_of( " \t\r\n" );
    if( b == -1 )
      return "";
    return s.substr( b, e-b+1 );
  }

  static vector<string> split( string str, const string& delim ){
    vector<string> vec;
    uint cutPos;
    while( (cutPos = str.find_first_of(delim)) != (unsigned)-1 ){
      if( cutPos > 0 ){
	vec.push_back(str.substr(0, cutPos));
      }
      str = str.substr(cutPos + 1);
    }
    if( str.length()>0 ){
      vec.push_back(str);
    }
    return vec;
  }

  static string join( const vector<string>& cols, const string& delim ){
    string str = cols.at(0);
    for( uint i=1 ; i<cols.size() ; ++i ){
      str += (delim + cols.at(i));
    }
    return str;
  }

  static string removeColon( const string& str ){
    vector<string> cols = split(str, ":");
    return join(cols, "^");
  }

  static string itoa( int i ){
    char c_str[MAX_STRING_LENGTH_UTIL];
    sprintf( c_str, "%d", i );
    return (string)c_str;
  }

  static string ltoa( long l ){
    char c_str[MAX_STRING_LENGTH_UTIL];
    sprintf( c_str, "%ld", l );
    return (string)c_str;
  }

  static string dtoa( double d ){
    char c_str[MAX_STRING_LENGTH_UTIL];
    sprintf( c_str, "%lf", d );
    return (string)c_str;
  }
  
  static uint max( const vector<uint>& vec ){
    uint max = vec.at(0);
    for( uint i=1 ; i<vec.size() ; ++i )
      if( vec.at(i)>max )
	max = vec.at(i);
    return max;
  }

  static void fileopen( ofstream& ofs, const string& filename, 
			const string& prefix="" ){
    try {
      ofs.open( filename.c_str() );
      if( ofs.fail() )
	throw (string)"File open error: " + filename;
    } catch(const string errorMessage){
      cerr << prefix << errorMessage << endl;
      exit(-1);
    }
  }

  static void fileopen( ifstream& ifs, const string& filename, 
			const string& prefix=""){
    try {
      ifs.open( filename.c_str() );
      if( ifs.fail() )
	throw (string)"File open error: " + filename;
    } catch(const string& errorMessage){
      cerr << prefix << errorMessage << endl;
      exit(-1);
    }
  }

  static void mkdir(const string& dir){
    system( ("mkdir -p "+dir).c_str() );
  }
};

#endif // UTIL
