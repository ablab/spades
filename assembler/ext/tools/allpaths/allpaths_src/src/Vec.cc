///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include <algorithm>

#include <fcntl.h>
#include <fstream>
#include <math.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "String.h"
#include "system/ErrNo.h"
#include "system/Types.h"
#include "Vec.h"
#include "VecTemplate.h"

void PrettyPrint(ostream& o, const vec<int>& v, int max_items, String terminator)
{    int chars_in_line = 0;
     for ( vec<int>::size_type i = 0; i < v.size( ); i++ )
     {    chars_in_line += int(floor(log10( static_cast<float>( v[i] )))) + 2;
          if ( (int) i == max_items - 1 )
          {    if ( chars_in_line + 3 > 80 ) o << "\n" << v[i];
               else o << v[i];
               o << " ...";
               break;    }
          if ( chars_in_line > 80 )
          {    o << "\n";
               chars_in_line = int(floor(log10( static_cast<float>( v[i] )))) + 2;    }
          o << v[i] << " ";    }
     o << terminator;    }

void PrettyPrint(ostream& o, const vec<longlong>& v, int max_items, 
     String terminator)
{    int chars_in_line = 0;
     for ( vec<longlong>::size_type i = 0; i < v.size( ); i++ )
     {    chars_in_line += int(floor(log10( static_cast<float>( v[i] )))) + 2;
          if ( (int) i == max_items - 1 )
          {    if ( chars_in_line + 3 > 80 ) o << "\n" << v[i];
               else o << v[i];
               o << " ...";
               break;    }
          if ( chars_in_line > 80 )
          {    o << "\n";
               chars_in_line = int(floor(log10( static_cast<float>( v[i] )))) + 2;    }
          o << v[i] << " ";    }
     o << terminator;    }

void PrettyPrint(ostream& o, const vec<double>& v, int max_items, String terminator)
{    int chars_in_line = 0;
     for ( vec<double>::size_type i = 0; i < v.size( ); i++ )
     {    String s = ToString( v[i] );
          chars_in_line += s.size( ) + 1;
          if ( (int) i == max_items - 1 )
          {    if ( chars_in_line + 3 > 80 ) o << "\n" << v[i];
               else o << v[i];
               o << " ...";
               break;    }
          if ( chars_in_line > 80 )
          {    o << "\n";
               chars_in_line = s.size( ) + 1;    }
          o << v[i] << " ";    }
     o << terminator;    }

void PrettyPrint(ostream& o, const vec<TraceInt>& v, int max_items, String terminator)
{    int chars_in_line = 0;
     for ( vec<TraceInt>::size_type i = 0; i < v.size( ); i++ )
     {    chars_in_line += int(floor(log10( static_cast<float>( v[i] )))) + 2;
          if ( (int) i == max_items - 1 )
          {    if ( chars_in_line + 3 > 80 ) o << "\n" << v[i];
               else o << v[i];
               o << " ...";
               break;    }
          if ( chars_in_line > 80 )
          {    o << "\n";
               chars_in_line = int(floor(log10( static_cast<float>( v[i] )))) + 2;    }
          o << v[i] << " ";    }
     o << terminator;    }


istream& operator>>(istream& s, vec<String>& v)
{    int n;
     s >> n;
     v.resize(n);
     char c;
     s.get(c);
     for ( vec<String>::size_type i = 0; i < v.size( ); i++ )
          getline( s, v[i] );
     return s;    }

// istream& operator>>(istream& s, vec<String>& v)
// {    int n;
//      s >> n;
//      v.resize(n);
//      char c;
//      s.get(c);
//      for ( vec<Sting>::size_type i = 0; i < v.size( ); i++ )
//           getline( s, v[i] );
//      return s;    }

ostream& operator<<(ostream& s, const vec<double>& v)
{    s << v.size( ) << "\n";
     for ( vec<double>::size_type i = 0; i < v.size( ); i++ )
          s << v[i] << "\n";
     return s;    }

ostream& operator<<(ostream& s, const vec<unsigned short>& v)
{    s << v.size( ) << "\n";
     for ( vec<unsigned short>::size_type i = 0; i < v.size( ); i++ )
          s << v[i] << "\n";
     return s;    }

ostream& operator<<(ostream& s, const vec<int>& v)
{    s << v.size( ) << "\n";
     for ( vec<int>::size_type i = 0; i < v.size( ); i++ )
          s << v[i] << "\n";
     return s;    }

ostream& operator<<(ostream& s, const vec<longlong>& v)
{    s << v.size( ) << "\n";
     for ( vec<longlong>::size_type i = 0; i < v.size( ); i++ )
          s << v[i] << "\n";
     return s;    }

ostream& operator<<(ostream& s, const vec<float>& v)
{    s << v.size( ) << "\n";
     for ( vec<float>::size_type i = 0; i < v.size( ); i++ )
          s << v[i] << "\n";
     return s;    }

ostream& operator<<(ostream& s, const vec<String>& v)
{    s << v.size( ) << "\n";
     for ( vec<String>::size_type i = 0; i < v.size( ); i++ )
          s << v[i] << "\n";
     return s;    }

// ostream& operator<<(ostream& s, const vec<String>& v)
// {    s << v.size( ) << "\n";
//      for ( vec<String>::size_type i = 0; i < v.size( ); i++ )
//           s << v[i] << "\n";
//      return s;    }

template <>
void WriteAppend( const String& f, const vec<String>& v )
{
    ForceAssert( !IsRegularFile( f + ".gz" ) );
    longlong max_size_bound = longlong(10000000) * longlong(100000000);
    if ( !IsRegularFile(f) )
    {
        ForceAssertLt( (longlong) v.size( ), max_size_bound );
        Ofstream( out, f );
        out << setfill( '0' ) << setw(15) << v.size( ) << "\n";
        for ( longlong i = 0; i < (longlong) v.size( ); i++ )
            out << v[i] << "\n";
    }
    else
    {
        longlong n;
        {
            Ifstream( in, f );
            in >> n;
        }
        ForceAssertLt( n + (longlong) v.size( ), max_size_bound );
        ostrstream osize;
        osize << setfill( '0' ) << setw(15) << n + v.size( ) << "\n";
        int fd = Open( f, O_WRONLY );
        WriteBytes( fd, osize.str( ), 16 );
        close(fd);
        ofstream out( f.c_str( ), ios::app );
        for ( longlong i = 0; i < (longlong) v.size( ); i++ )
            out << v[i] << "\n";
    }
}

void PrintTabular( ostream& out, const vec< vec<String> >& rows, int sep,
     String justify)
{    int nrows = rows.size( ), ncols = 0;
     for ( int i = 0; i < nrows; i++ )
          ncols = std::max( ncols, (int) rows[i].size( ) );
     vec<int> maxcol;
     maxcol.resize_and_set( ncols, 0 );
     for ( int i = 0; i < (int) rows.size( ); i++ )
     {    for ( int j = 0; j < (int) rows[i].size( ); j++ )
               maxcol[j] = std::max( maxcol[j], (int) rows[i][j].size( ) );    }
     for ( int i = 0; i < (int) rows.size( ); i++ )
     {    for ( int j = 0; j < (int) rows[i].size( ); j++ )
          {    if ( j < (int) justify.size( ) && justify[j] == 'r' )
               {    for ( int k = 0; k < maxcol[j] - (int) rows[i][j].size( ); k++ )
                         out << " ";
                    out << rows[i][j];
                    if ( j < ncols - 1 )
                    {    for ( int k = 0; k < sep; k++ )
                              out << " ";    }    }
               else
               {    out << rows[i][j];
                    if ( j < ncols - 1 )
                    {    for ( int k = 0; 
                              k < maxcol[j] - (int) rows[i][j].size( ) + sep; k++ )
                              out << " ";    }    }    }
          out << "\n";    }    }

void PrintCSV(ostream& out, const vec< vec<String> >& rows)
{
    for (unsigned int i = 0; i < rows.size(); i++)
    {
        for (unsigned int j = 0; j < rows[i].size(); j++)
        {
            out << "\"" << rows[i][j] << "\"";
            if (j != rows[i].size()-1) { out << ","; }
            else { out << "\n"; }
        }
    }
}

bool IsAsciiVec( const String &filename )
{
  ifstream in( filename.c_str() );
  char c;
  while ( in )
  {
    in.get( c );
    if ( ! isdigit( c ) && ! isspace( c ) )
      return false;

    if ( c == '\n' ) 
      break;
  }
  return true;
}

longlong AsciiOrBinary0VecSize( const String& filename )
{
  String ns;
  Ifstream( in, filename );
  in >> ns;  
  ForceAssert( ns.IsInt( ) );
  longlong n = ns.Int( );
  return n;
}


int GetBinary2Or3ElementSize(const String & filename, int version = -1)
{
  if( version == -1 ) version = WhichBinaryFormat(filename);
  ForceAssert( version == 2 || version == 3 );
  const String& b_header = (version == 2 ? binary2_first_line : binary3_first_line);
  int header_size = (version == 2 ? binary2_header_size : binary3_header_size);

  String header;
  header.resize(header_size);

  if ( true )
  { FileReader fr(filename.c_str());
    fr.read( &header[0], header_size ); }

  for ( unsigned int i = 0; i < b_header.size( ); i++ ) {    
    if ( header[i] != b_header[i] ) {    
      return -1;    
    }    
  }   
  String elements = header.substr(binary_element_count_begin, 
                                  binary_element_count_size);
  elements = elements.Before(" ");
  longlong n = elements.Int();
  longlong size = FileSize(filename) - header.size();
  AssertEq(size % n, 0);
  return size / n;
}  

int GetBinary2ElementSize(const String & filename, bool strict) {
  return GetBinary2Or3ElementSize( filename, (strict ? 2 : -1) );
}
int GetBinary3ElementSize(const String & filename, bool strict) {
  return GetBinary2Or3ElementSize( filename, (strict ? 3 : -1) );
}


int WhichBinaryFormat( const String& filename ) {
  String word1, word2;
  int version;

  Ifstream(file, filename);
  file >> word1 >> word2 >> version;
  
  ForceAssertEq( word1, "binary" );
  ForceAssertEq( word2, "format" );

  return version;
}


void BinaryCat2Or3( const String & target, const String & source, 
		    int version) {
  ForceAssert( WhichBinaryFormat(target) == version );
  ForceAssert( WhichBinaryFormat(source) == version );
  longlong ns = BinaryNumElements(source);
  longlong nt = BinaryNumElements(target);

  //adjust the number of records
  int fdt = OpenForWrite(target);
  lseek(fdt, binary_element_count_begin, SEEK_SET);
  BinaryWriteSize(fdt, ns + nt);
  lseek(fdt, 0, SEEK_END);

  //start reading at the beginning of the data.
  FileReader fr(source.c_str());
  int offset = version == 2 ? binary2_header_size : binary3_header_size;
  fr.seek(offset);

  //add in all the new records, 1024  bytes at a time so we don't overwhelm
  //the memory.
  const int BLOCK=1024;
  char c[BLOCK];
  while (true) {
    int nread = fr.readSome( c, BLOCK );
    if ( !nread ) break;
    WriteBytes(fdt, c, nread);
  }
  close(fdt);
}
  
void BinaryCat2( const String& target, const String & source ) {
  BinaryCat2Or3( target, source, 2);
}

void BinaryCat3( const String& target, const String & source ) {
  BinaryCat2Or3( target, source, 3);
}

longlong BinaryNumElements( const String & filename) {
  int version = WhichBinaryFormat(filename);
  ForceAssert(2 == version || 3 == version);

  Ifstream(is, filename);
  String junk;
  getline(is, junk);
  longlong elements;
  is >> elements;
  ForceAssert(is.good());

  return elements;
}

template< > vec<double>::vec( const String& s )
{    ForceAssert( s.Contains( "{", 0 ) );
     ForceAssert( s.Contains( "}", -1 ) );
     String t = s;
     t = t.After( "{" );
     t = t.Before( "}" );
     while(1)
     {    if ( t.Contains( "," ) )
          {    push_back( t.Before( "," ).Double( ) );
               t = t.After( "," );    }
          else 
          {    push_back( t.Double( ) );
               break;    }    }    }


template< > vec<int>::vec( const String& s )
{    ForceAssert( s.Contains( "{", 0 ) );
     ForceAssert( s.Contains( "}", -1 ) );
     String t = s;
     t = t.After( "{" );
     t = t.Before( "}" );
     while(1)
     {    if ( t.Contains( "," ) )
          {    push_back( t.Before( "," ).Int( ) );
               t = t.After( "," );    }
          else 
          {    push_back( t.Int( ) );
               break;    }    }    }


void BinaryWrite( int fd, const vec< String >& v ) {
  BinaryWriteComplex( fd, v );
}

void BinaryRead( int fd, vec< String >& v ) {
  BinaryReadComplex( fd, v );
}


BINARY2_DEF(char);
BINARY2_DEF(unsigned char);
BINARY2_DEF(int);
BINARY2_DEF(longlong);
BINARY2_DEF(float);

BINARY3_DEF(char);
BINARY3_DEF(unsigned char);
BINARY3_DEF(int);
BINARY3_DEF(longlong);
BINARY3_DEF(float);
BINARY3_DEF(double);
