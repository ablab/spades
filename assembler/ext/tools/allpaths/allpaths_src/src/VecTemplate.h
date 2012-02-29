///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef VEC_TEMPLATE
#define VEC_TEMPLATE

#include "String.h"
#include "system/Types.h"
#include "Vec.h"

const int binary2_header_size = 62;
const int binary3_header_size = 96;
const int binary_element_count_begin = 34;
const int binary_element_count_size = 13;
const longlong binary2_or3_max_size = 999999999999;  // 12 character long as ASCII
const unsigned int binary_element_count_end = binary_element_count_begin + binary_element_count_size;

const static String binary2_first_line( "binary format 2, header = 3 lines\n" );
const static String binary3_first_line( "binary format 3, header = 4 lines\n" );

const static String binary_little_endian_line( "\nlittle endian\n" );
const static String binary_big_endian_line(    "\nbig endian   \n" );

const static String binary3_padding_line( "padding to make long word aligned\n" );

inline
void BinaryWriteSize( const int fd, const longlong n ) {
  String length = ToString(n);
  int k = length.size( );
  ForceAssertLe( k, binary_element_count_size );
  length.resize( binary_element_count_size );
  for ( int i = k; i < binary_element_count_size; i++ )
    length[i] = ' ';
  WriteBytes( fd, length.c_str( ), binary_element_count_size );
}

template <class T>
void
BinaryWrite2Or3( const String& filename, const vec<T>& v,
                 const int version )
{    
  ForceAssert( version == 2 || version == 3 );
  longlong n = v.size( );
  ForceAssertLe(n, binary2_or3_max_size);
  static longlong ten_trillion = (longlong) 10000000 * (longlong) 10000000;
  ForceAssertLt( n, ten_trillion );
  Remove(filename);
  int fd = OpenForWrite(filename);
  if ( version == 2 )
    WriteBytes( fd, binary2_first_line.c_str( ), binary2_first_line.size( ) );
  else
    WriteBytes( fd, binary3_first_line.c_str( ), binary3_first_line.size( ) );
  BinaryWriteSize( fd, n );
  ForceAssertEq( binary_little_endian_line.size( ), binary_big_endian_line.size( ) );
#ifdef Little_Endian
  {    WriteBytes( fd, binary_little_endian_line.c_str( ), binary_little_endian_line.size( ) );    }
#else
  {    WriteBytes( fd, binary_big_endian_line.c_str( ), binary_big_endian_line.size( ) );    }
#endif
  if ( version == 3 )
    WriteBytes( fd, binary3_padding_line.c_str(), binary3_padding_line.size() );
  if ( n > 0 ) 
    WriteBytes( fd, &v[0], (longlong) sizeof(T) * n );
  close(fd);    
}

template <class T>
void
BinaryWrite2( const String& filename, const vec<T>& v ) {
  BinaryWrite2Or3( filename, v, 2 );
}

template <class T>
void
BinaryWrite3( const String& filename, const vec<T>& v ) {
  BinaryWrite2Or3( filename, v, 3 );
}


template <class T>
void
CheckHeader2Or3( FileReader const& fr, longlong& n, const int version )
{
  ForceAssert( version == 2 || version == 3 );

  String header;
  if ( version == 2 ) {
    header.resize( binary2_header_size );
    fr.read( &header[0], binary2_header_size );
    for ( unsigned int i = 0; i < binary2_first_line.size( ); i++ )
      if ( header[i] != binary2_first_line[i] )
        FatalErr( "Binary read 2 of " << fr.getFilename() << " failed: "
                  << "first line doesn't match expected value." );
  }
  else {
    header.resize( binary3_header_size );
    fr.read( &header[0], binary3_header_size );
    for ( unsigned int i = 0; i < binary3_first_line.size( ); i++ )
      if ( header[i] != binary3_first_line[i] )
        FatalErr( "Binary read 3 of " << fr.getFilename() << " failed: "
                  << "first line doesn't match expected value." );
  }

  Bool little = True, big = True;
  for ( unsigned int i = 0; i < binary_little_endian_line.size( ); i++ )
    if ( header.c_str()[ binary_element_count_end + i ] != binary_little_endian_line[i] ) {
      little = False;
      break;
    }
  for ( unsigned int i = 0; i < binary_big_endian_line.size( ); i++ )
    if ( header[ binary_element_count_end + i ] != binary_big_endian_line[i] ) {    
      big = False;
      break;
    }
  if ( !little && !big ) 
    FatalErr( "Binary read " << version << " of " << fr.getFilename()
             << " failed: can't determine endian setting from header." );
#ifdef Little_Endian
  if ( !little )
    FatalErr( "Binary read " << version << " of " << fr.getFilename()
              << " failed because file was written on a big endian "
              "architecture,\nand read back in on this little endian "
              "architecture.\nUnfortunately, this is not possible at present.");
#else
  if ( !big )
    FatalErr( "Binary read " << version << " of " << fr.getFilename()
              << " failed because file was written on a little endian "
              "architecture,\nand read back in on this big endian "
              "architecture.\nUnfortunately, this is not possible at present.");
#endif

  unsigned int d;
  for ( d = binary_element_count_begin; d < binary_element_count_end; d++ ) {
    if ( header.c_str()[d] == ' ' ) break;
    if ( !isdigit( header.c_str()[d] ) )
      FatalErr( "Binary read 2 of " << fr.getFilename() << " failed: "
                << "didn't find record count where it should be." );   
  }
  for ( unsigned int d2 = d + 1; d2 < binary_element_count_end; d2++ ) {
    if ( header[d] != ' ' )
      FatalErr( "Binary read 2 of " << fr.getFilename() << " failed: "
                << "didn't find white space where expected." );
  }
  String ns = header.substr( binary_element_count_begin, d - binary_element_count_begin );
  n = ns.Int( ); 
  longlong N = fr.getSize();
  if ( N != (longlong) header.size( ) + n * (longlong) sizeof(T) )
    FatalErr( "Binary read " << version << " of " << fr.getFilename()
              << " failed:\n" << "header size = " << header.size( )
              << ", filesize = " << N << ", record count = " << n
              << ", record size = " << sizeof(T) << "." ); 
}

template <class T>
void
CheckHeader2( FileReader const& fr, longlong& n ) {
  CheckHeader2Or3<T>( fr, n, 2 );
}

template <class T> 
void
CheckHeader3( FileReader const& fr, longlong& n ) {
  CheckHeader2Or3<T>( fr, n, 3 );
}


template <class T> 
void
BinaryRead2Or3( const String& filename, vec<T>& v, int version = -1, 
     const bool append = False )
{
  if( version == -1 ) version = WhichBinaryFormat(filename);
  ForceAssert( version == 2 || version == 3 );
  FileReader fr(filename.c_str());
  longlong n;
  if ( version == 2 )
    CheckHeader2<T>( fr, n );
  else
    CheckHeader3<T>( fr, n );
  longlong start = ( append ? v.size( ) : 0 );
  v.resize( append ? n + v.size( ) : n );
  if ( n > 0 ) fr.read( &v[start], (longlong) sizeof(T) * n );
}

template <class T> 
void
BinaryRead2( const String& filename, vec<T>& v, bool strict, const Bool append ) {
  BinaryRead2Or3<T>( filename, v, (strict ? 2 : -1), append );
}

template <class T> 
void
BinaryRead3( const String& filename, vec<T>& v, bool strict, const Bool append ) {
  BinaryRead2Or3<T>( filename, v, (strict ? 3 : -1), append );
}


template <class T> 
void
BinaryReadSubset2Or3( const String& filename, const vec<int>& ids, vec<T>& v, Bool append, 
                      int version = -1 )
{    
  if( version == -1 ) version = WhichBinaryFormat(filename);
  ForceAssert( version == 2 || version == 3 );
  FileReader fr(filename.c_str());
  longlong n;
  if ( version == 2 )
    CheckHeader2<T>( fr, n );
  else
    CheckHeader3<T>( fr, n );
  int newsize = ids.size( );
  if (append) newsize += v.size( );
  int start = ( append ? v.isize() : 0 );
  v.resize(newsize);
  const int header_size = ( version == 2 ? binary2_header_size : binary3_header_size );
  for ( int i = 0; i < ids.isize( ); i++ ) {
    fr.seek( header_size + ids[i] * sizeof(T) );
    fr.read( &v[start+i], sizeof(T) );
  }
}

template <class T> 
void
BinaryReadSubset2( const String& filename, const vec<int>& ids, vec<T>& v, 
		   Bool append, bool strict ) {
  BinaryReadSubset2Or3<T>( filename, ids, v, append, (strict ? 2 : -1) );
}

template <class T> 
void
BinaryReadSubset3( const String& filename, const vec<int>& ids, vec<T>& v, 
		   Bool append, bool strict ) {
  BinaryReadSubset2Or3<T>( filename, ids, v, append, (strict ? 3 : -1) );
}


template <class T>
void
BinaryReadRange2Or3( const String& filename, longlong from, longlong to, vec<T>& v, 
                     int version = -1 )
{    
  if( version == -1 ) version = WhichBinaryFormat(filename);
  ForceAssert( version == 2 || version == 3 );
  ForceAssert( from <= to );
  FileReader fr(filename.c_str());
  longlong n;
  if ( version == 2 )
    CheckHeader2<T>( fr, n );
  else
    CheckHeader3<T>( fr, n );
  v.resize( to - from );
  const int header_size = ( version == 2 ? binary2_header_size : binary3_header_size );
  fr.seek( header_size + from * sizeof(T) );
  if ( to > from ) fr.read( &v[0], (to - from) * sizeof(T) );
}

template <class T>
void
BinaryReadRange2( const String& filename, longlong from, longlong to, 
		  vec<T>& v, bool strict ) {
  BinaryReadRange2Or3<T>( filename, from, to, v, (strict ? 2 : -1) );
}

template <class T> 
void
BinaryReadRange3( const String& filename, longlong from, longlong to, 
		  vec<T>& v, bool strict ) {
  BinaryReadRange2Or3<T>( filename, from, to, v, (strict ? 3 : -1) );
}


template<class T> 
longlong 
BinarySize2Or3( const String& filename, int version = -1 )
{
  if( version == -1 ) version = WhichBinaryFormat(filename);
  ForceAssert( version == 2 || version == 3 );
  FileReader fr(filename.c_str());
  longlong n;
  if ( version == 2 )
    CheckHeader2<T>( fr, n );
  else
    CheckHeader3<T>( fr, n );
  return n;
}

template<class T> 
longlong 
BinarySize2( const String& filename, bool strict ) {
  return BinarySize2Or3<T>( filename, (strict ? 2 : -1) );
}

template<class T> 
longlong 
BinarySize3( const String& filename, bool strict ) {
  return BinarySize2Or3<T>( filename, (strict ? 3 : -1) );
}



template <typename T>
void Binary3Writer<T>::Open( const String& filename, bool keep_open ) {
  m_filename = filename;
  m_keep_open = keep_open;
  Remove( filename );
  vec<T> emptyVec;
  BinaryWrite3( filename, emptyVec );
  m_fd = ::Open( filename, O_WRONLY );
  off_t currpos = lseek( m_fd, 0, SEEK_END );
  ForceAssertEq( binary3_header_size, currpos );
  m_objectCount = 0;
  if( ! m_keep_open )
    ::Close(m_fd);
}

template <typename T>
Binary3Writer<T>::~Binary3Writer() {
  this->Close();
}

template <typename T>
void Binary3Writer<T>::Write( const T& object ) {
  if( ! m_keep_open ) {
    m_fd = ::Open( m_filename, O_WRONLY ^ O_APPEND );
  }

  WriteBytes( m_fd, (char*) &object, sizeof(T) );
  ++m_objectCount;

  if( ! m_keep_open ) {
    ::Close(m_fd);
  }
}

template <typename T>
void Binary3Writer<T>::WriteMultiple( const vec<T>& objects ) {
  if( objects.empty() ) return;

  if( ! m_keep_open ) {
    m_fd = ::Open( m_filename, O_WRONLY ^ O_APPEND );
  }

  WriteBytes( m_fd, (char*) &objects[0], sizeof(T) * longlong(objects.size()) );
  m_objectCount += objects.size();

  if( ! m_keep_open ) {
    ::Close(m_fd);
  }
}

template <typename T>
void Binary3Writer<T>::Close() {
  if ( m_fd < 0 ) 
    return;

  if( ! m_keep_open ) {
    m_fd = ::Open( m_filename, O_WRONLY );
    lseek( m_fd, 0, SEEK_END );    
  }
    
  off_t currpos = lseek( m_fd, binary_element_count_begin, SEEK_SET );
  ForceAssertEq( currpos, binary_element_count_begin );
  BinaryWriteSize( m_fd, m_objectCount );
  ::Close( m_fd );
  ForceAssertLe(m_objectCount, binary2_or3_max_size);
  m_fd = -1;
  m_objectCount = 0;
}

template <typename T>
Binary3Iter<T>::Binary3Iter( const String& filename, T* p_to_fill, 
			     longlong max_memory )
: mFR(filename.c_str())
{
  m_maxsize = max_memory / sizeof(T);

  CheckHeader3<T>( mFR, m_globalSize );

  m_globalIndex = m_localIndex = 0;

  FillBuffer();
  
  if( ! m_data.empty() )
    (*p_to_fill) = m_data[0];
}

template <typename T>
void Binary3Iter<T>::Next( T* p_to_fill ) {
  if( ++m_globalIndex >= m_globalSize ) return;
  if( ++m_localIndex == m_data.size() ) FillBuffer();
  (*p_to_fill) = m_data[m_localIndex];
}

template <typename T>
void Binary3Iter<T>::FillBuffer() {
  ForceAssertEq( m_localIndex, m_data.size() );
  m_data.resize( min( longlong(m_maxsize), m_globalSize - m_globalIndex ) );
  if( m_data.size() > 0 )
    mFR.read( &m_data[0], sizeof(T) * longlong(m_data.size()) );
  m_localIndex = 0;
}


#define BINARY2_DEF(T)                                                         \
     template void BinaryWrite2( const String& filename, const vec<T>& v );    \
     template void BinaryRead2( const String& filename, vec<T>& v, bool,       \
	  const Bool );                                                        \
     template void BinaryReadSubset2( const String& filename,                  \
          const vec<int>& ids, vec<T>& v, Bool, bool );                        \
     template void BinaryReadRange2( const String& filename,                   \
          longlong from, longlong to, vec<T>& v, bool );                       \
     template longlong BinarySize2<T>( const String& filename, bool )

#define BINARY3_DEF(T)                                                         \
     template void BinaryWrite3( const String& filename, const vec<T>& v );    \
     template void BinaryRead3( const String& filename, vec<T>& v, bool,       \
          const Bool  );                                                       \
     template void BinaryReadSubset3( const String& filename,                  \
          const vec<int>& ids, vec<T>& v, Bool, bool );                        \
     template void BinaryReadRange3( const String& filename,                   \
          longlong from, longlong to, vec<T>& v, bool );                       \
     template longlong BinarySize3<T>( const String& filename, bool );         \
     template Binary3Writer<T>::Binary3Writer( const String&, bool );          \
     template Binary3Writer<T>::~Binary3Writer();                              \
     template void Binary3Writer<T>::Write( const T& );                        \
     template void Binary3Writer<T>::WriteMultiple( const vec<T>& );           \
     template void Binary3Writer<T>::Open( const String&, bool );              \
     template void Binary3Writer<T>::Close();                                  \
     template Binary3Iter<T>::Binary3Iter( const String&, T*, longlong );      \
     template void Binary3Iter<T>::Next( T* );                                 \
     template void Binary3Iter<T>::FillBuffer()


#endif // #define VEC_TEMPLATE
