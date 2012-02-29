///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include <sys/types.h>

#include "CoreTools.h"
#include "FastIfstream.h"
#include "ShortVector.h"



/// Rewinds the currently open file back to start position;
/// has no effect if the stream is not associated with a file
/// (i.e. default constructor was used and no Open() was issued yet)

void fast_ifstream::Rewind() {
  fr_.seek(0);

  // discard all the data from the read buffer:
  buf_ptr_ = 0 ; 
  buf_top_ = 0 ;
  fail_ = False;
}

void fast_ifstream::get( char& c, Bool advance )
{    if (fail_) 
          FatalErr("get called on failed fast_ifstream for "
                      << fr_.getFilename() << "." );
     if ( buf_ptr_ == buf_top_ )
     {    buf_ptr_ = 0;
          buf_top_ = fr_.readSome(buf_, BUFFER_SIZE);
          if ( buf_top_ == 0 ) 
          {    fail_ = True;
               return;    }    }
     c = buf_[buf_ptr_];
     if (advance) ++buf_ptr_;    }

void fast_ifstream::peek( char& c )
{    get( c, False );    }


void getline( fast_ifstream& in, String& s )
{    
     if ( in.fail( ) )
          FatalErr( "getline called on failed fast_ifstream for " 
               << in.fr_.getFilename() << "." );

     // The fast case.

     for ( int i = in.buf_ptr_; i < in.buf_top_; i++ )
     {    if ( in.buf_[i] == '\n' )
          {    in.buf_[i] = 0;
               s.Set( in.buf_ + in.buf_ptr_, i - in.buf_ptr_ );
               in.buf_[i] = '\n';    
               in.buf_ptr_ = i + 1;
               return;    }    }

     // The slow case.
     char * start = 0;
     char * end = 0;
     in.ReadSlow(start, end);
     s.Set( start, end - start );    
}

void fast_ifstream::ReadSlow(char * & start, char * & end) {
  static avector<char> linebuf(100);
  char* x = linebuf.x;

  char* top = linebuf.x + linebuf.length;
  const int buf_mult = 2;
  while(1) {    
    if ( x == top ) {   
      int xpos = x - linebuf.x;
      linebuf.resize( linebuf.length * buf_mult );
      x = linebuf.x + xpos;
      top = linebuf.x + linebuf.length;    
    }
    get( *x );
    if ( fail( ) ) break;
    if ( *x == '\n' ) break;
    ++x;    }
  if ( x == top ) {    
    int xpos = x - linebuf.x;
    linebuf.resize( linebuf.length * buf_mult );
    x = linebuf.x + xpos;    
  }
  *x = 0;
  start = linebuf.x;
  end = x;
}

Bool getline_if_match( fast_ifstream& in, String& s, const String& begin )
{    
     if ( in.fail( ) )
          FatalErr( "getline_if_match called on failed fast_ifstream for " 
               << in.fr_.getFilename() << "." );

     // The fast case.

     for ( int i = in.buf_ptr_; i < in.buf_top_; i++ )
     {    if ( in.buf_[i] == '\n' )
          {    Bool match = True;
               unsigned int bp = 0;
               for ( int j = in.buf_ptr_; j < i; j++ )
               {    if ( bp == begin.size( ) ) break;
                    if ( in.buf_[j] != begin[bp] ) 
                    {    match = False;
                         break;    }
                    ++bp;    }
               if ( bp < begin.size( ) ) match = False;
               if (match)
               {    in.buf_[i] = 0;
                    s.Set( in.buf_ + in.buf_ptr_, i - in.buf_ptr_ );
                    in.buf_[i] = '\n';    }
               in.buf_ptr_ = i + 1;
               return match;    }    }

     // The slow case.
     char * start = 0;
     char * end = 0;
     in.ReadSlow(start, end);

     Bool match = True;
     unsigned int bp = 0;
     for ( int j = 0; j < end - start; j++ )
     {    if ( bp == begin.size( ) ) break;
          if ( start[j] != begin[bp] ) 
          {    match = False;
               break;    }
          ++bp;    }
     if ( bp < begin.size( ) ) match = False;
     if (match) s.Set( start, end - start );
     return match;    }

void fast_ifstream::fill_buffer( )
{    if ( buf_ptr_ > 0 )
     {    for ( int i = buf_ptr_; i < buf_top_; i++ )
               buf_[ i - buf_ptr_ ] = buf_[i];    }
     buf_top_ -= buf_ptr_;
     buf_ptr_ = 0;
     int answer = fr_.readSome(buf_ + buf_top_, BUFFER_SIZE - buf_top_);
     if ( answer == 0 )
     {    fail_ = True;
          return;    }
     buf_top_ += answer;
     return;    }

Bool get_to( fast_ifstream& in, String& s, const String& tail )
{    
     if ( in.fail( ) )
          FatalErr( "get_to called on failed fast_ifstream for " 
               << in.fr_.getFilename() << "." );

     static vec<char> linebuf;
     linebuf.resize(0);
     int nt = tail.size( ), ls;

     ForceAssert( nt < fast_ifstream::BUFFER_SIZE );
     if ( in.buf_top_ - nt < in.buf_ptr_ )
     {    in.fill_buffer( );
          if ( in.fail( ) ) return False;    }

     restart:

     for ( int i = in.buf_ptr_; i < in.buf_top_ - nt; i++ )
     {    int j;
          for ( j = 0; j < nt; j++ )
               if ( in.buf_[ i + j ] != tail[j] ) break;
          if ( j == nt )
          {    if ( linebuf.size( ) == 0 )
               {    char c = in.buf_[ i + nt ];
                    in.buf_[ i + nt ] = 0;
                    s.Set( in.buf_ + in.buf_ptr_, i + nt - in.buf_ptr_ );
                    in.buf_[ i + nt ] = c;    }
               else 
               {    int ls = linebuf.size( );
                    linebuf.resize( ls + i + nt - in.buf_ptr_ + 1 );
                    memcpy( &linebuf[ls], &in.buf_[in.buf_ptr_], 
                         i + nt - in.buf_ptr_ );
                    linebuf.back( ) = 0;
                    s.Set( &linebuf[0], linebuf.size( ) - 1 );    }
               in.buf_ptr_ += i + nt - in.buf_ptr_;
               return True;    }    }

     ls = linebuf.size( );
     linebuf.resize( ls + in.buf_top_ - nt - in.buf_ptr_ );
     memcpy( &linebuf[ls], &in.buf_[in.buf_ptr_], in.buf_top_ - nt - in.buf_ptr_ );
     for ( int i = 0; i < nt; i++ )
          in.buf_[i] = in.buf_[ in.buf_top_ - nt + i ];
     in.buf_ptr_ = 0;
     in.buf_top_ = nt;
     in.fill_buffer( );
     if ( in.fail( ) ) return False;
     goto restart;    }

fast_pipe_ifstream::fast_pipe_ifstream( String command )
{    command_ = command;
     file_ = popen( command.c_str( ), "r" );
     if ( file_ == 0 )
          FatalErr( "Couldn't open pipe for command \" << command << \"." );
     buf_ptr_ = 0;
     buf_top_ = 0;
     fail_ = False;    }

fast_pipe_ifstream::~fast_pipe_ifstream( )
{    if ( pclose(file_) < 0 )
          FatalErr( "Attempt to close " << command_ << " failed." );    }

void fast_pipe_ifstream::get( char& c, Bool advance )
{    if (fail_) 
          FatalErr( "get called on failed fast_pipe_ifstream for " 
               << command_ << "." );
     if ( buf_ptr_ == buf_top_ )
     {    buf_ptr_ = 0;
          buf_top_ = fread( buf_, 1, BUFFER_SIZE, file_ );
          if ( buf_top_ < 0 )
               FatalErr( "Attempt to read " << command_ << " failed." );
          if ( buf_top_ == 0 ) 
          {    fail_ = True;
               return;    }    }
     c = buf_[buf_ptr_];
     if (advance) ++buf_ptr_;    }

void fast_pipe_ifstream::peek( char& c )
{    get( c, False );    }

void getline( fast_pipe_ifstream& in, String& s )
{    
     if ( in.fail( ) )
          FatalErr( "getline called on failed fast_pipe_ifstream for " 
               << in.command_ << "." );

     // The fast case.

     for ( int i = in.buf_ptr_; i < in.buf_top_; i++ )
     {    if ( in.buf_[i] == '\n' )
          {    in.buf_[i] = 0;
               s.Set( in.buf_ + in.buf_ptr_, i - in.buf_ptr_ );
               in.buf_[i] = '\n';    
               in.buf_ptr_ = i + 1;
               return;    }    }

     // The slow case.

     static avector<char> linebuf(0);
     char* x = linebuf.x;
     char* top = linebuf.x + linebuf.length;
     const int buf_incr = 100;
     while(1)
     {    if ( x == top ) 
          {    linebuf.resize( linebuf.length + buf_incr );
               x = linebuf.x + linebuf.length - buf_incr;
               top = linebuf.x + linebuf.length;    }
          in.get( *x );
          if ( in.fail( ) ) break;
          if ( *x == '\n' ) break;
          ++x;    }
     if ( x == top ) 
     {    linebuf.resize( linebuf.length + buf_incr );
          x = linebuf.x + linebuf.length - buf_incr;    }
     *x = 0;
     s.Set( linebuf.x, x - linebuf.x );    }

Bool getline_if_match( fast_pipe_ifstream& in, String& s, const String& begin )
{    
     if ( in.fail( ) )
          FatalErr( "getline_if_match called on failed fast_pipe_ifstream for " 
               << in.command_ << "." );

     // The fast case.

     for ( int i = in.buf_ptr_; i < in.buf_top_; i++ )
     {    if ( in.buf_[i] == '\n' )
          {    Bool match = True;
               unsigned int bp = 0;
               for ( int j = in.buf_ptr_; j < i; j++ )
               {    if ( bp == begin.size( ) ) break;
                    if ( in.buf_[j] != begin[bp] ) 
                    {    match = False;
                         break;    }
                    ++bp;    }
               if ( bp < begin.size( ) ) match = False;
               if (match)
               {    in.buf_[i] = 0;
                    s.Set( in.buf_ + in.buf_ptr_, i - in.buf_ptr_ );
                    in.buf_[i] = '\n';    }
               in.buf_ptr_ = i + 1;
               return match;    }    }

     // The slow case.

     static avector<char> linebuf(0);
     char* x = linebuf.x;
     char* top = linebuf.x + linebuf.length;
     const int buf_incr = 100;
     while(1)
     {    if ( x == top ) 
          {    linebuf.resize( linebuf.length + buf_incr );
               x = linebuf.x + linebuf.length - buf_incr;
               top = linebuf.x + linebuf.length;    }
          in.get( *x );
          if ( in.fail( ) ) break;
          if ( *x == '\n' ) break;
          ++x;    }
     if ( x == top ) 
     {    linebuf.resize( linebuf.length + buf_incr );
          x = linebuf.x + linebuf.length - buf_incr;    }
     *x = 0;
     Bool match = True;
     unsigned int bp = 0;
     for ( int j = 0; j < x - linebuf.x; j++ )
     {    if ( bp == begin.size( ) ) break;
          if ( linebuf(j) != begin[bp] ) 
          {    match = False;
               break;    }
          ++bp;    }
     if ( bp < begin.size( ) ) match = False;
     if (match) s.Set( linebuf.x, x - linebuf.x );
     return match;    }

void fast_pipe_ifstream::fill_buffer( )
{    if ( buf_ptr_ > 0 )
     {    for ( int i = buf_ptr_; i < buf_top_; i++ )
               buf_[ i - buf_ptr_ ] = buf_[i];    }
     buf_top_ -= buf_ptr_;
     buf_ptr_ = 0;
     int answer = fread( buf_ + buf_top_, 1, BUFFER_SIZE - buf_top_, file_ );
     if ( answer < 0 ) FatalErr( "Attempt to read " << command_ << " failed." );
     if ( answer == 0 )
     {    fail_ = True;
          return;    }
     buf_top_ += answer;
     return;    }

Bool get_to( fast_pipe_ifstream& in, String& s, const String& tail )
{    
     if ( in.fail( ) )
          FatalErr( "get_to called on failed fast_pipe_ifstream for " 
               << in.command_ << "." );

     static vec<char> linebuf;
     linebuf.resize(0);
     int nt = tail.size( ), ls;

     ForceAssert( nt < fast_pipe_ifstream::BUFFER_SIZE );
     if ( in.buf_top_ - nt < in.buf_ptr_ )
     {    in.fill_buffer( );
          if ( in.fail( ) ) return False;    }

     restart:

     for ( int i = in.buf_ptr_; i < in.buf_top_ - nt; i++ )
     {    int j;
          for ( j = 0; j < nt; j++ )
               if ( in.buf_[ i + j ] != tail[j] ) break;
          if ( j == nt )
          {    if ( linebuf.size( ) == 0 )
               {    char c = in.buf_[ i + nt ];
                    in.buf_[ i + nt ] = 0;
                    s.Set( in.buf_ + in.buf_ptr_, i + nt - in.buf_ptr_ );
                    in.buf_[ i + nt ] = c;    }
               else 
               {    int ls = linebuf.size( );
                    linebuf.resize( ls + i + nt - in.buf_ptr_ + 1 );
                    memcpy( &linebuf[ls], &in.buf_[in.buf_ptr_], 
                         i + nt - in.buf_ptr_ );
                    linebuf.back( ) = 0;
                    s.Set( &linebuf[0], linebuf.size( ) - 1 );    }
               in.buf_ptr_ += i + nt - in.buf_ptr_;
               return True;    }    }

     ls = linebuf.size( );
     linebuf.resize( ls + in.buf_top_ - nt - in.buf_ptr_ );
     memcpy( &linebuf[ls], &in.buf_[in.buf_ptr_], in.buf_top_ - nt - in.buf_ptr_ );
     for ( int i = 0; i < nt; i++ )
          in.buf_[i] = in.buf_[ in.buf_top_ - nt + i ];
     in.buf_ptr_ = 0;
     in.buf_top_ = nt;
     in.fill_buffer( );
     if ( in.fail( ) ) return False;
     goto restart;    }

void ReadParagraphs( fast_ifstream& in, vec<String>& paragraphs )
{    paragraphs.clear( );
     String line, par;
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) break;
          if ( line.size( ) != 0 )
          {    if ( par.size( ) > 0 ) par += "\n";
               par += line;    }
          else if ( par.size( ) > 0 )
          {    paragraphs.push_back(par);
               par = "";    }    }
     if ( par.size( ) > 0 ) paragraphs.push_back(par);    }

void ReadParagraphs( fast_pipe_ifstream& in, vec<String>& paragraphs )
{    paragraphs.clear( );
     String line, par;
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) break;
          if ( line.size( ) != 0 )
          {    if ( par.size( ) > 0 ) par += "\n";
               par += line;    }
          else if ( par.size( ) > 0 )
          {    paragraphs.push_back(par);
               par = "";    }    }
     if ( par.size( ) > 0 ) paragraphs.push_back(par);    }
