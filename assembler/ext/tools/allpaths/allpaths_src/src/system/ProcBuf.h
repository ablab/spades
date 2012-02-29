///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#if ! defined( __GNUC__ ) || __GNUC__ > 2

#ifndef PROCBUF_H
#define PROCBUF_H

#include <iosfwd>
#include <streambuf>
#include <string>
#include <vector>

using namespace std;

template <class charT, class traits = char_traits<charT> >
class basic_procbuf : public basic_streambuf<charT,traits>
{
  public:
    typedef charT                      char_type;
    typedef typename traits::int_type  int_type;
    typedef typename traits::pos_type  pos_type;
    typedef typename traits::off_type  off_type;
    typedef traits                     traits_type;

    basic_procbuf();
    basic_procbuf( const char* command,
		   ios_base::openmode mode );
    ~basic_procbuf();

    bool is_open()  volatile
    {  return  (this->M_fd) >= 0;  }

    basic_procbuf<charT,traits> * open( const char* command,
					ios_base::openmode mode ) volatile;
    basic_procbuf<charT,traits> * close();

  protected:
    int_type overflow( int_type c = traits_type::eof() );
    int_type underflow();
    int_type pbackfail( int_type c );
    int_type sync();

  private:
    int M_fd;
    pid_t M_pid;

    char_type *M_internal_get_buffer;
    char_type *M_internal_get_buffer_end;
    char_type *M_internal_put_buffer;
    char_type *M_internal_put_buffer_end;
    // We might read an odd # of bytes in wide mode -
    // The extra bytes go here...
    char M_get_slop[sizeof(char_type)];
    int  M_n_slop;
    bool flush();
    bool fill();

    enum {
        DEFAULT_GET_BUFFER_SIZE=512,
	DEFAULT_PUT_BUFFER_SIZE=512
    };
    
    // prohibit copying and assignment
    basic_procbuf( const basic_procbuf & );
    basic_procbuf & operator= ( const basic_procbuf & );
};

typedef basic_procbuf<char> procbuf;

#endif

#endif
