///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
#include <cstdio>
#include <iostream>
#include <cerrno>
#include <unistd.h>
#include <sys/wait.h>
#include "system/Assert.h"
#include "system/ErrNo.h"
#include "system/ProcBuf.h"
#include "system/SysIncludes.h"

#if ! defined( __GNUC__ ) || __GNUC__ > 2 
#include <sstream>

template<class charT, class traits>
basic_procbuf<charT,traits>::basic_procbuf( )
    : M_fd( -1 ),
      M_internal_get_buffer(NULL),
      M_internal_get_buffer_end(NULL),
      M_internal_put_buffer(NULL),
      M_internal_put_buffer_end(NULL),
      M_n_slop(0)
{
}


template<class charT, class traits>
basic_procbuf<charT,traits>::basic_procbuf( const char *command,
					    ios_base::openmode mode )
    : M_fd( -1 ),
      M_internal_get_buffer(NULL),
      M_internal_get_buffer_end(NULL),
      M_internal_put_buffer(NULL),
      M_internal_put_buffer_end(NULL),
      M_n_slop(0)

{
    this->open( command, mode );
}


template<class charT, class traits>
basic_procbuf<charT,traits>::~basic_procbuf( )
{
    this->close();
    if (M_internal_get_buffer)
      delete [] M_internal_get_buffer;
    if (M_internal_put_buffer)
      delete [] M_internal_put_buffer;
}


// template<class charT, class traits>
// bool
// basic_procbuf<charT,traits>::is_open( )
// {
//     return M_fd >= 0;
// }


template<class charT, class traits>
basic_procbuf<charT,traits> *
basic_procbuf<charT,traits>::open( const char *command,
				   ios_base::openmode mode ) volatile
{
    int parent_in, parent_out, child_in, child_out;
    int pipe_fds[2];
    pid_t child_pid;
    
    if( this->is_open() )
    {
        printf( "Attempted open() on already open procbuf.\n" );
        TracebackThisProcess();
        return 0;
    }
    
    if( pipe (pipe_fds) < 0 )
    {
        ErrNo err;
        std::cout << "pipe() failed" << err << std::endl;
        TracebackThisProcess();
        return 0;
    }

    parent_in  = child_in  = pipe_fds[0];
    parent_out = child_out = pipe_fds[1];

    if (mode & ios_base::in) 
        mode = ios_base::in;
    else
        mode = ios_base::out;

    this->M_pid = child_pid = vfork ();
    if (child_pid == 0)
    {
        if ( mode == ios_base::in )
        {
          ::close( child_in );
          dup2( child_out, STDOUT_FILENO );
          //::close( child_out );
        }
        else
        {
          ::close( child_out );
          dup2( child_in, STDIN_FILENO );
          //::close( child_in );
        }

        char* shell = getenv("SHELL");
        if ( shell )
            execl ("/bin/sh", shell, "-c", command, (char *) 0);

        else 
            execl ("/bin/sh", "sh", "-c", command, (char *) 0);

        _exit (-1);
    }

    else if (child_pid < 0) 
    {
        ErrNo err;
        std::cout << "fork() failed" << err << std::endl;
        ::close( parent_in );
        ::close( parent_out );
        TracebackThisProcess();
        return 0;
    }

    if ( mode == ios_base::in )
    {
      ::close( parent_out );
      this->M_fd = parent_in;
    }
    else
    {
      ::close( parent_in );
      this->M_fd = parent_out;
    }

    return const_cast< basic_procbuf<charT,traits> * >(this);
}


template<class charT, class traits>
basic_procbuf<charT,traits> *
basic_procbuf<charT,traits>::close()
{
    int wstatus;
    basic_procbuf **ptr;
    pid_t wait_pid;
    int status = -1;

    flush();
    if ( ::close(this->M_fd) < 0 ) {
        this->M_fd=-1;
        return this;
    }
    this->M_fd=-1;

    /* While there is no child complete, or an interrupt, wait: */
    do 
    {
        wait_pid = waitpid (this->M_pid, &wstatus, 0);
    } 
    while (wait_pid == -1 && errno == EINTR);

    if (wait_pid == -1)
        return this;

    return this;
}


template<class charT, class traits>
typename traits::int_type
basic_procbuf<charT, traits>::sync ( )
{
    if (this->pptr() != this->pbase() && ! flush())
	return -1;
    return 0;
}


template<class charT, class traits>
typename traits::int_type
basic_procbuf<charT,traits>::overflow( int_type c )
{
    if ( ! traits_type::eq_int_type( c, traits_type::eof() ) )
    {
	if (this->pbase() == this->epptr()) {
	    // Empty buffer - for the first pass, we just use a fixed size
	    Assert(M_internal_put_buffer == NULL);
	    M_internal_put_buffer = new char_type [DEFAULT_PUT_BUFFER_SIZE];
	    M_internal_put_buffer_end = M_internal_put_buffer+DEFAULT_PUT_BUFFER_SIZE;
	    setp(M_internal_put_buffer,
		 M_internal_put_buffer_end);
	} else if (!flush())
	    return traits_type::eof();
	if (!traits_type::eq_int_type(c, traits_type::eof()))
	    return sputc(c);
	else
	    return traits_type::not_eof(c);
    }
    return traits_type::eof();
}


template<class charT, class traits>
typename traits::int_type
basic_procbuf<charT,traits>::underflow( )
{
    if (this->gptr() < this->egptr() || fill())
        return traits_type::to_int_type( *(this->gptr()) );
    else
        return traits_type::eof();
}


template<class charT, class traits>
typename traits::int_type
basic_procbuf<charT,traits>::pbackfail( int_type c )
{
    // If the character is EOF, then the caller isn't telling us what it is.
    // If we lost the buffer... return failure
    bool bKnowsChar = (!traits_type::eq_int_type(c, traits_type::eof()));

    // Check for room at the beginning of the buffer
    if (this->gptr() != this->eback()) {
	this->gbump(-1);
        if (bKnowsChar)
	    *(this->gptr()) = traits_type::to_char_type(c);
        return traits_type::not_eof(c);
    } else if (bKnowsChar && 
	       this->eback() == M_internal_get_buffer &&
	       this->egptr() < M_internal_get_buffer_end) {
	// Shift to make space, put character at current loc
	// (which is start of buffer)
	std::copy_backward(this->eback(), this->egptr(), this->egptr()+1);
	*(this->gptr()) = traits_type::to_char_type(c);
	setg(this->eback(), this->gptr(), this->egptr()+1);
        return traits_type::not_eof(c);
    } else
	return traits_type::eof();
}


extern void logSomething(const char *pszMessage);


template<class charT, class traits>
bool
basic_procbuf<charT,traits>::flush( )
{
    // Write the buffer to a pipe.

    int expected=(this->pptr() - this->pbase()) * sizeof(char_type);
    if (expected == 0)
      return true;
    int written=write(this->M_fd,
		      this->pbase(),
		      (this->pptr() - this->pbase()) * sizeof(char_type));

    if (written != expected) {
	// Did not make any room...probably error or closed pipe.
	// TBD: We could allow partial writes - but what happens if we're
	// doing wchars and we write an odd # of bytes?
	return false;
    } else {
	setp(this->pbase(), this->epptr());
	return true;
    }
}
	
template<class charT, class traits>
bool
basic_procbuf<charT,traits>::fill( )
{
    if (this->eback() == this->egptr()) {
	// The get buffer hasn't been built yet.
	Assert(M_internal_get_buffer == NULL);
	M_internal_get_buffer = new char_type[DEFAULT_GET_BUFFER_SIZE];
	M_internal_get_buffer_end = M_internal_get_buffer + DEFAULT_GET_BUFFER_SIZE;
	setg(M_internal_get_buffer, M_internal_get_buffer_end, M_internal_get_buffer_end);
    }
    // The "get" pointer should be at the end of the buffer - that's
    // why we need to fill it.
    Assert(this->gptr() == this->egptr());
    //
    // Put the "slop" bytes that we didn't use last time
    // into the start of the buffer.
    //
    char_type *p1=M_internal_get_buffer;
    char *ptr=reinterpret_cast<char *>(p1);
    Assert(size_t(M_n_slop) < sizeof(char_type));
    char *readptr= std::copy(M_get_slop, M_get_slop+M_n_slop, ptr);
    //
    // Read some - loop until we get at least 1 char_type
    //
    do {
	int numread=read( this->M_fd, 
			  readptr, 
			  DEFAULT_GET_BUFFER_SIZE * sizeof(char_type)-(readptr-ptr));
	if (numread <= 0) {
	    return false;
	}
	readptr+=numread;
    } while(size_t(readptr-ptr) < sizeof(char_type));

    //
    // compute and copy the new slop
    //
    M_n_slop = (readptr-ptr) % sizeof(char_type);
    int numusable=(readptr-ptr)-M_n_slop;
    std::copy(ptr+numusable, readptr, M_get_slop);
    //
    // reset the get pointers
    //
    setg(M_internal_get_buffer, 
	 M_internal_get_buffer, 
	 M_internal_get_buffer+numusable/sizeof(char_type));
    return true;
}


// instantiate specialization of basic_procbuf<char>

#define pbspec basic_procbuf<char>

template pbspec::basic_procbuf();
template pbspec::basic_procbuf( const char *command,
                                ios_base::openmode mode );
template pbspec::~basic_procbuf();

template bool pbspec::is_open() volatile;
template pbspec * pbspec::open( const char* command, 
                                ios_base::openmode mode ) volatile;
template pbspec * pbspec::close();

template pbspec::int_type pbspec::overflow( pbspec::int_type c = pbspec::traits_type::eof() );
template pbspec::int_type pbspec::underflow();
template pbspec::int_type pbspec::pbackfail( pbspec::int_type c );
template bool pbspec::fill();
template bool pbspec::flush();

#endif

