///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// This file defines class fast_ifstream and fast_pipe_ifstream, which has
// essentially identical member functions.

#ifndef FAST_IFSTREAM
#define FAST_IFSTREAM

#include "CoreTools.h"
#include "system/file/FileReader.h"

class fast_ifstream {

     public:

    fast_ifstream( const String& filename )
    : fr_(filename.c_str()), buf_ptr_(0), buf_top_(0), fail_(False) {}
    fast_ifstream( char const* filename )
    : fr_(filename), buf_ptr_(0), buf_top_(0), fail_(False) {}

    // compiler-supplied destructor is OK

     /// Rewinds the currently open file back to start position;
     /// has no effect if the stream is not associated with a file
     /// (i.e. default constructor was used and no Open() was issued yet)
     void Rewind() ;

     // get(c) reads a char and puts it in c

     void get( char& c, Bool advance = True );

     ///Read additional characters beyond the internal buffer.
     void ReadSlow(char * & start, char * & end);

     // peek(c) peeks at a char and puts it in c
     // If there is nothing left in the file, then fail is set.

     void peek( char& c );

     Bool fail( ) const
     {    return fail_;    }

     // Read a line from "in" and put it in "s".

     friend void getline( fast_ifstream& in, String& s );

     // Read a line form "in".  If it begins with "begin", put it in "s", return 
     // True.  Otherwise, return False.
     //
     // This is slightly faster than calling getline.

     friend Bool getline_if_match( fast_ifstream& in, String& s, 
          const String& begin );

     // Read until a string "tail" is encountered.  Set "s" to everything read,
     // including the tail.  If the tail is not found, return False;

     friend Bool get_to( fast_ifstream& in, String& s, const String& tail );

     friend void ReadParagraphs( fast_ifstream& in, vec<String>& paragraphs );

     private:

     fast_ifstream( fast_ifstream const& ); // unimplemented -- no copying
     fast_ifstream& operator=( fast_ifstream const& ); // unimplemented -- no copying

     static int const BUFFER_SIZE = 8192;

     FileReader fr_;
     char buf_[BUFFER_SIZE + 1];
     int buf_ptr_, buf_top_;
     Bool fail_;

     void fill_buffer( );

};

void getline( fast_ifstream& in, String& s );
Bool getline_if_match( fast_ifstream& in, String& s, const String& begin );

class fast_pipe_ifstream {

     public:

     fast_pipe_ifstream( String command );

     ~fast_pipe_ifstream( );

     void get( char& c, Bool advance = True );

     void peek( char& c );

     Bool fail( ) const
     {    return fail_;    }

     friend void getline( fast_pipe_ifstream& in, String& s );

     friend Bool getline_if_match( fast_pipe_ifstream& in, String& s, 
          const String& begin );

     friend Bool get_to( fast_pipe_ifstream& in, String& s, const String& tail );

     friend void ReadParagraphs( fast_pipe_ifstream& in, vec<String>& paragraphs );

     private:

     fast_pipe_ifstream( fast_pipe_ifstream const& ); // unimplemented -- no copying
     fast_pipe_ifstream& operator=( fast_pipe_ifstream const& ); // unimplemented -- no copying

     static int const BUFFER_SIZE = 8192;

     String command_;
     FILE* file_;
     char buf_[BUFFER_SIZE + 1];
     int buf_ptr_, buf_top_;
     Bool fail_;

     void fill_buffer( );

};

#endif
