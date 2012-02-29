///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file FileWriter.h
 * \author tsharpe
 * \date Jan 20, 2012
 *
 * \brief Handle writing to a file descriptor.
 */
#ifndef FILEWRITER_H_
#define FILEWRITER_H_

#include "system/file/FileReader.h"

class FileWriter : public FileReader
{
public:
    explicit FileWriter( char const* path )
    : FileReader(doOpen(path),path) { mMyFD = true; }

    /// Construct from something with a c_str() member (like string or String)
    template <class C>
    explicit FileWriter( C const& path,
                            char const* (C::*)() const = &C::c_str )
    : FileReader(doOpen(path.c_str()),path.c_str())
    { mMyFD = true; }

    FileWriter( int fd, char const* pseudoFilename )
    : FileReader(fd,pseudoFilename) {}

    // copying allowed if (or when) base class allows it
    // compiler-suppied destructor is OK

    FileWriter const& write( void const* buf, size_t len ) const;

private:
    static int doOpen( char const* path );
};

#endif /* FILEWRITER_H_ */
