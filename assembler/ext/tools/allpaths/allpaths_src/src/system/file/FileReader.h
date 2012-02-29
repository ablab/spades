///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file FileReader.h
 * \author tsharpe
 * \date Jan 12, 2012
 *
 * \brief Handle reading from a file descriptor.
 */
#ifndef SYSTEM_FILE_FILEREADER_H_
#define SYSTEM_FILE_FILEREADER_H_

#include "system/file/File.h"
#include <string>
#include <cstddef>

/// Handle reading from a file descriptor.
class FileReader
{
public:
    explicit FileReader( char const* path )
    : mPath(path) { doOpen(); }

    /// Construct from something with a c_str() member (like string or String)
    template <class C>
    explicit FileReader( C const& path,
                             char const* (C::*)() const = &C::c_str )
    : mPath(path.c_str()) { doOpen(); }

    FileReader( int fd, char const* pseudoFilename )
    : mFD(fd), mMyFD(false), mPath(pseudoFilename) {}

    ~FileReader() { close(); }

    std::string const& getFilename() const { return mPath; }

    /// Reads as much as is immediately available.  Returns 0 at EOF.
    size_t readOnce( void* buf, size_t len ) const;

    /// Reads bytes from the file into the buffer, returning the number read.
    /// This will be < len only if we hit EOF, all other errors are fatal.
    size_t readSome( void* buf, size_t len ) const;

    /// Reads exactly len bytes from the file into the buffer.
    /// If we hit EOF before getting that many, it's a fatal error.
    FileReader const& read( void* buf, size_t len ) const;

    /// SEEK_SET to this offset
    FileReader const& seek( size_t off ) const;

    /// SEEK_CUR with this offset.
    FileReader const& seekRel( long off ) const;

    /// SEEK_END with this offset.
    FileReader const& seekEnd( long off ) const;

    struct stat getStat() const;

    /// Return the file's size.
    size_t getSize() const { return getStat().st_size; }

    /// Memory-map the file.
    void* map( size_t offset, size_t len, bool readOnly=false );

    bool isOpen() const { return mFD != -1; }
    void close()
    { if ( !mMyFD ) mFD = -1;
      else if ( mFD != -1 ) doClose(); }

private:
    FileReader( FileReader const& ); // not implemented -- no copying
    FileReader& operator=( FileReader const& ); // not implemented -- no copying
    // note that we could implement copying by using ::dup(), if necessary

    void doOpen();
    void doClose();
    void wait() const;

    int mFD;
    bool mMyFD;
    std::string mPath;

    static size_t const MAX_IO_LEN = (1ul<<31) - (1ul<<12);

    friend class FileWriter;
};

#endif /* SYSTEM_FILE_FILEREADER_H_ */
