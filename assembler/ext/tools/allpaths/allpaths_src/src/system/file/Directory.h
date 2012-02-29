/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
/**
 * \file Directory.h
 * \author aaron
 * \date Thursday, January 22, 2009
 *
 * \brief A File that contains other Files.
 */
#ifndef SYSTEM_FILE_DIRECTORY_H_
#define SYSTEM_FILE_DIRECTORY_H_

#include "system/file/File.h"
#include <sys/types.h>
#include <dirent.h>

/**
 * \class Directory
 * \brief A Directory is a type of file that allows listing of its contents.
 *
 */
class Directory : public File
{
public:
    class const_iterator : public std::iterator<std::input_iterator_tag,File>
    {
    public:
        const_iterator() : mpDIR(0), mPos(-1L) {}

        explicit const_iterator( std::string const& path );

        const_iterator( const_iterator const& that )
        : mDir(that.mDir), mpDIR(0), mPos(-1L) { *this = that; }

        ~const_iterator() { if ( mpDIR ) endStream(); }

        const_iterator& operator=( const_iterator const& that );

        File operator*() const { return mFile; }

        File* operator->() { return &mFile; }

        const_iterator& operator++() { nextEntry(); return *this; }

        const_iterator operator++(int)
        { const_iterator tmp(*this); nextEntry(); return tmp; }

        friend bool operator==( const_iterator const& itr1,
                                const_iterator const& itr2 )
        { return itr1.mPos == itr2.mPos; }

        friend bool operator!=( const_iterator const& itr1,
                                const_iterator const& itr2 )
        { return !(itr1 == itr2); }

    private:
        void checkStream();
        void endStream();
        void nextEntry();

        std::string mDir;
        DIR* mpDIR;
        off_t mPos;
        File mFile;
    };

    /// Default constructor refers to the root directory.
    Directory() : File(File::ROOT) {}
    Directory( char const* path ) : File(path) {}

    /// Construct from something with a c_str() member (like string or String)
    template <class C>
    explicit Directory( C const& path,
                            char const* (C::*)() const = &C::c_str )
    : File(path.c_str()) {}

    // compiler-supplied copying and destructor are OK

    /// Is this an existing directory?
    bool isValid() const { return isDir(); }

    /// Create the directory associated with this path, if necessary.
    /// Returns false if the directory already existed.
    /// It's a fatal error if it can't be done, e.g., if the file already
    /// exists as a regular file.
    /// If recursive==true, also create parent dirs, as necessary.
    bool create( bool recursive=false, int mode=0777 ) const;

    /// Removes the directory.
    /// Directory must be empty.
    void remove() const;

    /// Get a path for some file relative to this directory.
    File file( std::string filename ) const
    { return File(toString()+'/'+filename); }

    /// Get a path for some sub-directory relative to this directory.
    Directory subdir( std::string filename ) const
    { return Directory(toString()+'/'+filename); }

    const_iterator begin() const { return const_iterator(toString()); }
    const_iterator end() const { return const_iterator(); }
};

#endif /* SYSTEM_FILE_DIRECTORY_H_ */
