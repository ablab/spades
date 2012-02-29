/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
/**
 * \file File.h
 * \author aaron
 * \date Wednesday, December 03, 2008
 *
 * \brief A file, as described by a path.
 */
#ifndef SYSTEM_FILE_FILE_H_
#define SYSTEM_FILE_FILE_H_

#include "system/Assert.h"
#include <cstddef>
#include <ostream>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>

class Directory;
class SymLink;

/**
 * \class File
 * \brief A file, as described by a path.
 *
 * A stat on the file is done on-demand, only when this information is needed.
 */
class File
{
public:
    // make it look container-ish
    typedef std::string::value_type value_type;
    typedef std::string::const_iterator const_iterator;
    typedef std::string::const_reverse_iterator const_reverse_iterator;

    /// Default constructor creates an empty path.
    File() : mType(UNKNOWN), mpStat(0) {}

    /// Construct from a path specified as a c-style string.
    File( char const* path ) : mPath(path), mType(UNKNOWN), mpStat(0)
    { patchPath(); }

    /// Construct from something with a c_str() member (like string or String)
    template <class C>
    explicit File( C const& path, char const*(C::*)() const = &C::c_str )
    : mPath(path.c_str()), mType(UNKNOWN), mpStat(0)
    { patchPath(); }

    /// copy constructor
    File( File const& file ) : mType(UNKNOWN), mpStat(0) { *this = file; }

    /// destructor
    ~File()
    { delete mpStat; }

    /// copy-by-assignment operator
    File& operator=( File const& that )
    { if ( this != &that )
      { mPath = that.mPath; mType = that.mType;
        if ( !that.mpStat ) { delete mpStat; mpStat = 0; }
        else if ( !mpStat ) mpStat = new struct stat(*that.mpStat);
        else *mpStat = *that.mpStat; }
      return *this; }

    /**
     ** Standard container iterators
     **/
    const_iterator begin() const { return mPath.begin(); }
    const_iterator cbegin() const { return mPath.begin(); }
    const_iterator end() const { return mPath.end(); }
    const_iterator cend() const { return mPath.end(); }
    const_reverse_iterator rbegin() const { return mPath.rbegin(); }
    const_reverse_iterator crbegin() const { return mPath.rbegin(); }
    const_reverse_iterator rend() const { return mPath.rend(); }
    const_reverse_iterator crend() const { return mPath.rend(); }

    /// complete path as a C-string
    char const* c_str() const { return mPath.c_str(); }

    /**
     ** Operations that have to do with the File's status
     **/

    /// Test to see if the file exists (i.e., it can be stat'd).
    /// Call with freshen==true to reexecute the call to stat rather than
    /// relying on previously-retrieved values.
    bool exists( bool freshen = false ) const
    { if ( mType==UNKNOWN || freshen ) setStat();
      return mpStat; }

    /// Evaluate the file's status.
    /// Returns null if the file doesn't exist.
    /// Call with freshen==true to reexecute the call to stat rather than
    /// relying on previously-retrieved values.
    struct stat const* stat( bool freshen = false ) const
    { if ( mType==UNKNOWN || freshen ) setStat();
      return mpStat; }

    /// File types enumerated.
    enum FILETYPE
    {
        NOT_A_FILE = 0,
        UNKNOWN = 1,
        FIFO = S_IFIFO,
        CHARACTER_DEVICE = S_IFCHR,
        DIRECTORY = S_IFDIR,
        BLOCK_DEVICE = S_IFBLK,
        REGULAR = S_IFREG,
        SYM_LINK = S_IFLNK,
        SOCKET = S_IFSOCK
    };

    /// Return the filetype.
    /// Returns NOT_A_FILE if stat fails.  No freshen option, because files
    /// don't change type very often.  (You can always overwrite with a copy.)
    FILETYPE type() const
    { if ( mType == UNKNOWN ) setStat();
      return mType; }

    /// Is it a directory?
    bool isDir() const { return type() == DIRECTORY; }

    /// Get a directory object with the same path.
    Directory asDir() const;

    /// Is it a symbolic link?
    bool isLink() const { return type() == SYM_LINK; }

    /// Get a symbolic link object with the same path.
    SymLink asLink() const;

    /// Return the size of the file.
    /// Returns -1 if stat fails.
    off_t filesize( bool freshen=false ) const
    { if ( mType==UNKNOWN || freshen ) setStat();
      return mpStat ? mpStat->st_size : -1L; }

    /// Last access time.
    /// Returns 0 if stat fails.
    time_t accessTime( bool freshen=false ) const
    { if ( mType==UNKNOWN || freshen ) setStat();
      return mpStat ? mpStat->st_atime : 0; }

    /// Last modification time.
    /// Returns 0 if stat fails.
    time_t modifyTime( bool freshen=false ) const
    { if ( mType==UNKNOWN || freshen ) setStat();
      return mpStat ? mpStat->st_mtime : 0; }

    /// Last status change time.
    time_t createTime( bool freshen=false ) const
    { if ( mType==UNKNOWN || freshen ) setStat();
      return mpStat ? mpStat->st_ctime : 0; }

    /// Checks to see if two (possibly different) paths refer to same the file.
    bool isSameFile( File const& file ) const
    { struct stat const* pStat1 = stat();
      struct stat const* pStat2 = file.stat();
      return pStat1 && pStat2 &&
              pStat1->st_ino == pStat2->st_ino &&
              pStat1->st_dev == pStat2->st_dev; }

    /**
     ** Operations on the File's path
     **/

    bool isAbsolute() const
    { return mPath.size() && mPath[0]=='/'; }

    bool isExplicitlyRelative() const
    { return mPath.size()>1 && mPath[0]=='.' && mPath[1]=='/'; }

    /// Return the file's path as a string.
    std::string const& toString() const { return mPath; }

    /// Return the filename portion of the path.
    /// E.g., "/home/user/boat/file.txt" would return "file.txt".
    /// Anomaly:  The filename of the root directory is empty.
    std::string filename() const
    { return mPath.substr(filenameOffset()); }

    /// Return the directory portion of the path.
    /// E.g., "/home/user/boat/file.txt" would return "/home/user/boat/",
    /// and barefile.txt would return "./".
    /// Anomaly:  The dirname of the root directory is empty.
    std::string dirname() const
    { size_t pos = filenameOffset();
      return pos ? mPath.substr(0,pos-1) : CWD; }

    /// The directory that contains this file.
    /// I.e., the directory described by dirname().
    Directory directory() const;

    /// Return the final extension (or an empty string, if there is none).
    /// E.g., "/home/user/boat/file.txt" would return ".txt".
    std::string extension() const
    { std::string result;
      size_t start = mPath.find_last_of('.');
      if ( start != std::string::npos && start > filenameOffset() )
          result = mPath.substr(start);
      return result; }

    /// Return a new File with its final extension (if any) removed.
    /// E.g., "/home/user/boat/file.txt" would return "/home/user/boat/file".
    File removeExtension() const
    { size_t end = mPath.find_last_of('.');
      return end != std::string::npos && end > filenameOffset() ?
                      File(mPath.substr(0,end)) : *this; }

    /// Return a new File with a path obtained by appending a new extension.
    /// If the new extension doesn't start with a '.', one is added, unless
    /// the new extension is empty, in which case the path is unaltered.
    File addExtension( std::string newExtension ) const
    { AssertEq(newExtension.find_first_of('/'),std::string::npos);
      size_t xSize = newExtension.size();
      std::string newPath; newPath.reserve(mPath.size()+1+xSize);
      newPath = mPath;
      if ( xSize )
      { if ( newExtension[0] != '.' ) newPath += '.';
        newPath += newExtension; }
      return File(newPath); }

    /// Return a new File with the final extension changed.
    /// Equivalent to removeExtension + addExtension.
    /// E.g., "/user/boat/file.txt" would return "/user/boat/file.newExtension".
    File changeExtension( std::string newExtension ) const
    { AssertEq(newExtension.find_first_of('/'),std::string::npos);
      size_t end = mPath.find_last_of('.');
      if ( end == std::string::npos || end <= filenameOffset() )
          end = mPath.size();
      size_t xSize = newExtension.size();
      std::string newPath; newPath.reserve(end+1+xSize);
      newPath = mPath.substr(0,end);
      if ( xSize )
      { if ( newExtension[0] != '.' ) newPath += '.';
        newPath += newExtension; }
      return File(newPath); }

    /**
     ** Operations on the File itself.
     **/

    /// Remove the file from its directory.
    /// It's a fatal error if it doesn't happen (so consider calling the
    /// exists() method first).
    void remove() const;

    /// Give the file a new name.  If you're crossing filesystems, this might
    /// actually cause the file to be copied (which might be expensive).
    void rename( File const& file ) const;

    /**
     *** Some friends.
     **/

    /// Tests whether the paths are identical.
    friend bool operator==( File const& file1, File const& file2 )
    { return file1.toString() == file2.toString(); }

    /// Tests whether the paths are different.
    friend bool operator!=( File const& file1, File const& file2 )
    { return !(file1 == file2); }

    /// Compares paths lexicographically.
    friend bool operator<( File const& file1, File const& file2 )
    { return file1.toString() < file2.toString(); }

    friend std::ostream& operator<<( std::ostream& s, File const& file )
    { s << file.toString(); return s; }

    static std::string CWD; // i.e., "."
    static std::string ROOT;

private:
    friend class Directory;
    friend class SymLink;

    void patchPath();
    size_t filenameOffset() const
    { return mPath.find_last_of('/')+1; } // trick: +1 turns npos into 0

    void clearStat() const { mType = UNKNOWN; delete mpStat; mpStat = 0; }
    void setStat() const;

    std::string mPath;
    mutable FILETYPE mType;
    mutable struct stat* mpStat;
};

#endif /* SYSTEM_FILE_FILE_H_ */
