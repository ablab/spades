/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
/*
 * \file SymLink.h
 * \author tsharpe
 * \date Feb 17, 2009
 *
 * \brief Class for reading and creating symbolic links.
 * It's a type of File.
 */
#ifndef SYSTEM_FILE_SYMLINK_H_
#define SYSTEM_FILE_SYMLINK_H_

#include "system/file/File.h"

/// Class for reading and creating symbolic links.
class SymLink : public File
{
public:
    /// Default constructor refers to empty nonsense path.
    SymLink() {}
    SymLink( char const* path ) : File(path) {}

    /// Construct from something with a c_str() member (like string or String)
    template <class C>
    explicit SymLink( C const& path,
                            char const* (C::*)() const = &C::c_str )
    : File(path.c_str()) {}

    // compiler-supplied copying and destructor are OK

    bool isValid() const { return isLink() && stat()!=0; }

    /// Read the target of the symlink.
    File target() const;

    /// Make this symlink point to the specified target.
    /// It's a fatal error if it can't be done for some reason.
    /// If force==true, it first removes any existing symlink.
    void setTarget( File const& target, bool force=false ) const;
};

#endif /* SYSTEM_FILE_SYMLINK_H_ */
