/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/**
    Header file: Utils.h

    Miscellaneous programming utilities.

    @file
 */

#ifndef __INCLUDE_system_Utils_h
#define __INCLUDE_system_Utils_h

#include <stdlib.h>
#include "system/Types.h"
#include "String.h"

/**
    Class: NonCopyable

    Derive your class from this to prevent it from being copied.
    Useful when the default copy constructor is incorrect and
    you don't want to write a correct one.
 */
class NonCopyable {
protected:
  NonCopyable() {}

private:
  NonCopyable( NonCopyable const& );
  NonCopyable& operator=( NonCopyable const& );
};

// Function: FromEnv
// Return an integer value from the environment,
// returning the specified default if the environment var
// is not defined.
inline int FromEnv( const char *s, int dflt ) {
  return getenv( s ) ? atoi( getenv( s ) ) : dflt;
}


/**
   Function: IntervalsOverlap

   Tests whether two closed discrete intervals [x_lower,x_upper] and
   [y_lower,y_upper] share at least one point.
*/
template <class T> inline
Bool IntervalsOverlap( T x_lower, T x_upper, T y_lower, T y_upper ) {
  return x_lower <= y_lower && y_lower <= x_upper ||
         y_lower <= x_lower && x_lower <= y_upper;
}

void StrongWarning( String msg );


#endif
// #ifndef __INCLUDE_system_Utils_h
