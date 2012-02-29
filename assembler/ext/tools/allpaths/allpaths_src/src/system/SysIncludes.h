///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef SYS_INCLUDES_H
#define SYS_INCLUDES_H

// -----------------------------------------------------------------------------

// This file is for common system includes, not includes of STL files
// (in general).  The main purpose of having this file is to centralize 
// platform-specific content.  There is no need to put all system includes here.

// -----------------------------------------------------------------------------

#include <string.h>
#include <fcntl.h>
#include <stdlib.h>
#include <limits>

// -----------------------------------------------------------------------------

// Note that the order of the following two includes is required under darwin.

#include <sys/time.h>
#include <sys/resource.h>

// -----------------------------------------------------------------------------

// The following includes are here because we need them to come before the
// the defines at the beginning of System.h.

#include <sys/stat.h>
#include <unistd.h>
#include <vector>

// -----------------------------------------------------------------------------

// We supply a copy of the standard procbuf.h (not supplied under all platforms).

#include "system/ProcBuf.h"

// -----------------------------------------------------------------------------

#endif
