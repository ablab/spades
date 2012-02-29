/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef PARALLEL_H
#define PARALLEL_H

#include "CoreTools.h"

// Execute commands in parallel, running at most max_threads of them
// at once. If not null, use the pointer stream pOut to show progress.

void RunCommandsInParallel( const vec<String>& commands,
			    const int max_threads,
			    ostream *pOut = NULL );

#endif
