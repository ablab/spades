///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef MAIN_TOOLS_H
#define MAIN_TOOLS_H

// MainTools.h extends CoreTools.h with includes needed by almost all
// main program modules.  It should not be included into .cc files
// that do not define a main(), to keep compile times down.

#include "CoreTools.h"
#include "system/ParsedArgs.h"
#include "system/RunTime.h"
#include "system/SysConf.h"

#endif
