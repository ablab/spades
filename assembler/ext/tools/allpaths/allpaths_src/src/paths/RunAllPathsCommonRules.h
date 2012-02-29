///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/**
   Contains generalized rules used in RunAllPaths and RunAllPathsLG.
*/

#ifndef RUN_ALLPATHS_COMMON_RULES_H
#define RUN_ALLPATHS_COMMON_RULES_H

#include "String.h"
#include "system/MiscUtil.h"


void AddRule_MakeLookupTable( MakeMgr& make,
			      String refFile, String refDir,
			      String comment = "");


void AddRule_SQLAlign( MakeMgr& make,
		       String refFile, String queryFile,
		       String refDir, String queryDir, String outDir,
		       int maxErrors, Bool fwrc,
		       int maxIndelLen, int errDiff,
		       String comment = "");

void AddRule_LookupAndSQLAlign( MakeMgr& make,
				String refFile, String queryFile,
				String refDir, String queryDir, String outDir,
				int maxErrors, Bool fwrc,
				int maxIndelLen, int errDiff,
				String comment = "");


void AddRule_CommonPather(MakeMgr& make,
			  int nThreads,
			  int K, String outDir,
			  String dir1, String file1,
			  String dir2 = "", String file2 = "",
			  String dir3 = "", String file3 = "");


void AddRule_CommonPather(MakeMgr& make,
			  int K, String outDir,
			  String dir1, String file1,
			  String dir2 = "", String file2 = "",
			  String dir3 = "", String file3 = "");


void AddRule_MergeReadSets(MakeMgr& make,
			   bool quals, bool pairs, bool paths, bool dist,
			   bool tracker, bool repath, int K, int nThreads,
                           const String& dir, const String& headOut,
			   const String& headIn1, const String& headIn2, 
			   const String& headIn3 = "", const String& headIn4 = "",
			   const String& headIn5 = "", const String& headIn6 = "" );


#endif
