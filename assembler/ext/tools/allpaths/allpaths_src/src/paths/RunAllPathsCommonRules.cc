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

#include "String.h"
#include "paths/RunAllPathsCommonRules.h"
#include "system/MiscUtil.h"


void AddRule_MakeLookupTable( MakeMgr& make,
			      String refFile, String refDir,
			      String comment)
{
  String refHead = refFile.SafeBefore(".fastb");
  String lookupFile = refHead + ".lookup";

// MakeDepend: dependency MakeLookupTable
  make.AddRule( FilesIn( refDir, lookupFile),
		FilesIn( refDir, refFile ),
		"MakeLookupTable"
		" SOURCE=" + refDir + "/" + refFile +
		" OUT_HEAD=" + refDir + "/" + refHead +
		" LOOKUP_ONLY=True",
		comment);
}


void AddRule_SQLAlign( MakeMgr& make,
		       String refFile, String queryFile,
		       String refDir, String queryDir, String outDir,
		       int maxErrors, Bool fwrc,
		       int maxIndelLen, int errDiff,
		       String comment)
{
  String queryHead = queryFile.SafeBefore(".fastb");
  String refHead = refFile.SafeBefore(".fastb");
  String lookupFile = refHead + ".lookup";

  if (queryDir == "")
    queryDir = refDir;
  if (outDir == "")
    outDir = queryDir;
 
  String alignsHead = queryHead + ".aligned." + refHead;
  String alignsFile = alignsHead + ".sql.qltout";

// MakeDepend: dependency ShortQueryLookup
  make.AddRule( FilesIn( outDir, alignsFile ),
		JoinVecs( FilesIn( refDir, lookupFile),
			  FilesIn( queryDir, queryFile)),
		"ShortQueryLookup "
		" SEQS=" + queryDir + "/" + queryFile +
		" LOOKUP_TABLE=" + refDir + "/" + lookupFile +
		" OUT_PREFIX=" + outDir + "/" + alignsHead +
		" OUT_SUFFIX=.sql.qltout"
		" FWRC=" +ToStringBool(fwrc) +
		" MAX_INDEL_LEN=" + ToString(maxIndelLen) +
		" MAX_ERRS=" + ToString(maxErrors) +
		" ERR_DIFF=" + ToString(errDiff) +
		" CHUNK=500000 READABLE=False",
		
		comment );
}

void AddRule_LookupAndSQLAlign( MakeMgr& make,
				String refFile, String queryFile,
				String refDir, String queryDir, String outDir,
				int maxErrors, Bool fwrc,
				int maxIndelLen, int errDiff,
				String comment)
{

  // Make Lookup Table

  AddRule_MakeLookupTable( make, refFile, refDir,
			   "Build lookup table.");
    
  // Perform Alignment

  AddRule_SQLAlign( make, refFile, queryFile, refDir, queryDir, outDir,
		    maxErrors, fwrc, maxIndelLen, errDiff, comment);


}

void AddRule_CommonPather(MakeMgr& make,
			  int K, String outDir,
			  String dir1, String file1,
			  String dir2, String file2,
			  String dir3, String file3) {

  AddRule_CommonPather(make, 1, K, outDir, dir1, file1, dir2, file2, dir3, file3);
}

void AddRule_CommonPather(MakeMgr& make,
			  int nThreads,
			  int K, String outDir,
			  String dir1, String file1,
			  String dir2, String file2,
			  String dir3, String file3) {


  vec< filename_t > sourceFiles, targetFiles;
  String readsIn  = "\"{";
  String pathsOut = "\"{";
  String KS = "k" + ToString(K);

  // Add file1 to list to path
  String file1Head = file1.SafeBefore(".fastb");
  readsIn  += dir1 +   "/" + file1;
  pathsOut += outDir + "/" + file1Head + ".paths";
  sourceFiles.push_back(dir1 + "/" + file1);
  targetFiles.push_back(outDir + "/" + file1Head + ".paths." + KS);
  String comment = file1;

  // Add file2 to list to path
  if (file2 != "") {
    String file2Head = file2.SafeBefore(".fastb");
    readsIn  += "," + dir2 +   "/" + file2;
    pathsOut += "," + outDir + "/" + file2Head + ".paths";
    sourceFiles.push_back(dir2 + "/" + file2);
    targetFiles.push_back(outDir + "/" + file2Head + ".paths." + KS);
    comment += ", " + file2;
  }

  // Add file3 to list to path
  if (file3 != "") {
    String file3Head = file3.SafeBefore(".fastb");
    readsIn  += "," + dir3 +   "/" + file3;
    pathsOut += "," + outDir + "/" + file3Head + ".paths";
    sourceFiles.push_back(dir3 + "/" + file3);
    targetFiles.push_back(outDir + "/" + file3Head + ".paths." + KS);
    comment += ", " + file3;
  }

  readsIn  += "}\"";
  pathsOut += "}\"";

// MakeDepend: dependency CommonPather
  make.AddRule( targetFiles, sourceFiles,
		"CommonPather K=" + ToString(K) +
		" NUM_THREADS=" + ToString(nThreads) + 
		" READS_IN=" + readsIn +
		" PATHS_OUT=" + pathsOut,
		
		"Pathing reads: " + comment );
}


void AddRule_MergeReadSets(MakeMgr& make,
			   bool quals, bool pairs, bool paths, bool dist,
			   bool tracker, bool repath, int K, int nThreads,
                           const String& dir, const String& headOut,
			   const String& headIn1, const String& headIn2, 
			   const String& headIn3, const String& headIn4,
			   const String& headIn5, const String& headIn6 ) {

  // Decide which file types are required as input and output

  String suffixIn  = "{fastb";
  String suffixOut = "{fastb";

  if (quals) { suffixIn += ",qualb";   suffixOut += ",qualb";  } 
  if (pairs) { suffixIn += ",pairs";  suffixOut += ",pairs"; }
  if (dist)  { suffixOut += ",distribs"; }
  if (paths && !repath) suffixIn  += ",paths.kN";
  if (paths ||  repath) suffixOut += ",paths.kN";
  if (tracker) suffixOut += ",readtrack";
  
  suffixIn  += "}";
  suffixOut += "}";

  // Create a list of all source files required

  vec<String> headsIn;
  
  headsIn.push_back(headIn1.SafeBefore(".fastb"));
  headsIn.push_back(headIn2.SafeBefore(".fastb"));
  if (headIn3 != "" ) headsIn.push_back(headIn3.SafeBefore(".fastb"));
  if (headIn4 != "" ) headsIn.push_back(headIn4.SafeBefore(".fastb"));
  if (headIn5 != "" ) headsIn.push_back(headIn5.SafeBefore(".fastb"));
  if (headIn6 != "" ) headsIn.push_back(headIn6.SafeBefore(".fastb"));


  String headsInArg  = "{" + headsIn[0];
  vec< String > sourceFiles(FilesIn(dir, headsIn[0] + "." + suffixIn));

  for (size_t i = 1; i < headsIn.size(); ++i) {
    headsInArg  += "," + headsIn[i];
    sourceFiles.append(FilesIn(dir, headsIn[i] + "." + suffixIn));
  }  
					
  headsInArg  += "}";

  // Create a list of all target files to be created

  String headOutArg = headOut.SafeBefore(".fastb");
  vec< String > targetFiles = FilesIn(dir, headOutArg + "." + suffixOut);


// MakeDepend: dependency MergeReadSets
  make.AddRule( targetFiles, sourceFiles,
		"MergeReadSets "
		" DIR=" + dir +
		" READS_IN=" + "\"" + headsInArg + "\""
		" READS_OUT=" + headOutArg +
		" K=" + ToString(K) +
		" REPATH=" + ToStringBool(repath) +
		" NUM_THREADS=" + ToString(nThreads) +
		" MERGE_QUALS=" + ToStringBool(quals) +
		" MERGE_PAIRS=" + ToStringBool(pairs) +
		" MERGE_PATHS=" + ToStringBool(paths) +
		" MERGE_DISTRIBS=" + ToStringBool(dist) + 
		" FORCE_MERGE_DISTRIBS=" + ToStringBool(dist) +
		" TRACK_READS=" + ToStringBool(tracker),
		
		"Merging reads: " + headsInArg);

}
