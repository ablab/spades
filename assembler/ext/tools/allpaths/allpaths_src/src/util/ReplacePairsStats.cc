///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

const char *DOC =
  "Replaces library stats with those from another pairs file.";

#include "MainTools.h"
#include "Basevector.h"
#include "PairsManager.h"
#include "lookup/LibInfo.h"
#include "Vec.h"


int main( int argc, char *argv[] )
{
  RunTime( );
  BeginCommandArguments;
  CommandDoc(DOC);
  CommandArgument_String_Doc(PAIRS_IN,
    "Pairs file with original stats");
  CommandArgument_String_OrDefault_Doc(PAIRS_OUT, "",
    "Updated pairs file (now containing replacement stats)");
  CommandArgument_String_OrDefault_Doc(STATS_IN, "",
    "Pairs file containing replacement stats.");
  EndCommandArguments;

  if (PAIRS_OUT=="") {
    PAIRS_OUT=PAIRS_IN;
    cout << "Warning - about to overwrite original pairs file." << endl;
  }

  // Strip .pairs from filenames

  PAIRS_IN = PAIRS_IN.SafeBefore(".pairs");
  PAIRS_OUT = PAIRS_OUT.SafeBefore(".pairs");
  STATS_IN = STATS_IN.SafeBefore(".pairs");
  
  // Load pairing information

  cout << Date() << ": Loading original pairing info from:" << endl;
  cout << "  " << PAIRS_IN + ".pairs" << endl;
  PairsManager pairs;
  pairs.Read(PAIRS_IN + ".pairs");


  cout << Date() << ": Loading replacement pairing stats from:" << endl;
  cout << "  " << STATS_IN + ".pairs" << endl;
  longlong nreads; 
  vec<String> lib_names;
  vec<int> lib_sep, lib_sd; 
  ReadPairsManagerLibInfo(STATS_IN + ".pairs", nreads, lib_names, lib_sep, lib_sd );

  cout << Date() << ": Original Stats:" << endl;
  pairs.printLibraryStats( cout );

  for (size_t lib_id = 0; lib_id < pairs.nLibraries(); ++lib_id) {
    int pos = Position(lib_names, pairs.getLibraryName(lib_id));
    if (pos != -1)
      pairs.changeLibrarySepSd(lib_id, lib_sep[pos], lib_sd[pos]);
    else
      FatalErr("Could not find replacement stats for library " 
	       + pairs.getLibraryName(lib_id));
  }

  cout << Date() << ": Replacement Stats:" << endl;
  pairs.printLibraryStats( cout );

   // Write corrected file
  cout << Date() << ": Writing new pairs file: " << endl 
       << "  " << PAIRS_OUT << ".pairs" << endl;
  pairs.Write(PAIRS_OUT + ".pairs");

}
