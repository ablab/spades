///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

const char *DOC =
  "Displays library statistics for a set of paired reads as defined in the "
  "associated .pairs or .pairto file. It does not compute library stats.";

#include "MainTools.h"
#include "Basevector.h"
#include "PairsManager.h"
#include "feudal/VirtualMasterVec.h"


int main( int argc, char *argv[] )
{
  RunTime( );
  BeginCommandArguments;
  CommandDoc(DOC);
  CommandArgument_String_Doc(READS,
    "Fastb file with associated .pairs file.");
  CommandArgument_Bool_OrDefault_Doc( QUICK, False,
    "Quick peek at library stats without loading the entire file." );
  CommandArgument_Bool_OrDefault_Doc( READ_STATS, True,
    "Compute per library read statistics." );
  EndCommandArguments;

  // Strip .fastb from filenames

  READS = READS.SafeBefore(".fastb");
  READS = READS.SafeBefore(".pairs");
  if (READS.EndsWith("."))
    READS = READS.RevBefore(".");

  if (QUICK) {  // 'n' dirty
    cout << "Loading pairing information (basic library stats only) from:" << endl;
    if (IsRegularFile(READS + ".pairs")) {
      cout << "  " << READS + ".pairs" << endl;

      longlong nreads; 

      vec<PM_LibraryStats> stats;
      ReadPairsManagerLibInfo(READS + ".pairs", nreads, stats);

      size_t nlibs = stats.size();

      cout << endl;
      cout << "Libraries     : " << nlibs << endl;
      cout << "Total reads   : " << nreads << endl;
      cout << endl;

      writeLibraryStats(cout, stats);

    } else {
      FatalErr("Could not find .pairs file");
    }
    exit(0);
  }

  // Load pairing information

  cout << "Loading pairing information from:" << endl;
  PairsManager pairs;
  if (IsRegularFile(READS + ".pairs")) {
    cout << "  " << READS + ".pairs" << endl;
    pairs.Read(READS + ".pairs");
    if (IsRegularFile(READS + ".fastb")) {
      size_t nreads = MastervecFileObjectCount(READS + ".fastb");
      if (nreads != pairs.nReads()) {
	cout << "Warning: Inconsistency between .pairs and .fastb files..." << endl;
	cout << "         Fastb file contains " << nreads << " reads" << endl;
	cout << "         Pairs file reports  " << pairs.nReads() << " reads" << endl;
      }
    }
  } else if (IsRegularFile(READS + ".pairto")) {
    cout << "  " << READS + ".pairto" << endl;
    size_t nreads = MastervecFileObjectCount(READS + ".fastb");
    pairs.ReadFromPairtoFile(READS + ".pairto", nreads);
  } else {
    FatalErr("Could not find either a .pairs or .pairto file");
  }

  // Display library stats

  size_t npairs = pairs.nPairs();
  size_t nlibs = pairs.nLibraries();
  size_t nreads = pairs.nReads();

  cout << endl;
  cout << "Libraries      : " << nlibs << endl;
  cout << "Pairs          : " << npairs << endl;
  cout << "Total reads    : " << nreads << endl;
  cout << "Paired reads   : " << npairs * 2 << endl;
  cout << "Unpaired reads : " << nreads - (npairs * 2) << endl;
  cout << endl;

  if (READ_STATS) 
    pairs.printLibraryStats( cout , READS + ".fastb");
  else
    pairs.printLibraryStats( cout );
   
}
