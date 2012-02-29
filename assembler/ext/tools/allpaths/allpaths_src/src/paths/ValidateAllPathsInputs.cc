///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

const char *DOC =
  "Validates input data for ALLPATHS-LG. Checks that supplied libraries fit "
  "our requirements.";
  
#include "MainTools.h"
#include "String.h"
#include "Vec.h"
#include "Basevector.h"
#include "feudal/BinaryStream.h"
#include "feudal/VirtualMasterVec.h"
#include "PairsManager.h"


enum LibType {frag,jump,long_jump};
  
void ErrorMessage(const String message) {
  cout << "[ERROR] " << message << endl; 
}

bool IsFile(const String filename) {
  bool found = IsRegularFile(filename);
  if (!found)
    ErrorMessage("Cannot find file: " + filename);
  return found;
}

bool ErrorIfTrue(bool test, const String message) {
  if (test)
    ErrorMessage(message);
  return test;
}

bool ErrorIfFalse(bool test, const String message) {
  return ErrorIfTrue(!test, message);
}


// Performs various tests on a set of libraries and returns true if no problems were found

bool ValidateLibraries(const LibType lib_type, const String file_head, const int K, ostream &log ) {
  bool error = false;

  error |= !IsFile(file_head + ".qualb");
  error |= !IsFile(file_head + ".pairs"); 
  error |= !IsFile(file_head + ".fastb");

  if (error) // Missing files, unable to continue
    return false;
		   
  PairsManager pairs;
  pairs.Read(file_head + ".pairs");
  
  size_t npairs = pairs.nPairs();
  size_t nlibs = pairs.nLibraries();
  size_t nreads = pairs.nReads();
  
  if (ErrorIfTrue( nlibs == 0, "No libraries founds.") )  // Empty files, unable to continue
    return false;

  vec<PM_LibraryStats> stats = pairs.getLibraryStats(file_head + ".fastb");

  // Display Library Stats

  cout << "Libraries     : " << nlibs << endl;
  cout << "Pairs         : " << npairs << endl;
  cout << "Total reads   : " << nreads << endl;
  cout << "Paired reads  : " << npairs * 2 << endl;
  cout << "Unpaired reads: " << nreads - (npairs * 2) << endl;
  cout << endl;

  writeLibraryStats( cout, stats );

  // Save output as a csv table.

  String str_type;
  if ( lib_type == frag )      str_type = "frag";
  else if ( lib_type == jump ) str_type = "jump";
  else                         str_type = "long_jump";

  log << "lib_type,lib_name,tot_bases,min_read_len,max_read_len,mean_read_len\n";
  for (size_t lib_id=0; lib_id<stats.size( ); lib_id++) {
    const PM_LibraryStats &lib = stats[lib_id];
    log << str_type << ","
	<< lib.name << ","
	<< lib.n_bases << ","
	<< lib.min_len << ","
	<< lib.max_len << ","
	<< lib.mean_len << "\n";
  }
  log << endl;

  // Basic consistency tests
  
  error |= ErrorIfTrue(MastervecFileObjectCount(file_head + ".fastb") != nreads,
				   "Inconsistency found between pairs manager and fastb");

  error |= ErrorIfTrue(MastervecFileObjectCount(file_head + ".qualb") != nreads,
		   "Inconsistency found between pairs manager and qualb");
  
  // General per library tests
  for ( size_t i = 0; i < nlibs; i++ ) {
    String libname =  stats[i].name;
    if (lib_type == frag) {
      error |= ErrorIfTrue(2*stats[i].max_len + stats[i].sep + 2 * stats[i].sd <= static_cast<unsigned>(K),
			   "Library " + libname + " contains no reads that overlap to create a super read larger K (" + ToString(K) + ").") ;
    }
  }
  
  // Fragment library tests
  if (lib_type == frag) {
    bool found_good_sep = false;
    for ( size_t i = 0; i < nlibs; i++ ) 
      found_good_sep |= (stats[i].sep - 2 * stats[i].sd <= 0 );
    error |= ErrorIfFalse(found_good_sep, "Could not find a fragment library whose reads overlapped.");
  }

  // Jumping library tests
  if (lib_type == jump) {
    bool found_good_sep = false;
    for ( size_t i = 0; i < nlibs; i++ ) 
      found_good_sep |= (stats[i].sep > 1000 && stats[i].sep < 10000 );
    error |= ErrorIfFalse(found_good_sep, "Could not find a jumping library with a separation > 1000 and < 10000 bases.");
  }

  // Long jumping library tests
  if (lib_type == long_jump) {
    bool found_good_sep = false;
    for ( size_t i = 0; i < nlibs; i++ ) 
      found_good_sep |= (stats[i].sep > 5000 );
    error |= ErrorIfFalse(found_good_sep, "Could not find a long jumping library with a separation > 5000.");
  }

  return !error;
}


int main( int argc, char *argv[] )
{
  RunTime( );
  BeginCommandArguments;
  CommandDoc(DOC);
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_Int(K);
  CommandArgument_Bool_OrDefault_Doc( WARN_ONLY, False,
    "Don't stop executation of pipeline if a problem is found.")
  CommandArgument_Bool_OrDefault_Doc( LONG_JUMPS, False,
    "Validate long jumping read libraries.")
  CommandArgument_String_OrDefault_Doc( REPORT, "ValidateAllPathsInputs_core_stats.out",
    "Save to file core stats (used by reporting code).")
  EndCommandArguments;

  String data_dir = PRE + "/" + DATA;

  String frag_reads_head = data_dir + "/frag_reads_orig";
  String jump_reads_head = data_dir + "/jump_reads_orig";
  String long_jump_reads_head = data_dir + "/long_jump_reads_orig";
  String core_stats_log_file = data_dir + "/" + REPORT;

  ofstream log( core_stats_log_file.c_str( ) );
  PrintCommandPretty( log );
  
  bool error = false;

  // Validate fragment reads

  cout << "Validating Fragment Libraries" << endl
       << "=============================" << endl << endl;

  error |= !ValidateLibraries(frag, frag_reads_head, K, log); 


  cout << endl
       << "Validating Jumping Libraries" << endl
       << "============================" << endl << endl;

  error |= !ValidateLibraries(jump, jump_reads_head, K, log); 


  if (LONG_JUMPS) {
    cout << endl
	 << "Validating Long Jumping Libraries" << endl
	 << "=================================" << endl << endl;
    
    error |= !ValidateLibraries(long_jump, long_jump_reads_head, K, log); 
  }

  cout << endl;
  log.close( );
  
  if (!error)
    cout << "Validation succeeded." << endl;
  if (error && WARN_ONLY) 
    cout << "Validation failed, but proceeding anyway." << endl;
  else if (error) {
    cout << "Validation failed. Aborting RunAllPathsLG." << endl
	 << "Please fix problems noted above and try again." << endl
	 << "For library requirements see the ALLPATHS-LG manual." << endl;
    exit(1);
  }
}
