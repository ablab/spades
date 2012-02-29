///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/**
   Contains additional code for RunAllPathsLG.
*/

#include "paths/RunAllPathsCommon.h"

#include "MainTools.h"
#include "Basevector.h"
#include "math/Functions.h"
#include "system/SysConf.h"


void UpdateRefGenomeSize(String ref_dir) {
  if (!IsRegularFile(ref_dir + "/genome.size") || 
      IsOlder(ref_dir + "/genome.size", ref_dir + "/genome.fastb") ) {
    longlong genome_size = 0;
    vecbasevector genome(ref_dir + "/genome.fastb");
    for ( size_t i = 0; i < genome.size( ); i++ )
      genome_size += genome[i].size( ); 
    Ofstream( genomeSizeFile, ref_dir + "/genome.size" );
    genomeSizeFile << genome_size;
  }
}


// Set module thread count to global value or override with module value

void SetThreadValue(int& module_value , const int global_value) {
  module_value = (-1 == module_value ? global_value : module_value);
}


// Sets threading values for each module

int ParseThreadArgs(const String& THREADS,
		    int& LR_THREADS, 
		    int& CP_THREADS, 
		    int& FF_THREADS,
		    int& ECJ_THREADS,
		    int& KP_THREADS,
		    int& CUG_THREADS, 
		    int& RDR_THREADS, 
		    int& FE_THREADS,
		    int &MN_THREADS,
		    int &PR_THREADS,
		    int &PC_THREADS) 
{
   // Determine global thread count and set special cases

  int global_thread_count = 1;
  if (THREADS == "max") {      // use all cpus/cores
    global_thread_count = processorsOnline();

  } else if (THREADS.IsInt())       // used fixed value
    global_thread_count = THREADS.Int();

  else {
    cout << "Input Error: invalid value for THREADS: " << THREADS << endl;
    exit(1);
  }

  // Set all other thread counts to global value

  SetThreadValue(LR_THREADS, global_thread_count);
  SetThreadValue(CP_THREADS, global_thread_count);
  SetThreadValue(FF_THREADS, global_thread_count);
  SetThreadValue(ECJ_THREADS, global_thread_count);
  SetThreadValue(KP_THREADS, global_thread_count);
  SetThreadValue(CUG_THREADS, global_thread_count);
  SetThreadValue(RDR_THREADS, global_thread_count);
  SetThreadValue(FE_THREADS, global_thread_count);
  SetThreadValue(MN_THREADS, global_thread_count);
  SetThreadValue(PR_THREADS, global_thread_count);
  SetThreadValue(PC_THREADS, global_thread_count);
  
  return global_thread_count;
}


void DisplaySystemStatus() {

  const int giga = 1024 * 1024 * 1024;
  cout << endl;
  cout << "Hostname                  : " << StringOfOutput( "hostname" ) << endl; 
  cout << "Available processors      : " << processorsOnline() << endl;
  cout << "Total physical memory     : " << physicalMemory() / giga << " GB" << endl;
  cout << endl;

}
