//////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// PatcherCottage.  Employee of UnipathPatcher and PostPatcher.  
// Does jobs on demand, at home in his cottage.

// MakeDepend: dependency QueryLookupTable

#include "MainTools.h"
#include "ParseSet.h"
#include "feudal/BinaryStream.h"
#include "feudal/FeudalTools.h"
#include "paths/PatcherCottageCore.h"
#include "paths/UnipathFixerTools.h"
#include "system/WorklistMP.h"
// MakeDepend: cflags OMP_FLAGS

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_Int(K);
     CommandArgument_String(JOINDATA);
     CommandArgument_String(ROOT);
     CommandArgument_String(data_dir);
     CommandArgument_Int(MAX_READS);
     CommandArgument_Int(MAX_JOINS);
     CommandArgument_Int(MIN_OVERLAP_END);
     CommandArgument_String(LOG);
     EndCommandArguments;

     // Define logging options.

     vec<String> log_options;
     ParseStringSet( LOG, log_options );
     #define REQUESTED(x) Member( log_options, String(x) )
     Bool log_all = REQUESTED( "ALL" );
     Bool log_some = log_all || REQUESTED( "SOME" );
     Bool log_align_reads = log_all || REQUESTED( "ALIGN_READS" );
     Bool log_aligns = log_some || REQUESTED( "ALIGNS" );
     Bool log_aligns_all = log_some || REQUESTED( "ALIGNS_ALL" );
     Bool log_assembly = log_some || REQUESTED( "ASSEMBLY" );
     Bool log_assembly1 = log_some || REQUESTED( "ASSEMBLY1" );
     Bool log_assembly2 = log_all || REQUESTED( "ASSEMBLY2" );
     Bool log_assembly3 = log_all || REQUESTED( "ASSEMBLY3" );
     Bool log_correct = log_all || REQUESTED( "CORRECT" );
     double verbosity = 0;
     if (log_assembly) verbosity = 0.1;
     if (log_assembly1) verbosity = 1;
     if (log_assembly2) verbosity = 2;
     if (log_assembly3) verbosity = 3;
     PatcherCottage_LogParams log_params( verbosity, log_correct, log_align_reads, 
          log_aligns, log_aligns_all );

     BinaryReader jdRdr(JOINDATA.c_str(),false);
     WorklistMP<PCottageWhichJoin,PCottageResults>::Client client;

     while ( !client.isQuittingTime() )
     {
         PCottageWhichJoin iMsg = client.getWork();

         jdRdr.seek(iMsg.offset);
         PCottageJoinData joinData;
         jdRdr.read(&joinData);

         PCottageResults results;
         results.itemNumber = iMsg.itemNumber;

         if ( !SameSizes(joinData.reads,joinData.quals) )
         {
             std::cout << "Cottage ignoring item number " << iMsg.itemNumber
                       << " at offset " << iMsg.offset
                       << ": reads and quals have different shapes."
                       << std::endl;
             results.report = "FAILED";
         }
         else
         {
             PatcherCottageCore( joinData.L, joinData.R, joinData.sep, joinData.dev,
                                 joinData.reads, joinData.quals, joinData.pairs,
                                 results.report, log_params, data_dir,
                                 results.startStop, K, MAX_READS, MAX_JOINS,
                                 MIN_OVERLAP_END, ROOT, iMsg.itemNumber,
                                 iMsg.u1, iMsg.u2, False );

             joinData.reads.swap(results.reads);
         }
         client.reportResults(results);
     }
}
