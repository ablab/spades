///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "lookup/QueryLookupTableCore.h"

int main( int argc, char *argv[] )
{    
     // Define interrupt handling.

     Bool TRACEBACK_ON_INTERRUPT = False;
     for ( int i = 1; i < argc; i++ )
     {    if ( String(argv[i]) == String("TRACEBACK_ON_INTERRUPT=True") )
               TRACEBACK_ON_INTERRUPT = True;    }
     RunTime(1,
             (TRACEBACK_ON_INTERRUPT
              ? &arachne_signal_handler_standard 
              : &arachne_signal_handler_no_ctrlc_traceback));

     // Run QueryLookupTable.

     QueryLookupTableCore( argc, argv );    }
