/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// DumpHyper
//
// Write out a HyperKmerPath in a text format independent of our code base.
// The edge sequences are written to a .fasta file, and the graph structure
// to a .graphml file.

#include "Basevector.h"
#include "MainTools.h"
#include "lookup/LookAlign.h"
#include "math/HoInterval.h"
#include "paths/AlignHyperKmerPath.h"
#include "paths/HyperKmerPath.h"
#include "paths/KmerBaseBroker.h"

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_String_OrDefault(INSTANCE, "");
     CommandArgument_String(SUBDIR);
     CommandArgument_String_OrDefault(WRUN, "run");
     CommandArgument_String_OrDefault(HYPER, "hyper");
     CommandArgument_String_OrDefault(FASTA, "");
     CommandArgument_String_OrDefault(FASTB, "");
     CommandArgument_String_OrDefault(GRAPHML, "");
     
     EndCommandArguments;

     // Set up directories.

     SUBDIR += INSTANCE;
     String data_dir = PRE + "/" + DATA;
     String run_dir = PRE + "/" + DATA + "/" + RUN;
     String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;
     String wdata_dir = sub_dir;
     String wrun_dir = sub_dir + "/" + WRUN;

     // Load data.

     String hkpFname = sub_dir + "/" + HYPER;
     cout << Date( ) << ": Loading HyperKmerPath..." << endl;
     HyperKmerPath h( hkpFname );
     int K = h.K( );
     cout << Date( ) << ": Loading KmerBaseBroker..." << endl;
     KmerBaseBroker* kbb = new KmerBaseBroker( wrun_dir, K );

     if ( FASTA == "" ) FASTA = hkpFname + ".fasta";
     if ( FASTB == "" ) FASTB = hkpFname + ".fastb";
     if ( GRAPHML == "" ) GRAPHML = hkpFname + ".graphml";

     cout << Date( ) << ": Dumping edge sequence to " << FASTB << " ..." << endl;
     h.DumpFastb( FASTB, *kbb );
     cout << Date( ) << ": Dumping edge sequence to " << FASTA << " ..." << endl;
     h.DumpFasta( FASTA, *kbb );
     cout << Date( ) << ": Dumping graph structure to " << GRAPHML << " ..." << endl;
     h.DumpGraphML( GRAPHML );
     cout << Date( ) << ": Done." << endl;
}

