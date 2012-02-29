/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// Program: PrepareGenome
//
// Make the unipathed genome for analysis purposes.
//
// Input file: genome.fastb
//
// Output files: genome.paths{,_rc,db_big}.k*
//               genome.unipaths{,db_big}.k*
//               genome.unipaths.placements.k*
//               genome.unipaths.true_count.k*
//
// Part of <reference genome analysis>.

#include "MainTools.h"
#include "CommonSemanticTypes.h"
#include "system/MiscUtil.h"
#include "paths/KmerPath.h"
#include "paths/KmerBaseBroker.h"

int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String(PRE);
  CommandArgument_String(DATA);
  CommandArgument_Int(K);
  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
    "Number of threads to use (use all available processors if set to 0)");
  CommandArgument_Bool_OrDefault_Doc(CANONICALIZE, False,
				     "canonicalize the genome paths" );
  EndCommandArguments;
  
  // Thread control
   
  NUM_THREADS = configNumThreads(NUM_THREADS);

  String datadir = PRE + "/" + DATA;
  String kpd = " K=" + ToString(K) + " PRE=" + PRE + " DATA=" + DATA + " ";
  String kpdr = kpd + "RUN=. ";
  String dotK = ".k" + ToString(K);

  // Run needed executables
  // Yes, this could be done by a makefile, thank you for asking.

// MakeDepend: dependency GenomeToPaths
  Make( FilesIn( datadir, "genome.paths" + dotK ),
	FilesIn( datadir, "genome.fastb" ),
	"GenomeToPaths HEAD="+datadir+"/genome K=" + ToString(K)
	+ " CANONICALIZE=" + ( CANONICALIZE ? "True" : "False" ) 
	+ " NUM_THREADS= " + ToString(NUM_THREADS) );

// MakeDepend: dependency MakeRcDb
  Make( FilesIn( datadir, "genome.paths_rc" + dotK, "genome.pathsdb_big" + dotK ),
	FilesIn( datadir, "genome.paths" + dotK ),
	"MakeRcDb" + kpdr + "READS=genome FORCE_BIG=True" );

// MakeDepend: dependency UnipatherBig
  Make( FilesIn( datadir, "genome.unipaths" + dotK, "genome.unipathsdb_big" + dotK ),
	FilesIn( datadir, "genome.paths" + dotK, "genome.paths_rc" + dotK,
		 "genome.pathsdb_big" + dotK ),
	"UnipatherBig" + kpdr + "READS=genome" );

// MakeDepend: dependency TrueUnipathCoverage
  Make( FilesIn( datadir, "genome.unipaths.true_count" + dotK ),
	FilesIn( datadir, "genome.unipaths" + dotK ),
	"TrueUnipathCoverage" + kpdr + "READS=genome" );

// MakeDepend: dependency PathsToLocs
  Make( FilesIn( datadir, "genome.unipaths.placements" + dotK ),
	FilesIn( datadir, "genome.unipaths" + dotK ),
	"PathsToLocs" + kpdr 
	+ "PATHS=" + datadir + "/genome.unipaths" + dotK
	+ " PATHS_ABSOLUTE=True "
	+ " SAVE_PLACEMENTS_TO=" + datadir 
	+ "/genome.unipaths.placements" + dotK
	+ " SHOW_PLACEMENTS=True PRINT_ADJACENCY_GRAPH=False" );

  // Create genome unibases
  if ( NeedToMake( FilesIn( datadir, "genome.unibases" + dotK ),
		   FilesIn( datadir, "genome.unipaths" + dotK,
			    "genome.unipathsdb_big" + dotK ) ) ) {
    cout << Date( ) << ": Making genome unibases" << endl;
    BREAD2( datadir + "/genome.pathsdb_big" + dotK,
	    vec<big_tagged_rpint>, genome_pathsdb );
    vecKmerPath genome_unipaths( datadir + "/genome.unipaths" + dotK );
    KmerBaseBrokerBig kbb( datadir, K, "genome" );
    
    // Reserve memory
    vecbasevector genome_unibases;
    size_t n_paths = genome_unipaths.size( );
    longlong n_bases = 0;
    for ( size_t i = 0; i < n_paths; i++ )
      n_bases += genome_unipaths[i].KmerCount( ) + K - 1;
    genome_unibases.Reserve( n_bases / 16 + n_paths, n_paths );
    
    for ( size_t i = 0; i < n_paths; i++ )
      genome_unibases.push_back( kbb.Seq( genome_unipaths[ i ] ) );
    genome_unibases.WriteAll( datadir + "/genome.unibases" + dotK );
    cout << Date( ) << ": Made genome unibases" << endl;
  }

  CpIfNewer( datadir + "/genome.unibases" + dotK, datadir + "/genome.unibases" + dotK + ".fastb" );
  
// MakeDepend: dependency MakeLookupTable
  Make( FilesIn( datadir, "genome.unibases" + dotK + ".lookup" ),
	FilesIn( datadir, "genome.unibases" + dotK + ".fastb" ),
	"MakeLookupTable SOURCE=" + datadir + "/genome.unibases" + dotK + ".fastb" + 
	" OUT_HEAD=" + datadir + "/genome.unibases" + dotK + " LOOKUP_ONLY=True" );
  
}


