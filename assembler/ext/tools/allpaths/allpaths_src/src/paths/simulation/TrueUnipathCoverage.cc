/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// TrueUnipathCoverage.  
// Count the number of times each unipath appears in the genome,
// then make a MasterVec< SerfVec< pair<copy_num_t,prob_t> > > in the same
// format as UnipathCoverage (the pair <n,p> means that the probability that the
// copy # is n is p), except that here we know the true copy number so the
// probability for the true copy # is 1 while the probability for any other
// copy # is 0.
// Save as reads.unipaths.true_count.k*, to which predicted_count
// may be symlinked.  (Or maybe we should add in a gaussian blur?)

#include "MainTools.h"
#include "CommonSemanticTypes.h"
#include "feudal/MasterVec.h"
#include "feudal/SerfVec.h"
#include "feudal/OuterVecDefs.h"
#include "feudal/SmallVecDefs.h"
#include "paths/KmerPath.h"

typedef MasterVec< SerfVec< pair<copy_num_t,prob_t> > > VecPairVec;

int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String(PRE);
  CommandArgument_String(DATA);
  CommandArgument_String(RUN); 
  CommandArgument_Int(K);
  CommandArgument_String_OrDefault(READS, "genome");
  EndCommandArguments;
  
  // Set up directories.

  String datadir = PRE + "/" + DATA;
  String run_dir = PRE + "/" + DATA + "/" + RUN;

  // Read in pathsdb and unipaths.

  BREAD2( run_dir + "/"+READS+".pathsdb_big.k" + ToString(K),
          vec<big_tagged_rpint>, pathsdb );
  vecKmerPath unipaths( run_dir + "/"+READS+".unipaths.k" + ToString(K) );
  int nuni = unipaths.size();

  // The thing we're going to generate:
  VecPairVec pdf(nuni);

  vec<unipath_interval_id_t> hits;
  for( int i=0; i<nuni; i++ ) {
    if( unipaths[i].IsEmpty() )
      hits.clear();
    else
      Contains( pathsdb, unipaths[i].Start(0), hits );

    pdf[i].push_back( pair<copy_num_t,prob_t>( hits.size(), 1.0 ) );
  }

  pdf.WriteAll( run_dir + "/"+READS+".unipaths.true_count.k" 
		+ ToString(K) );    
}
