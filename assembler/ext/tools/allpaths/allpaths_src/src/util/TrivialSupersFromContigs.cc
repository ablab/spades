/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "FetchReads.h"
#include "Loader.h"
#include "MainTools.h"
#include "Superb.h"

/**
 * GenerateTrivialSupers
 *
 * Generate a trivial supercontigs structure (each super a singleton) in
 * a given SUBDIR. Supers are guaranteed to inherit the id from their only
 * contig. It only needs the mergedcontigs.fastb file in SUBDIR. Notice
 * that it will silently overwrite an existing super structure.
 *
 * MC: it will load MC.fastb and save MC.superb* (in PRE/DATA/SUBDIR)
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( HEAD );
  EndCommandArguments;

  // Dir and file names.
  String contigs_file = HEAD + ".contigs.fasta";
  String superb_file = HEAD + ".superb";
  
  // Load contigs.
  cout << Date( ) << ": loading contigs" << endl;
  vecbvec contigs;
  FetchReads( contigs, 0, contigs_file );
  
  // Fill supers.
  vec<superb> supers( contigs.size( ) );
  for (int ii=0; ii<(int)supers.size( ); ii++) {
    supers[ii].PlaceFirstTig(ii, contigs[ii].size());
  }

  // Save.
  WriteSuperbs( superb_file, supers );

}
