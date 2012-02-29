///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "String.h"
#include "SupersHandler.h"
#include "agp/AgpFile.h"
#include "system/ParsedArgs.h"
#include "system/RunTime.h"
#include "system/System.h"

/**
 * SuperbsToAgp
 *
 * It generates the AGP files from a given assembly. This is similar
 * to agp/SupersToAgp, but it only requires a superb structure. Output
 * is saved either as OUTPUT or as SUPERBS.agp
 *
 * OUTPUT: full path name for output file (default is <SUPERBS>.agp)
 * ACCESSIONS: if given, then use these as contig names
 * GAP_FLOOR: minimum gap size
 */
int main( int argc, char *argv[] )
{ 
  RunTime();
  
  BeginCommandArguments;
  CommandArgument_String( SUPERBS );
  CommandArgument_String_OrDefault( OUTPUT, "" );
  CommandArgument_String_OrDefault( ACCESSIONS, "" );
  CommandArgument_UnsignedInt_OrDefault( GAP_FLOOR, 100 );
  EndCommandArguments;

  // Output file.
  String out_file = OUTPUT == "" ? SUPERBS + ".agp" : OUTPUT;
  ofstream out( out_file.c_str( ) );
  
  // Load.
  cout << Date( ) << ": loading scaffolds" << endl;
  shandler supers( -1, SUPERBS );
  
  vec<String> acc_names;
  if ( ACCESSIONS != "" ) {
    cout << Date( ) << ": loading accession names" << endl;
    READX( ACCESSIONS, acc_names );
  }
  
  // Loop over the supers.
  cout << Date( ) << ": loop over " << supers.Size( ) << " scaffolds" << endl;
  for (int ii=0; ii<supers.Size( ); ii++) {
    agp_chromosome agp( "_" + ToString( ii ) ) ;
    
    // Loop over the contigs.
    for (int jj=0; jj<supers[ii].Ntigs( ); jj++) {

      // Add contig.
      String cg_name;
      if ( acc_names.size( ) > 0 ) cg_name = acc_names[ supers[ii].Tig( jj ) ];
      else cg_name = "contig_" + ToString( supers[ii].Tig( jj ) );

      int cg_len = supers[ii].Len( jj );
      int cg_start = 0;
      int cg_stop = cg_len - 1;
      bool cg_RC = false;
      agp_contig new_contig( cg_name, cg_len, cg_start, cg_stop, cg_RC );
      new_contig.SetType( agp_contig::wgs_contig );
      agp.AddContig ( new_contig );

      // Add gap.
      if ( jj < supers[ii].Ntigs( ) - 1 ) {
	String gap_type = "fragment";
	bool is_bridged = true;
	int gap_len = supers[ii].Gap( jj );
	if ( gap_len < (int)GAP_FLOOR ) gap_len = GAP_FLOOR;
	agp_gap new_gap( gap_type, gap_len, is_bridged );
	agp.AddGap( new_gap );
      }
    }
    
    // Print info for super.
    String base_name = "scaffold";
    agp.Print( out, &base_name );
  }
  
  // Done.
  cout << Date( ) << ": done" << endl;
  out.close( );
  
}
