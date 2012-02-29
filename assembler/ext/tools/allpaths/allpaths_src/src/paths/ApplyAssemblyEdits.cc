///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// ApplyAssemblyEdits.
// Read multiple suggested edits of an assembly and perform the edits.
// Output the assembly to efasta format.

#include "MainTools.h"
#include "Basevector.h"
#include "Fastavector.h"
#include "efasta/EfastaTools.h"
#include "paths/AssemblyEdit.h"

int main(int argc, char *argv[])
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String(PRE);
  CommandArgument_String(DATA);
  CommandArgument_String(RUN);
  CommandArgument_String_OrDefault(SUBDIR, "test");
  CommandArgument_String_OrDefault(SCAFFOLDS_IN, "linear_scaffolds0.patched");
  CommandArgument_String_OrDefault(SCAFFOLDS_OUT, "linear_scaffolds0.patched.fixed");
  CommandArgument_StringSet_OrDefault(EDITS_IN,
      "linear_scaffolds0.patched.fixed.edits, linear_scaffolds0.patched.fixed2.edits");
  CommandArgument_Bool_OrDefault(VERBOSE, False);
  EndCommandArguments;

  // Begin.
  double clock = WallClockTime();

  // Define directories, etc.
  String data_dir = PRE + "/" + DATA, run_dir = PRE + "/" + DATA + "/" + RUN;
  String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;
  String out_head = sub_dir + "/" + SCAFFOLDS_OUT;
  String scaffolds_tig_file = sub_dir + "/" + SCAFFOLDS_IN + ".contigs.vecfasta";

  // Read all contigs.
  const size_t n_tigs = MastervecFileObjectCount(scaffolds_tig_file);
  FastaVecVec contigs;
  contigs.ReadAll(scaffolds_tig_file);
  ForceAssertEq( n_tigs, contigs.size() );

  // Read all edits, those that overlap with previous edits will not be performed

  vec< int > left_trim( n_tigs, 0 ), right_trim( n_tigs, 0 ); // record all trims
  vec< vec<assembly_edit> > all_edits( n_tigs );
  vec< String > input_edits ( EDITS_IN.size() );
  for ( size_t iFile = 0; iFile < input_edits.size(); iFile ++ )
  {
    String filename =  sub_dir + "/" + EDITS_IN[iFile];
    ForceAssert( IsRegularFile( filename ) );
    // read the edits
    vec<assembly_edit> editsi;
    BinaryReader::readFile( filename.c_str( ), &editsi );
    cout << " read edits from " << EDITS_IN[iFile] << endl;
    for ( int j = 0; j < editsi.isize( ); j++ )
    { const assembly_edit& e = editsi[j];
      ForceAssert( e.Internal( ) );
      int tigId = e.Tig( );
      // check conflicts
      bool overLap = false;
      for ( size_t i = 0; i < all_edits[tigId].size(); i++ )
      {    if ( Overlap( all_edits[tigId][i], e ) )
	   {
	     overLap = true;
	     break;
	   }
      }
      if ( overLap ) {
	if (VERBOSE) 
	  cout << " Overlapping edits skipped: " << e << endl;
      } else {
	if ( VERBOSE) 
	  cout << " read edits: " << e << endl;
	if ( e.PureDeletion( ) && e.Start( ) == 0 ) { // left trimming
	  left_trim[tigId] = e.Stop( );
	}
	if ( e.PureDeletion( ) && (unsigned) e.Stop( ) == contigs[tigId].size() ) { // right trimming
	  right_trim[tigId] = e.Stop( ) - e.Start( );
	}
	all_edits[tigId].push( e );
      }
    }
  }

  // Make the edits for each contigs
  vec<EFastaVec> econtigs_new(n_tigs);
  for ( size_t tig = 0; tig < contigs.size(); tig++ )
  {
    vec<assembly_edit> edits = all_edits[tig];
    Sort( edits );
    FastaVec& contig = contigs[tig];
    size_t contig_size = contig.size(), ie = 0;
    for ( size_t i = 0; i < contig_size; i++ ) 
    { 
      if ( ie < edits.size() && i == (unsigned)edits[ie].Start( ) )
      {    
        efasta e( edits[ie].Reps( ) );
        econtigs_new[tig].append(e);
	i += edits[ie].Stop( ) - edits[ie].Start( ) - 1;
	ie++;    
      }
      else  
	if (  contig[i] == 'A' || contig[i] == 'C' || contig[i] == 'G' || contig[i] == 'T' )
	  econtigs_new[tig].push_back( contig[i] );    
	else
	  econtigs_new[tig].append( ExpandAmbCode(contig[i]) );  
    }
  }

  // Write output

  cout << "\n" << Date( ) << ": Writing fixed efasta files to : " 
    << out_head << ".*" << endl;
  vec<FastaVec> flattened_fasta(econtigs_new.size());
  vec<FastaVec> flattened_fasta_max(econtigs_new.size());
  vecbasevector flattened_fastb(econtigs_new.size());
  {    
    Ofstream(out_e, out_head + ".contigs.efasta");
    Ofstream(out_a, out_head + ".contigs.fasta");
    Ofstream(out_m, out_head + ".contigs.max.fasta");
    for (size_t it = 0; it < n_tigs; it++) 
    {    
      econtigs_new[it].FlattenTo(flattened_fasta[it]);
      econtigs_new[it].FlattenNMaxTo(flattened_fasta_max[it]);
      econtigs_new[it].FlattenTo(flattened_fastb[it]);
      econtigs_new[it].Print(out_e, ToString(it));
      flattened_fasta[it].Print(out_a, ToString(it));   
      flattened_fasta_max[it].Print(out_m, ToString(it));  
    }
    flattened_fastb.WriteAll( out_head + ".contigs.fastb" );     
  }
  vecfastavector flattened_fasta2(econtigs_new.size());
  for (size_t it = 0; it < n_tigs; it++)
    flattened_fasta2.push_back_reserve( flattened_fasta[it] );
  flattened_fasta2.WriteAll( out_head + ".contigs.vecfasta" );
  vec<superb> scaffolds;
  ReadSuperbs(sub_dir + "/" + SCAFFOLDS_IN + ".superb", scaffolds);
  for (int ii=0; ii<scaffolds.isize( ); ii++) 
  { 
    superb& S = scaffolds[ii];
    for ( int jj=0; jj < S.Ntigs( ); jj++ ) 
    {
      int tig_id = S.Tig( jj );
      int tig_len = flattened_fastb[tig_id].size( );
      S.SetLen( jj, tig_len );    
      if ( jj > 0 ) S.SetGap( jj-1, S.Gap(jj-1) + left_trim[tig_id] );
      if ( jj < S.Ngaps( ) - 1 ) 
	S.SetGap( jj, S.Gap(jj) + right_trim[tig_id] );    
    }    
  }
  WriteScaffoldedEFasta(out_head + ".assembly.efasta", econtigs_new, scaffolds);
  WriteScaffoldedFasta(out_head + ".assembly.fasta", flattened_fasta, scaffolds);
  WriteScaffoldedFasta(out_head + ".assembly.max.fasta", flattened_fasta_max, scaffolds);
  WriteSuperbs( out_head + ".superb", scaffolds );
  WriteSummary( out_head + ".summary", scaffolds );
  // Done.
  cout << "time used = " << TimeSince(clock) << endl;    
}
