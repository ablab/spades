/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2012) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

const char *DOC =
  "Takes an IN_HEAD.contigs.<efasta|fastb|fasta> file and it's associated IN_HEAD.superb "
  "file and generates the following new files:\n"
  "OUT_HEAD.contigs.{efasta,fasta,fastb}\n"
  "OUT_HEAD.assembly.{efasta,fasta}\n"
  "OUT_HEAD.superb";


#include "CoreTools.h"
#include "MainTools.h"
#include "system/System.h"

#include "Fastavector.h"
#include "efasta/EfastaTools.h"


int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandDoc(DOC);
  CommandArgument_String(IN_HEAD);
  CommandArgument_String(OUT_HEAD);
  CommandArgument_String_OrDefault_Doc(CONTIG_PREFIX, "contig_",
    "Add this string to the contig number to create the contig name");
  CommandArgument_Bool_OrDefault_Doc(REORDER, False,
    "If True, then reorder contigs to match scaffolds - removing any orphan contigs in the process.");
  EndCommandArguments;

  // Load data from efasta or fastb or fasta
  
  vec<efasta> econtigs;
  if ( IsRegularFile( IN_HEAD + ".contigs.efasta" ) ) {           // Use efasta

    cout << Date() << " : Loading efasta..." << endl;
    LoadEfastaIntoStrings( IN_HEAD + ".contigs.efasta", econtigs );

  } else if (IsRegularFile( IN_HEAD + ".contigs.fasta" ) ) {      // Use fasta

    cout << Date() << " : Loading fasta..." << endl;
    vec<fastavector> fastas;
    LoadFromFastaFile(  IN_HEAD + ".contigs.fasta", fastas );
    econtigs.resize( fastas.size() );
    for ( size_t i = 0; i < fastas.size(); i++ )
      econtigs[i] = efasta( fastas[i] );

  } else if (IsRegularFile( IN_HEAD + ".contigs.fastb" ) ) {      // Use fastb

    cout << Date() << " : Loading fastb..." << endl;
    vecbasevector bases( IN_HEAD + ".contigs.fastb");
    econtigs.resize( bases.size() );
    for ( size_t i = 0; i < bases.size( ); i++ )
      econtigs[i] = fastavector(bases[i]);
    
  } else 
    FatalErr("Could not find an efasta, fastb or fasta file: " + IN_HEAD + ".contigs.<efasta|fastb|fasta>" );


  cout << Date() << " : Loading superb..." << endl; 
  vec<superb> scaffolds;
  ReadSuperbs(IN_HEAD + ".superb", scaffolds);

  size_t n_scaffolds = scaffolds.size();
  size_t n_contigs = 0;
  for (size_t i = 0; i < scaffolds.size( ); i++)
    n_contigs += scaffolds[i].Ntigs( );
  cout << Date() << " : Assembly contains " << n_contigs << " contigs in " 
       << n_scaffolds << " scaffolds." << endl; 

  // Sanity Check
  if (n_contigs != econtigs.size())
    cout << Date() << " : WARNING: superb and efasta contig count inconsistency found" << endl;

  vec<uint32_t> used(econtigs.size(),0);
  vec<efasta> econtigs_reorder;
  if (REORDER) {
    // Renumber contigs according to scaffold position and reorder efasta accordingly

    cout << Date() << " : Re-ordering contigs to match scaffold" << endl;

    econtigs_reorder.resize(n_contigs);
    size_t new_index = 0;
    for (size_t scaffold_index = 0; scaffold_index < n_scaffolds; scaffold_index++) {
      for (size_t tig_index = 0; tig_index < static_cast<size_t>(scaffolds[scaffold_index].Ntigs( ) ); tig_index++) {
	size_t old_index = scaffolds[scaffold_index].Tig(tig_index);
	scaffolds[scaffold_index].SetTig(tig_index, new_index);
	econtigs_reorder[new_index] = econtigs[old_index];
	used[old_index]++;
	new_index++;
      }
    }
  }

  vec<efasta>& econtigs_out = (REORDER ? econtigs_reorder : econtigs);

  // Write output
  {
    cout << Date() << " : Writing new files to : " 
	 << OUT_HEAD << ".*" << endl;
    
    vec<FastaVec> flattened_fasta(econtigs_out.size());
    vecbasevector flattened_fastb(econtigs_out.size());
    
    {
      Ofstream(out_e, OUT_HEAD + ".contigs.efasta");
      Ofstream(out_a, OUT_HEAD + ".contigs.fasta");
      for (size_t it = 0; it < econtigs_out.size(); it++) {
	econtigs_out[it].FlattenTo(flattened_fasta[it]);
	econtigs_out[it].FlattenTo(flattened_fastb[it]);
	econtigs_out[it].Print(out_e, CONTIG_PREFIX + ToString(it));
	flattened_fasta[it].Print(out_a, CONTIG_PREFIX + ToString(it));
      }
      flattened_fastb.WriteAll( OUT_HEAD + ".contigs.fastb" );  
    }

    WriteScaffoldedEFasta( OUT_HEAD + ".assembly.efasta", econtigs_out, scaffolds );
    WriteScaffoldedFasta( OUT_HEAD + ".assembly.fasta", flattened_fasta, scaffolds );
    WriteSuperbs( OUT_HEAD + ".superb", scaffolds );
    WriteSummary( OUT_HEAD + ".summary", scaffolds );
  }
 
  cout << Date( ) << " : Done." << endl;
}

