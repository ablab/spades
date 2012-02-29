///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "paths/CLinFastavec.h"
#include "paths/FlattenHyperFastavector.h"
#include "paths/SquashHyperFastavector.h"

void FlattenHyperFastavector( ostream &log,
			      HyperFastavector &hfv,
			      const String base_out,
			      const String hyper_inter_file,
			      const String dump_hfv_head,
			      const bool NEW_ALGORITHM,
			      const bool INITIAL_SCAFFOLD_PER_CONTIG,
			      const int MAX_CELL_SIZE,
			      const int MIN_EDGE_TO_SAVE,
			      const int NUM_THREADS )
{
  String scaffolds_file = base_out + ".superb";
  String contigs_file = base_out + ".contigs.fasta";
  String assembly_file = base_out + ".assembly.fasta";

  CLinFastavec flatter( hfv, &log );
  if ( NEW_ALGORITHM ) flatter.SetNewAlgorithm( );
  if ( hyper_inter_file != "" ) flatter.SetBaseInterSave( hyper_inter_file );
  
  if ( dump_hfv_head != "" )
  {    Ofstream( fout, dump_hfv_head + ".old.fasta" );
       for ( int v = 0; v < hfv.N( ); v++ )
       {    for ( size_t j = 0; j < hfv.From(v).size(); j++ )
            {    int e = hfv.EdgeObjectIndexByIndexFrom( v, j );
                 int w = hfv.From(v)[j];
                 hfv.EdgeObject(e).Print( fout, "edge_" + BaseAlpha(e) + "[vert_"
                      + ToString(v) + "-->vert_" + ToString(w)
                      + "]" );    }    }
       Ofstream( dout, dump_hfv_head + ".old.dot" );
       hfv.PrintSummaryDOT0w( dout, True, False, True, NULL, True );    }

  // Iteratively flatten the HyperFastavector.

  log << Date() << ": Flattening HyperFastavector" << endl;
  while ( 1 ) {
    while ( 1 ) {
      int npops = flatter.FlattenCells( MAX_CELL_SIZE );
      if ( npops < 1 ) break;
    }

    int npopsB = flatter.FlattenBubbles( );

    int npopsF = flatter.FlattenFrayedEnds( );

    if ( npopsB < 1 && npopsF < 1 ) {
      flatter.ReportBubbles( );
      break;
    }
  }

  if ( dump_hfv_head != "" )
  {    Ofstream( fout, dump_hfv_head + ".new.fasta" );
       for ( int v = 0; v < hfv.N( ); v++ )
       {    for ( size_t j = 0; j < hfv.From(v).size(); j++ )
            {    int e = hfv.EdgeObjectIndexByIndexFrom( v, j );
                 int w = hfv.From(v)[j];
                 hfv.EdgeObject(e).Print( fout, "edge_" + BaseAlpha(e) + "[vert_"
                      + ToString(v) + "-->vert_" + ToString(w)
                      + "]" );    }    }
       Ofstream( dout, dump_hfv_head + ".new.dot" );
       hfv.PrintSummaryDOT0w( dout, True, False, True, NULL, True );    }
  
  log << Date() << ": Converting HyperFastavector to scaffolds" << endl;
  vec<superb> scaffolds;
  hfv.ConvertToSuperbs( scaffolds, NUM_THREADS, MIN_EDGE_TO_SAVE );
  int largest_contig = 0;
  for (size_t ii=0; ii<scaffolds.size( ); ii++)
    for (int jj=0; jj<scaffolds[ii].Ntigs( ); jj++)
      largest_contig = Max( largest_contig, scaffolds[ii].Tig( jj ) );

  if ( base_out == "" ) {
    log << Date( ) <<": Not saving output" << endl;
    return;
  }
  log << Date() << ": Writing output files" << endl;

  int ntigs = hfv.EdgeObjectCount();
  vec<Bool> rctig( ntigs, False );

  if (INITIAL_SCAFFOLD_PER_CONTIG) {
    // Temporary hack: break scaffolds into singleton scaffolds (one
    //  contig per scaffold). This is done to fix a problem with
    //  ConvertToSuperbs, which caused some contigs to be joined
    //  incorrectly across a repeat. We still need to run
    //  ConvertToSuperb to remove the tangles from the assembly.
  
    int n_select = 0;
    for (size_t ii=0; ii<scaffolds.size( ); ii++)
      n_select += scaffolds[ii].Ntigs( );

    vec<int> select;
    select.reserve( n_select );
    for (size_t ii=0; ii<scaffolds.size( ); ii++)
      for (int jj=0; jj<scaffolds[ii].Ntigs( ); jj++)
	select.push_back( scaffolds[ii].Tig( jj ) );
    sort( select.begin( ), select.end( ) );
    
    vec<superb> singletons( select.size( ) );
    for (int ii=0; ii<singletons.isize( ); ii++) {
      int contig_id = select[ii];
      int contig_len = hfv.EdgeObject( contig_id ).size( );
      singletons[ii].PlaceFirstTig( contig_id, contig_len );
    }
    WriteSuperbs( scaffolds_file, singletons );
    hfv.WriteScaffoldsFasta( assembly_file, singletons, rctig );
  } else {
    WriteSuperbs( scaffolds_file, scaffolds );
    hfv.WriteScaffoldsFasta( assembly_file, scaffolds, rctig );
  }
  
  Ofstream( contig_fasta, contigs_file );
  for (int j = 0; j <= largest_contig; ++j)
    hfv.EdgeObject(j).Print( contig_fasta, "contig_" + ToString(j) );
  contig_fasta.close();
  
}

