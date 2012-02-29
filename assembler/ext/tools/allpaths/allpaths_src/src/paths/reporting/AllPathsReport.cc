///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "FastIfstream.h"
#include "MainTools.h"
#include "Superb.h"
#include "efasta/EfastaTools.h"
#include "math/Functions.h"
#include "reporting/PerfStat.h"

/**
 * AllPathsReport
 *
 * Gather and/or generate statistics on an ALLPATHS assembly.
 *
 * MIN_CONTIG: min contig length (used in some statistics)
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_String_OrDefault( SUBDIR, "test" );
  CommandArgument_String_OrDefault( ASSEMBLY, "final" );
  CommandArgument_Int_OrDefault( MIN_CONTIG, 1000 );
  EndCommandArguments;
  
  // Dir and file names.
  String run_dir = PRE + "/" + DATA + "/" + RUN;
  String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;

  String accuracy_report_file = sub_dir + "/assembly_accuracy.report";
  String supers_file = sub_dir + "/" + ASSEMBLY + ".superb";
  
  // Super stats.
  vec<superb> supers;
  ReadSuperbs( supers_file, supers );

  int ntigs = 0;
  for (int ii=0; ii<supers.isize( ); ii++)
    ntigs += supers[ii].Ntigs( );

  vec<int> selected_clens;
  vec<int> ungapped_slens;
  vec<int> gapped_slens;
  selected_clens.reserve( ntigs );
  ungapped_slens.reserve( supers.size( ) );
  gapped_slens.reserve( supers.size( ) );

  for (int ii=0; ii<supers.isize( ); ii++) {
    int ungapped_super_len = 0;
    for (int cpos=0; cpos<supers[ii].Ntigs( ); cpos++) {
      if ( supers[ii].Len( cpos ) < MIN_CONTIG ) continue;
      selected_clens.push_back( supers[ii].Len( cpos ) );
      ungapped_super_len += supers[ii].Len( cpos );
    }
    if ( ungapped_super_len > 0 )
      ungapped_slens.push_back( ungapped_super_len );
    gapped_slens.push_back( supers[ii].TrueLength( ) );
  }
  sort( selected_clens.begin( ), selected_clens.end( ) );
  sort( ungapped_slens.begin( ), ungapped_slens.end( ) );
  sort( gapped_slens.begin( ), gapped_slens.end( ) );
  
  PerfStat::log( ) << PerfStat( "ap_report_min_contig",
				"contig minimum size for reporting",
				MIN_CONTIG );

  PerfStat::log( ) << PerfStat( "n_contigs",
				"number of contigs",
				selected_clens.size( ) );
  
  longlong est_genome_size = BigSum( gapped_slens );
  double c_ratio = double( selected_clens.size( ) )  / double( est_genome_size );
  PerfStat::log( ) << std::fixed << std::setprecision( 1 )
		   << PerfStat( "contigs_per_Mb", 
				"number of contigs per Mb",
				1000000.0 * c_ratio );

  PerfStat::log( ) << PerfStat( "n_scaffolds", 
				"number of scaffolds",
				supers.size( ) );

  PerfStat::log( ) << std::fixed << std::setprecision( 0 )
		   << PerfStat( "contig_length",
				"total contig length",
				BigSum( selected_clens ) );

  PerfStat::log( ) << std::fixed << std::setprecision( 0 )
		   << PerfStat( "scaff_length_gap",
				"total scaffold length, with gaps",
				BigSum( gapped_slens ) );
  
  double clens_n50 = N50( selected_clens );
  PerfStat::log( ) << std::fixed << std::setprecision( 1 )
		   << PerfStat( "N50_contig",
				"N50 contig size in kb",
				clens_n50 / 1000.0 );
  
  double ungapped_super_n50 = N50( ungapped_slens );
  PerfStat::log( ) << std::fixed << std::setprecision( 0 )
		   << PerfStat( "N50_scaffold",
				"N50 scaffold size in kb",
				ungapped_super_n50 / 1000.0 );

  double gapped_super_n50 = N50( gapped_slens );
  PerfStat::log( ) << std::fixed << std::setprecision( 0 )
		   << PerfStat( "N50_scaff_gap",
				"N50 scaffold size in kb, with gaps",
				gapped_super_n50 / 1000.0 );
  
  double d_ratio = double( supers.size( ) )  / double( est_genome_size );
  PerfStat::log( ) << std::fixed << std::setprecision( 2 )
		   << PerfStat( "scaff_per_Mb", 
				"number of scaffolds per Mb",
				1000000.0 * d_ratio );
  
  // Gap stats.
  int tot_pos_gap = 0;
  int neg_gap_5 = 0;
  vec<int> gap_size;
  vec<int> gap_dev;
  gap_size.reserve( ntigs );
  gap_dev.reserve( ntigs );
  for (int ii=0; ii<supers.isize( ); ii++) {
    for (int cpos=0; cpos<supers[ii].Ntigs( )-1; cpos++) {
      int gap = supers[ii].Gap( cpos );
      int dev = supers[ii].Dev( cpos );
      gap_size.push_back( gap  );
      gap_dev.push_back( dev );
      if ( gap + 5 * dev < 0 ) neg_gap_5 += gap;
      if ( gap > 0 ) tot_pos_gap += gap;
    }
  }
  sort( gap_size.begin( ), gap_size.end( ) );
  sort( gap_dev.begin( ), gap_dev.end( ) );
  
  bool have_gaps = gap_size.nonempty();
  int median_gap_size = (have_gaps ? Median(gap_size) : 0);
  int median_gap_dev  = (have_gaps ? Median(gap_dev)  : 0);

  PerfStat::log( ) << std::fixed << std::setprecision( 0 )
		   << PerfStat( "median_gap",
				"median size of gaps in scaffolds",
				median_gap_size );

  PerfStat::log( ) << std::fixed << std::setprecision( 0 )
		   << PerfStat( "median_gap_dev",
				"median dev of gaps in scaffolds",
				median_gap_dev );
  
  double tot_slen = BigSum( gapped_slens );
  d_ratio = (double)tot_pos_gap / tot_slen;
  PerfStat::log( ) << std::fixed << std::setprecision( 2 )
		   << PerfStat( "frac_captured_gaps",
				"% of bases in captured gaps",
				d_ratio * 100.0 );
  
  d_ratio = double( - neg_gap_5 )  / tot_slen;
  PerfStat::log( ) << std::fixed << std::setprecision( 2 )
		   << PerfStat( "frac_negative_gaps",
				"% of bases in negative gaps (after 5 devs)",
				d_ratio * 100.0 );
  
     // Count ambiguous bases.

     String efasta_file = sub_dir + "/" + ASSEMBLY + ".contigs.efasta";
     vec<efasta> tigse;
     LoadEfastaIntoStrings( efasta_file, tigse );
     longlong total_bases = 0, ambiguous_bases = 0, ambiguities = 0;
     for ( size_t i = 0; i < tigse.size( ); i++ )
     {    const efasta& A = tigse[i];
          if ( A.Length1( ) < (int) MIN_CONTIG ) continue;
          total_bases += A.Length1( );
          ambiguous_bases += A.AmbCount( );
          ambiguities += A.Ambiguities( );    }
     double amb_base_frac = double(ambiguous_bases) / double(total_bases);
     PerfStat::log( ) << std::fixed << std::setprecision(2)
          << PerfStat( "amb_base_frac", 
               "%% of ambiguous bases", 10000.0 * amb_base_frac );
     double ambiguity_frac = double(ambiguities) / double(total_bases);
     PerfStat::log( ) << std::fixed << std::setprecision(2)
          << PerfStat( "ambiguity_frac", 
               "ambiguities per 10,000 bases", 10000.0 * ambiguity_frac );    }
