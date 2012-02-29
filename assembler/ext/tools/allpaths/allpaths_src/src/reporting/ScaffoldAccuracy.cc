///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// ScaffoldAccuracy.  Given an assembly consisting of scaffolds, and a reference
// sequence for the genome, assess the accuracy of the scaffolds by determining the
// probability that positions at given distance D (default 100,000) apart are
// correctly connected.
//
// Also generate stats about assembly contiguity.
//
// READ THIS:
// - The assembly can be in fasta format, in which case gaps are
//   identified by lower case Ns, e.g. nnnnnnnnnnnnnnnnnnnn.
// - Or the assembly can be in efasta format.  This is determined by seeing if
//   the filename ends in .efasta.

#include "Alignment.h"
#include "Basevector.h"
#include "Bitvector.h"
#include "FetchReads.h"
#include "FetchReadsAmb.h"
#include "MainTools.h"
#include "Superb.h"
#include "efasta/EfastaTools.h"
#include "lookup/PerfectLookup.h"
#include "math/Functions.h"
#include "random/Random.h"
#include "reporting/PerfStat.h"
#include "reporting/ScaffoldAccuracyCore.h"

int main(int argc, char **argv)
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_Doc(ASSEMBLY, 
          "fasta file for assembly, with ns for gaps");
     CommandArgument_String_Doc(REFHEAD, 
          "files REFHEAD.fasta and REFHEAD.lookup must exist" );
     CommandArgument_Int_OrDefault(D, 100000);
     CommandArgument_Int_OrDefault(MAX_ERR, 25000);
     CommandArgument_Int_OrDefault(SAMPLE, 100000);
     CommandArgument_Int_OrDefault_Doc(MIN_CONTIG, 1000,
          "minimum contig size to be used in contiguity estimates");
     CommandArgument_Bool_OrDefault(VERBOSE, False);
     CommandArgument_Bool_OrDefault(PERF_STATS, False);
     CommandArgument_Bool_OrDefault(ALLOW_N, False);
     EndCommandArguments;

     // Check for existence of files.
     
     ForceAssert( IsRegularFile(ASSEMBLY) );
     ForceAssert( IsRegularFile( REFHEAD + ".fasta" ) );
     ForceAssert( IsRegularFile( REFHEAD + ".lookup" ) );
     String supers_file = ASSEMBLY;
     if ( ASSEMBLY.Contains( ".assembly.efasta" ) )
       supers_file = supers_file.Before( ".assembly.efasta" );
     else if ( ASSEMBLY.Contains( ".assembly.fasta" ) )
       supers_file = supers_file.Before( ".assembly.fasta" );
     supers_file += ".superb";
     
     // Load the reference and assembly.

     vecbasevector ref, assembly;
     vecbitvector ref_amb, assembly_amb;
     FetchReads( ref, 0, REFHEAD + ".fasta" );
     FetchReadsAmb( ref_amb, REFHEAD + ".fasta" );
     if ( !ASSEMBLY.Contains( ".efasta", -1 ) )
     {    FetchReads( assembly, 0, ASSEMBLY );
          if ( !ALLOW_N ) FetchReadsAmb( assembly_amb, ASSEMBLY, AMB_EQ_n );
          else FetchReadsAmb( assembly_amb, ASSEMBLY, AMB_EQ_Nn );    }
     else
     {    LoadEfastaFlat( ASSEMBLY, assembly );
          LoadEfastaFlatGaps( ASSEMBLY, assembly_amb );    }
     if (PERF_STATS && IsRegularFile( supers_file ))
     {    longlong bases_in_genome = 0;
          vec<superb> supers;
	  ReadSuperbs( supers_file, supers );
          for ( size_t t = 0; t < supers.size( ); t++ )
	    bases_in_genome += supers[t].TrueLength( );
          PerfStat::log( ) << std::fixed << std::setprecision(1)
               << PerfStat( "scaff_per_Mb", 
               "number of scaffolds per Mb",
               1000000.0 * double( assembly.size( ) ) 
               / double(bases_in_genome) );    }
     int nref = ref.size( );
     size_t A = assembly.sumSizes();

     // Select test sequences.

     const int n = 100;
     vecbasevector query;
     // TODO: potentially dangerous truncation of index by all these ints
     vec<int> t1s, t2s, start1s, start2s;
     SelectTestSequences( SAMPLE, D, assembly, assembly_amb, query,
          t1s, t2s, start1s, start2s );
     
     // Align the test sequences.

     vec<look_align> aligns;
     PerfectLookup( 12, query, REFHEAD + ".lookup", aligns, FW_OR_RC );
     vec< vec<int> > aligns_index(2*SAMPLE);
     for ( int i = 0; i < aligns.isize( ); i++ )
          aligns_index[ aligns[i].query_id ].push_back(i);

     // Report granular accuracy results.
     
     int invalids = 0, total = 0;
     for ( int i = 0; i < SAMPLE; i++ )
     {    if ( !aligns_index[2*i].solo( ) ) continue;
          if ( !aligns_index[2*i+1].solo( ) ) continue;
          const look_align& a1 = aligns[ aligns_index[2*i][0] ];
          const look_align& a2 = aligns[ aligns_index[2*i+1][0] ];
          Bool ambig = False;
          const bitvector &r1 = ref_amb[ a1.target_id ]; 
          const bitvector &r2 = ref_amb[ a2.target_id ];
          for ( int j = a1.pos2( ); j < a1.Pos2( ); j++ )
               if ( r1[j] ) ambig = True;
          for ( int j = a2.pos2( ); j < a2.Pos2( ); j++ )
               if ( r2[j] ) ambig = True;
          if (ambig) continue;

          // Make sure we're far from the end.

          /*
          int DP = 2*D;
          int dist_to_end1 
               = Min( a1.pos2( ), ref[ a1.target_id ].isize( ) - a1.Pos2( ) );
          int dist_to_end2 
               = Min( a2.pos2( ), ref[ a2.target_id ].isize( ) - a2.Pos2( ) );
          if ( dist_to_end1 < DP || dist_to_end2 < DP ) continue;
          */

          ++total;
          if (VERBOSE)
          {    cout << "\n#" << i
                    << "\nreference " << a1.target_id << ".";
               if ( a1.Fw1( ) ) cout << a1.pos2( ) << "-" << a1.Pos2( );
               else cout << a2.pos2( ) << "-" << a2.Pos2( );
               cout << ( a1.Rc1( ) ? "[rc]" : "[fw]" )
                    << ", " << a2.target_id << ".";
               if ( a2.Fw1( ) ) cout << a2.pos2( ) << "-" << a2.Pos2( );
               else cout << a1.pos2( ) << "-" << a1.Pos2( );
               cout << ( a2.Rc1( ) ? "[rc]" : "[fw]" ) << "\n"
                    << "assembly " << t1s[i] << "." << start1s[i] << "-"
                    << start1s[i] + n << "[fw], " << t2s[i] << "." << start2s[i] 
                    << "-" << start2s[i] + n << "[fw]" << endl;    }

          // Check for different scaffolds or different orientation.

          if ( a1.target_id != a2.target_id )
          {    if (VERBOSE) cout << "INVALID: different reference scaffolds\n";
               ++invalids;    }
          else if ( a1.Rc1( ) != a2.Rc1( ) ) 
          {    if (VERBOSE) cout << "INVALID: orientation\n";
               ++invalids;    }

          // Same scaffold, consistent orientations.

          else
          {    int delta1, delta2;
               if ( a1.Fw1( ) )
               {    delta1 = a2.pos2( ) - a1.pos2( );
                    delta2 = a2.pos2( ) 
                         + ref[ a1.target_id ].isize( ) - a1.pos2( );    }
               else
               {    delta1 = a1.pos2( ) - a2.pos2( );
                    delta2 = a1.pos2( ) 
                         + ref[ a2.target_id ].isize( ) - a2.pos2( );    }
               int err1 = Abs( delta1 - D ), err2 = Abs( delta2 - D );
               int err = Min( err1, err2 );
               if (VERBOSE) cout << "err = " << err;
               if ( err > MAX_ERR ) 
               {    if (VERBOSE) cout << " INVALID: separation";
                    ++invalids;    }    
               if (VERBOSE) cout << "\n";    }    }    

     // Report summary statistics.

     if ( MIN_CONTIG <= 0 ) MIN_CONTIG = 1;
     if (VERBOSE) cout << "\n";
     
     vec<int> contig_sizes, scaffold_sizes;
     for ( size_t i = 0; i < assembly.size( ); i++ )
     {    int scaffold = 0, count = 0;
          for ( int j = 0; j < assembly[i].isize( ); j++ )
          {    if ( assembly_amb[i][j] && count > 0 )
               {    if ( count >= MIN_CONTIG )
                    {    contig_sizes.push_back(count);
                         scaffold += count;    }
                    count = 0;    }
               if ( !assembly_amb[i][j] ) ++count;    }
          if ( count >= MIN_CONTIG ) 
          {    contig_sizes.push_back(count);
               scaffold += count;    }
          if ( scaffold > 0 ) scaffold_sizes.push_back(scaffold);    }
     Sort(contig_sizes), Sort(scaffold_sizes);
     if ( contig_sizes.size( ) < 1 ) {
       cout << "\nno contig found >= " << MIN_CONTIG << "\nExit." << endl;
       return 0;
     }
     cout << "total bases in contigs = " 
          << ToStringAddCommas( BigSum(contig_sizes) ) << endl;
     cout << "N50 contig size = " << N50(contig_sizes) 
	  << " (" << contig_sizes.isize() << " qualified contigs)" << endl;
     double N50_contig = int( round( double( N50(contig_sizes) ) / 1000.0 ) );
     cout << "N50 scaffold size (excluding gaps) = " << N50(scaffold_sizes) 
	  << " (" << scaffold_sizes.isize() << " scaffolds)" << endl;
     double N50_scaffold = int( round( double( N50(scaffold_sizes) ) / 1000.0 ) );
     PRINT2( invalids, total );
     cout << "validation rate = " << PERCENT_RATIO( 4, total-invalids, total ) 
          << endl;
     if (PERF_STATS)
     {    if ( invalids == 0 )
          {    PerfStat::log( ) << PerfStat( "validation_rate", 
                    "% validation rate at 100 kb", 100 );    }
          else
          {    double validation_rate 
                    = 100.0 * double(total-invalids) / double(total);
               PerfStat::log( ) << std::fixed << std::setprecision(2)
                    << PerfStat( "validation_rate", 
                    "% validation rate at 100 kb", validation_rate );    }    }

     cout << Date() << ": Done!" << endl;
}
