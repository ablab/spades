///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// EvaluateGaps.  Use reference to assess accuracy of gaps.  
//
// Don't run multiple instances on the same sub_dir.  To allow this we would need
// to use better tmp file names.
//
// This code will run much faster the second time, after the genome lookup table 
// has been cached in memory.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS
// MakeDepend: dependency QueryLookupTable

#include <omp.h>

#include "Basevector.h"
#include "FastIfstream.h"
#include "MainTools.h"
#include "Superb.h"
#include "lookup/LookAlign.h"
#include "math/Functions.h"
#include "random/Random.h"
#include "reporting/PerfStat.h"

int main(int argc, char **argv)
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_String_OrDefault(SUBDIR, "test");
     CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
	   "Number of threads to use (use all available processors if set to 0)");
     CommandArgument_String_OrDefault(SCAFFOLDS_IN, "linear_scaffolds0.patched");
     CommandArgument_Int_OrDefault_Doc(NGAPS, 1000, 
          "upper limit on number of gaps to be assessed");
     CommandArgument_Bool_OrDefault(VERBOSE, False);
     EndCommandArguments;

     // Check arguments.

     ForceAssertGe( NGAPS, 1 );

     // Define directories.

     String data_dir = PRE + "/" + DATA, run_dir = PRE + "/" + DATA + "/" + RUN;
     String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;
     String tmp_dir = sub_dir + "/tmp";
     Mkdir777(tmp_dir);

     // Thread control

     NUM_THREADS = configNumThreads(NUM_THREADS);
     omp_set_num_threads( NUM_THREADS );

     // Load scaffolds.

     vec<superb> S;
     ReadSuperbs( sub_dir + "/" + SCAFFOLDS_IN + ".superb", S );
     int ngaps = 0;
     for ( int i = 0; i < S.isize( ); i++ )
          ngaps += S[i].Ngaps( );

     // Get basic stats.

     vec<int> gaps, devs;
     for ( int i = 0; i < S.isize( ); i++ )
     {    for ( int j = 0; j < S[i].Ngaps( ); j++ )
          {    gaps.push_back( S[i].Gap(j) );
               devs.push_back( S[i].Dev(j) );    }    }
     Sort(gaps), Sort(devs);
     int median_gap = 0, median_dev = 0;
     if ( gaps.nonempty( ) )
     {    median_gap = gaps[ngaps/2], median_dev = devs[ngaps/2];    }

     // Select gaps for alignment.

     vec<Bool> use( ngaps, True );
     int gcount = ngaps;
     while( gcount > NGAPS )
     {    int x = randomx( ) % ngaps;
          if ( use[x] )
          {    use[x] = False;
               gcount--;    }    }
     size_t ntigs = MastervecFileObjectCount( 
          sub_dir + "/" + SCAFFOLDS_IN + ".contigs.fastb" );
     vec<Bool> tigs_to_use( ntigs, False );
     gcount = 0;
     for ( int i = 0; i < S.isize( ); i++ )
     {    for ( int j = 0; j < S[i].Ngaps( ); j++ )
          {    if ( use[gcount] )
               {    tigs_to_use[ S[i].Tig(j) ] = True;
                    tigs_to_use[ S[i].Tig(j+1) ] = True;    }
               gcount++;    }    }

     // Align the contigs.

     for ( size_t j = 0; j < NUM_THREADS; j++ )
     {    Ofstream( iout, tmp_dir + "/AssessGaps.ind" + ToString(j) );
          int count = 0;
          for ( size_t i = 0; i < ntigs; i++ )
          {    if ( tigs_to_use[i] ) 
               {    if ( count % NUM_THREADS == j ) iout << i << "\n";
                    count++;    }    }    }
     #pragma omp parallel for
     for ( size_t j = 0; j < NUM_THREADS; j++ )
     {    SystemSucceed( "QueryLookupTable K=12 MM=12 MC=0.15 "
			 " TMP_DIR=" + tmp_dir +" SEQS=" + sub_dir + "/" 
			 + SCAFFOLDS_IN + ".contigs.fastb" + " L=" + data_dir 
			 + "/genome.lookup PARSEABLE=True SEQS_TO_PROCESS=@"
			 + tmp_dir + "/AssessGaps.ind" + ToString(j)+ " > " + tmp_dir 
			 + "/AssessGaps.aligns" + ToString(j) );    }

     // Load the alignments.

     vec<look_align> aligns;
     vec< vec<int> > aligns_index(ntigs);
     String line;
     for ( size_t j = 0; j < NUM_THREADS; j++ )
     {    fast_ifstream in( tmp_dir + "/AssessGaps.aligns" + ToString(j) );
          while(1)
          {    getline( in, line );
               if ( in.fail( ) ) break;
               if ( !line.Contains( "QUERY", 0 ) ) continue;
               look_align la;
               la.ReadParseable(line);
               aligns_index[ la.query_id ].push_back( aligns.size( ) );
               aligns.push_back(la);    }
          Remove( tmp_dir + "/AssessGaps.ind" + ToString(j) );
          Remove( tmp_dir + "/AssessGaps.aligns" + ToString(j) );    }

     // Survey the gaps.

     int gaps_assayed = 0, gap3 = 0, gap4 = 0, gap5 = 0;
     gcount = -1;
     for ( int i = 0; i < S.isize( ); i++ )
     {    for ( int j = 0; j < S[i].Ngaps( ); j++ )
          {    gcount++;
               if ( !use[gcount] ) continue;
               int m1 = S[i].Tig(j), m2 = S[i].Tig(j+1);
               if ( !aligns_index[m1].solo( ) || !aligns_index[m2].solo( ) )
                    continue;
               const look_align& la1 = aligns[ aligns_index[m1][0] ];
               const look_align& la2 = aligns[ aligns_index[m2][0] ];
               if ( la1.target_id != la2.target_id || la1.Fw1( ) != la2.Fw1( ) )
                    continue;
               int actual_gap;
               if ( la1.Fw1( ) )
               {    if ( la1.Pos1( ) != (int) la1.query_length ) continue;
                    if ( la2.pos1( ) != 0 ) continue;
                    actual_gap = la2.pos2( ) - la1.Pos2( );    }
               else
               {    if ( la2.Pos1( ) != (int) la2.query_length ) continue;
                    if ( la1.pos1( ) != 0 ) continue;
                    actual_gap = la1.pos2( ) - la2.Pos2( );    }
               int gap = S[i].Gap(j), dev = S[i].Dev(j);
               double offby = double( gap - actual_gap ) / double(dev);
	       String str_offby = ToString( offby, 2 );
	       if ( Abs(offby) >= 5.0 ) {
		 ++gap5;
		 str_offby += "  #####";
	       }
               if (VERBOSE) PRINT6( m1, m2, actual_gap, gap, dev, str_offby );
	       if ( Abs(offby) >= 4.0 ) ++gap4;
	       if ( Abs(offby) >= 3.0 ) ++gap3;

               gaps_assayed++;    }    }

     // Summarize results.

     cout << "median gap = " << median_gap << endl;
     cout << "median dev = " << median_dev << endl;
     double gap3_rate 
          = 100.0 * ( gaps_assayed == 0 ? 0 : double(gap3)/double(gaps_assayed) );
     double gap4_rate 
          = 100.0 * ( gaps_assayed == 0 ? 0 : double(gap4)/double(gaps_assayed) );
     double gap5_rate 
          = 100.0 * ( gaps_assayed == 0 ? 0 : double(gap5)/double(gaps_assayed) );
     cout << "evaluated gaps = " << gaps_assayed << endl;
     cout << "evaluated gaps off from real value by at least 3 devs = " << gap3 
          << " (" << fixed << setprecision(1) << gap3_rate << "%)" << endl;
     cout << "evaluated gaps off from real value by at least 4 devs = " << gap4 
          << " (" << fixed << setprecision(1) << gap4_rate << "%)" << endl;
     cout << "evaluated gaps off from real value by at least 5 devs = " << gap5 
          << " (" << fixed << setprecision(1) << gap5_rate << "%)\n" << endl;

     if ( SCAFFOLDS_IN == "linear_scaffolds0" ) 
         PerfStat::log( ) << std::fixed << std::setprecision(1) << PerfStat( 
          "gap5_rate_init", "percent of init gaps off by 5+ devs", gap5_rate );    
     else
         PerfStat::log( ) << std::fixed << std::setprecision(1) << PerfStat( 
          "gap5_rate", "percent of gaps off by 5+ devs", gap5_rate );    }
