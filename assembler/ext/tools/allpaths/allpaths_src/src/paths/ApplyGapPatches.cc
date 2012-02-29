///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// ApplyGapPatches.  Read multiple suggested gaps patches of an assembly and 
// perform the patches.  This code uses the particular properties of the different
// patch programs to try to make the best choice.

#include "Basevector.h"
#include "Fastavector.h"
#include "FetchReads.h"
#include "MainTools.h"
#include "Superb.h"
#include "efasta/EfastaTools.h"
#include "math/Functions.h"
#include "paths/AssemblyCleanupTools.h"
#include "paths/AssemblyEdit.h"
#include "paths/PostPatcherBridgeGap.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_String_OrDefault(SUBDIR, "test");
     CommandArgument_String_OrDefault(SCAFFOLDS_IN, "linear_scaffolds0.clean");
     CommandArgument_String_OrDefault(SCAFFOLDS_OUT, 
          "linear_scaffolds0.clean.applied");
     CommandArgument_StringSet_OrDefault(EDITS_IN,
          "linear_scaffolds0.clean.patched.edits, "
          "linear_scaffolds0.clean.longread.fixed.edits, "
          "linear_scaffolds0.clean.kpatch.edits");
     CommandArgument_Bool_OrDefault(VERBOSE, False);
     CommandArgument_Double_OrDefault_Doc(BRIDGEGAP_MAX_SIGMA, 5.0, "If positive, "
          "max standard deviations a bridge may vary from estimated gap size");
     CommandArgument_Int_OrDefault_Doc(BRIDGEGAP_FLAGS, BRIDGEGAP_BESTFIT,
	  "1: ONESIZE, 2: BESTFIT, 4: MAXCOUNT");
     CommandArgument_Bool_OrDefault(WRITE, True);
     EndCommandArguments;

     // Begin.

     double clock = WallClockTime();

     // Define directories, etc.

     String data_dir = PRE + "/" + DATA, run_dir = PRE + "/" + DATA + "/" + RUN;
     String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;
     String out_head = sub_dir + "/" + SCAFFOLDS_OUT;
     String tigsa_file = sub_dir + "/" + SCAFFOLDS_IN + ".contigs.fasta";
     String efasta_file = sub_dir + "/" + SCAFFOLDS_IN + ".contigs.efasta";
     String supers_file = sub_dir + "/" + SCAFFOLDS_IN + ".superb";

     // Read assembly.

     vec<fastavector> tigsa;
     LoadFromFastaFile( tigsa_file, tigsa );
     int n_tigs = tigsa.size( );
     vec<efasta> tigse;
     LoadEfastaIntoStrings( efasta_file, tigse );
     vec<superb> scaffolds;
     ReadSuperbs( supers_file, scaffolds );

     // Read the edits.

     cout << Date( ) << ": Importing patches from:" << endl;
     vec< vec< vec<assembly_edit> > > edits( EDITS_IN.size( ) );
     for ( size_t i = 0; i < EDITS_IN.size( ); i++ )
     {    edits[i].resize(n_tigs);
          String filename =  sub_dir + "/" + EDITS_IN[i];
	  cout << filename << endl;
          vec<assembly_edit> editsi;
          BinaryReader::readFile( filename.c_str( ), &editsi );
          for ( size_t j = 0; j < editsi.size( ); j++ )
               edits[i][ editsi[j].Tig1( ) ].push_back( editsi[j] );    }
     int patched_index = -1, longread_index = -1, kpatch_index = -1;
     for ( size_t i = 0; i < EDITS_IN.size( ); i++ )
     {    if ( EDITS_IN[i].Contains( ".patched.edits", -1 ) ) patched_index = i;
          if ( EDITS_IN[i].Contains( ".longread.fixed.edits", -1 ) )
               longread_index = i;
          if ( EDITS_IN[i].Contains( ".kpatch.edits", -1 ) ) kpatch_index = i;    }

     // Install the patches, working backwards through each scaffold so that edits
     // don't step on each other.  This code is quadratic in the scaffold size, 
     // since each edit causes the entire scaffold to be edited.  

     for ( size_t s = 0; s < scaffolds.size( ); s++ )
     {    superb& S = scaffolds[s];
          int last_m1 = -1, last_start1 = -1;
          for ( int p = S.Ngaps( ) - 1; p >= 0; p-- )
          {    int m1 = S.Tig(p), m2 = S.Tig(p+1);
               Bool use_patched
                    = patched_index >= 0 && edits[patched_index][m1].nonempty( );
               Bool use_longread
                    = longread_index >= 0 && edits[longread_index][m1].nonempty( );
               Bool use_kpatch
                    = kpatch_index >= 0 && edits[kpatch_index][m1].nonempty( );
               if ( VERBOSE && ( use_patched || use_longread || use_kpatch ) )
               {    cout << "\n===================================================="
                         << "================================\n";
                    cout << "\nconsidering patches from " << m1 
                         << "[l=" << tigsa[m1].size( ) << "] to " << m2 
                         << "[l=" << tigsa[m2].size( ) << "], gap = "
                         << S.Gap(p) << " +/- " << S.Dev(p) << "\n";    }
               if ( VERBOSE && use_patched ) cout << "have patch from PostPatcher\n";
               if ( VERBOSE && use_longread )
                    cout << "have patch from LongReadPostPatcher\n";
               if ( VERBOSE && use_kpatch ) cout << "have patch from KPatch\n";
               fastavector b;
               Bool b_defined = False;

               // If we're using longreads, KPatch trumps.  Otherwise PostPatcher
               // trumps it.

               if ( longread_index >= 0 && use_kpatch )
               {    use_patched = False;
                    use_longread = False;    }
               if ( longread_index < 0 && use_kpatch && use_patched )
                    use_kpatch = False;

               // Decide if PostPatcher patch is OK.  Note that PostPatcher could
               // have contributed two patches, first from PostPatcher per se, and
               // second from its MERGE_CONSECUTIVE_CONTIGS phase.  If the first one
               // fails, we use the second one (if it exists).

               int pid = 0;
               if (use_patched)
               {    const assembly_edit& e = edits[patched_index][m1][0];
                    vec<int> gaps;
                    for ( int j = 0; j < e.Nreps( ); j++ )
                    {    gaps.push_back( e.Start1( ) + e.Rep(j).isize( )
                              - e.Stop2( ) - (int) tigsa[m1].size( ) );    }
	            if ( BridgeGap( e.Reps( ), gaps, S.Gap(p), S.Dev(p), 
                         BRIDGEGAP_MAX_SIGMA, BRIDGEGAP_FLAGS, VERBOSE, b ) ) 
                    {    b_defined = True;    }
                    else
                    {    if ( edits[patched_index][m1].size( ) == 1 )
                              use_patched = False;    
                         else pid = 1;    }    }

               // If the gap is closed by both PostPatcher and LongReadPostPatcher,
               // pick one.

               if ( use_patched && use_longread )
               {    const assembly_edit& e = edits[patched_index][m1][pid];
                    vec<double> abs_offby( e.Nreps( ) );
                    for ( int j = 0; j < e.Nreps( ); j++ )
                    {    int gap = e.Start1( ) + e.Rep(j).isize( )
                              - e.Stop2( ) - (int) tigsa[m1].size( );
                         abs_offby[j] = Abs( double( gap - S.Gap(p) ) 
                              / double( S.Dev(p) ) );    }
                    Sort(abs_offby);
                    const double max_offby = 4.0;
                    const double min_dist = 2.0;
                    if ( abs_offby[0] > max_offby || ( abs_offby.size( ) > 1 
                         && abs_offby[1] - abs_offby[0] <= min_dist ) )
                    {    use_patched = False;
                         b_defined = False;    }
                    else use_longread = False;    }

               // Patches from PostPatcher and LongReadPostPatcher use fasta
               // coordinates.  Translate them to efasta coordinates.

               for ( int pass = 1; pass <= 2; pass++ )
               {    int index = ( pass == 1 ? patched_index : longread_index );
                    if ( index < 0 ) continue;
                    for ( int j = 0; j < edits[index][m1].isize( ); j++ )
                    {    assembly_edit& e = edits[index][m1][j];
                         int start1 = e.Start1( ), stop2 = e.Stop2( );
                         String s1 = tigsa[m1].ToString( ).substr( 
                              start1, (int) tigsa[m1].size( ) - start1 );
                         String sexp1 = ExpandAmbCode(s1);
                         e.SetStart1( tigse[m1].isize( ) - sexp1.isize( ) );
                         String s2 = tigsa[m2].ToString( ).substr( 0, stop2 );
                         e.SetStop2( ExpandAmbCode(s2).size( ) );    }    }

               // Install the patch;

               if ( use_patched || use_longread || use_kpatch )
               {    fastavector M1 = tigsa[m1], M2 = tigsa[m2];
                    efasta E1 = tigse[m1], E2 = tigse[m2];
                    const assembly_edit& e = ( use_patched ? 
                         edits[patched_index][m1][pid] :
                         ( use_kpatch ? edits[kpatch_index][m1][0] :
                         edits[longread_index][m1][0] ) );
                    if ( !b_defined && !use_kpatch ) b = fastavector( e.Rep(0) );
                    efasta evec;
                    if ( !use_kpatch ) evec = efasta(b);
                    else evec = efasta( e.Reps( ) );
                    if (VERBOSE) 
                    {    cout << "patch has effective size "
                              << evec.size( ) - ( E1.size( ) - e.Start1( ) ) 
                              - e.Stop2( ) << endl;    }
               
                    // Check for collision with last patch.

                    if ( m2 == last_m1 && e.Stop2( ) > last_start1 )
                    {    if (VERBOSE)
                         {    cout << "collision, can't use patch from " << m1 
                                   << " to " << m2 << "\n";    }
                         continue;    }
                    last_m1 = m1;
                    last_start1 = e.Start1( );

                    // Go ahead and do it.

                    ForceAssertLe( e.Start1( ), (int) E1.size( ) );
                    ForceAssertLe( e.Stop2( ), (int) E2.size( ) );
  	            E1.resize( e.Start1( ) );
                    E2 = E2.substr( e.Stop2( ), E2.size( ) - e.Stop2( ) );
	            efasta E1B = E1 + evec;
	            efasta newtig = E1B + E2;
	            tigse[m1] = newtig, tigse[m2].resize(0);
	            S.SetLen( p, newtig.Length1( ) );
                    if (VERBOSE) 
                    {    cout << "new contig has length " << newtig.Length1( ) 
                              << "\n";    }
	            int gap = 0, dev = 0;
	            bool nextgap = ( p+1 < S.Ngaps( ) );
	            if (nextgap) 
                    {    gap = S.Gap(p+1), dev = S.Dev(p+1);    }
	            S.RemoveTigByPos(p+1);
	            if (nextgap) 
                    {    S.SetGap( p, gap ), S.SetDev( p, dev );    }    }    }    }

     // Write results.

     if (WRITE)
     {    cout << Date( ) << ": writing assembly" << endl;
          Assembly A( scaffolds, tigse );
          A.remove_unused_contigs( );
          A.check_integrity( );
          A.WriteAll( sub_dir + "/" + SCAFFOLDS_OUT );    }
     cout << Date( ) << ": done, time used = " << TimeSince(clock) << endl;    }
