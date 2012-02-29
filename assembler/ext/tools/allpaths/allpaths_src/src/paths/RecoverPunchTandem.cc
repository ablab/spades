///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// RecoverPunchTandem.  This is run in sequence as PunchTandemHoles,
// LongReadPostPatcher, RecoverPunchTandem.  It reads in the cleaned long read 
// patches and either stuffs them back into the assembly or else uses the original 
// content that PunchTandemHoles saved.

#include "Basevector.h"
#include "Fastavector.h"
#include "MainTools.h"
#include "PackAlign.h"
#include "Superb.h"
#include "efasta/EfastaTools.h"
#include "pairwise_aligners/SmithWatFree.h"
#include "paths/AssemblyCleanupTools.h"
#include "paths/AssemblyEdit.h"

int SmithWatFreeSymPenalize( const basevector& b1, const basevector& b2, align& a )
{    alignment al;
     int best_loc, errs;
     if ( b1.size( ) <= b2.size( ) )
     {    errs = SmithWatFree( b1, b2, best_loc, al, true, true, 1, 1 );
          a.UnpackFrom(al);    }
     else
     {    errs = SmithWatFree( b2, b1, best_loc, al, true, true, 1, 1 );
          a.UnpackFrom(al);
          a.Flip( );    }    
     return errs;    }

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(SCAFFOLDS_IN);
     CommandArgument_String_OrDefault(SCAFFOLDS_OUT, SCAFFOLDS_IN + ".recover");
     CommandArgument_Bool_OrDefault(VERBOSE, False);
     CommandArgument_Bool_OrDefault(WRITE, True);
     EndCommandArguments;

     // Load assembly.

     cout << Date( ) << ": loading data" << endl;
     String supers_file = SCAFFOLDS_IN + ".superb";
     String tigs_file = SCAFFOLDS_IN + ".contigs.fastb";
     String efasta_file = SCAFFOLDS_IN + ".contigs.efasta";
     vec<superb> scaffolds;
     ReadSuperbs( supers_file, scaffolds );
     vecbasevector tigs(tigs_file);
     int ntigs = tigs.size( );
     vec<efasta> tigse;
     LoadEfastaIntoStrings( efasta_file, tigse );

     // Load edits.

     cout << Date( ) << ": loading edits" << endl;
     vec< vec< vec<assembly_edit> > > edits(2);
     for ( int i = 0; i < 2; i++ )
     {    edits[i].resize(ntigs);
          String filename = SCAFFOLDS_IN;
          if ( i == 0 ) filename += ".edits";
          else filename += ".longread.edits";
          vec<assembly_edit> editsi;
          BinaryReader::readFile( filename.c_str( ), &editsi );
          for ( size_t j = 0; j < editsi.size( ); j++ ) 
            if (editsi[j].Rep(0).size() > 0) 
              edits[i][ editsi[j].Tig1( ) ].push_back( editsi[j] );    }

     // Install the patches, working backwards through each scaffold so that edits
     // don't step on each other.  This code is quadratic in the scaffold size, 
     // since each edit causes the entire scaffold to be edited.  

     cout << Date( ) << ": installing patches" << endl;
     for ( size_t s = 0; s < scaffolds.size( ); s++ )
     {    superb& S = scaffolds[s];
          int last_m1 = -1, last_start1 = -1;
          for ( int p = S.Ntigs( ) - 1; p >= 0; p-- )
          {    int m1 = S.Tig(p), m2 = S.Tig(p+1);
               if ( edits[0][m1].empty( ) ) continue;
               if (VERBOSE)
               {    cout << "\nm1 = " << m1 << ", gap = " 
                         << S.Gap(p) << " +/- " << S.Dev(p) << endl;
                    PRINT( efasta( edits[0][m1][0].Reps( ) ) );
                    cout << "number of long read patches = "
                         << edits[1][m1].size( ) << endl;
                    for ( int j = 0; j < edits[1][m1].isize( ); j++ )
                         PRINT2( j, edits[1][m1][j].Reps( ) );    }

               // If there is no long read patch, we use the original patch.

               assembly_edit ed;
               if ( edits[1][m1].empty( ) ) ed = edits[0][m1][0];

               // Otherwise, we take the original patch and expand or contract the
               // tandem repeat until it most closely matches the long read patch.
               // We call that the answer.

               else
               {    basevector lr_patch = edits[1][m1][0].Rep(0);
                    vec<basevector> bx = edits[0][m1][0].Reps( );

                    // Test for markup.

                    if ( bx.solo( ) ) ed = edits[1][m1][0];

                    // Tandem repeat case.

                    else
                    {    efasta e(bx);
                         int start1 = edits[1][m1][0].Start1( );
                         int stop2 = edits[1][m1][0].Stop2( );
                         basevector left_ext( tigs[m1], start1, -1 );
                         if (VERBOSE)
                         {    PRINT2( start1, tigs[m1].size( ) );
                              PRINT2( left_ext, lr_patch );    }
                         basevector right_ext( tigs[m1], 0, stop2 );
                         String left = left_ext.ToString( ) + e.Before( "{" ); 
                         String right = e.After( "}" ) + right_ext.ToString( );
                         efasta middle = "{" + e.Between( "{", "}" ) + "}";
                         vec<basevector> b;
                         middle.ExpandTo(b);
                         sort( b.begin( ), b.end( ), LtBySize<basevector>( ) );
                         String repeat, b1 = b[1].ToString( );
                         for ( int div = b[1].isize( ); div >= 1; div-- )
                         {    if ( b1.isize( ) % div != 0 ) continue;
                              int n = b1.isize( ) / div;
                              repeat = b1.substr( 0, n );
                              String whole;
                              for ( int j = 0; j < div; j++ )
                                   whole += repeat;
                              if ( whole == b1 ) break;    }
                         while ( left.Contains( repeat, -1 ) )
                         {    left.resize( left.isize( ) - repeat.isize( ) );    }
                         while ( right.Contains( repeat, 0 ) )
                         {    right = right.substr( repeat.isize( ),
                                   right.isize( ) - repeat.isize( ) );   }
                         vec<int> errs;
                         String nrepeat;
                         while(1)
                         {    basevector b( left + nrepeat + right );
                              align a;
                              int e = SmithWatFreeSymPenalize( b, lr_patch, a );
                              if (VERBOSE) PRINT3( b, lr_patch, e );
                              if ( errs.nonempty( ) && e > errs.back( ) ) break;
                              errs.push_back(e);
                              nrepeat += repeat;    }
                         nrepeat.resize( nrepeat.size( ) - repeat.size( ) );
                         ed = edits[1][m1][0];
                         ed.Rep(0) = basevector( left + nrepeat + right );    
                         if (VERBOSE) PRINT( ed.Rep(0) );    }    }

               // Check for collision with last patch.

               if ( m2 == last_m1 && ed.Stop2( ) > last_start1 )
               {    if (VERBOSE)
                    {    cout << "collision, can't use patch from " << m1 
                              << " to " << m2 << "\n";    }
                    continue;    }
               last_m1 = m1;
               last_start1 = ed.Start1( );

               // Install the patch.  This is complicated by the fact that
               // LongReadPostPatcher isn't written to handle efasta.  So we have
               // to be excruciatingly careful not to accidentally create nonsense.
               // This involves backing up over brackets.

               efasta M1 = tigse[m1], M2 = tigse[m2];
               int e1 = M1.Index1Alt( ed.Start1() ), e2 = M2.Index1Alt( ed.Stop2() );
               int last_left = -1, brack_count = 0;
               String left_extra, right_extra;
               for ( int j = 0; j < e1; j++ )
               {    if ( M1[j] == '{' ) 
                    {    brack_count++;
                         last_left = j;    }
                    if ( M1[j] == '}' ) brack_count--;    }
               if ( brack_count != 0 )
               {    left_extra = M1.substr( last_left+1, e1 - (last_left+1) );
                    e1 = last_left;    }
               brack_count = 0;
               for ( int j = 0; j < e2; j++ )
               {    if ( M2[j] == '{' ) 
                    {    brack_count++;
                         last_left = j;    }
                    if ( M2[j] == '}' ) brack_count--;    }
               if ( brack_count != 0 )
               {    int comma_pos = -1, rbrack_pos = -1;
                    for ( int j = last_left; j < M2.isize( ); j++ )
                    {    if ( comma_pos < 0 && M2[j] == ',' ) comma_pos = j;
                         if ( M2[j] == '}' )
                         {    rbrack_pos = j;
                              break;    }    }
                    ForceAssertGe( comma_pos, 0 );
                    ForceAssertGe( rbrack_pos, 0 );
                    right_extra = M2.substr( e2, comma_pos - e2 );
                    e2 = rbrack_pos + 1;    }
               if (VERBOSE)
               {    PRINT3( M1.size( ), M2.size( ), tigs[m1].size( ) );
                    PRINT2( ed.Start1( ), ed.Stop2( ) );
                    cout << "trimming off ";
                    for ( int j = e1; j < M1.isize( ); j++ )
                         cout << M1[j];
                    cout << "\n";    }
               M1.resize(e1);
               M1 += left_extra;
               M2 = M2.substr( e2, M2.isize( ) - e2 );
               M2 = right_extra + M2;
               efasta newtig = M1 + efasta( ed.Reps( ) ) + M2;
               tigse[m1] = newtig, tigse[m2].resize(0);
               S.SetLen( p, newtig.Length1( ) );
               if (VERBOSE) 
                    cout << "new contig has length " << newtig.size( ) << "\n";
               int gap = 0, dev = 0;
               bool nextgap = ( p+1 < S.Ngaps( ) );
               if (nextgap) 
               {    gap = S.Gap(p+1), dev = S.Dev(p+1);    }
               S.RemoveTigByPos(p+1);
               if (nextgap) 
               {    S.SetGap( p, gap ), S.SetDev( p, dev );    }    }    }

     // Write results.

     if (WRITE)
     {    cout << Date( ) << ": writing assembly" << endl;
          Assembly A( scaffolds, tigse );
          A.remove_unused_contigs( );
          A.check_integrity( );
          A.WriteAll( SCAFFOLDS_OUT );    }    }
