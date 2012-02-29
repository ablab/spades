///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// DisplayLocs.  Display read locations on a contig.  MakeReadLocs must have been
// run first.  This is automatic in the RunAllPathsLG pipeline if 
// PATCH_SCAFFOLDS=True.
//
// For SHOW_ALIGNS and SHOW_PILEUPS to work correctly, the bandwidth entries in
// the read locations have to be set.  If we wanted to work around this we could
// pass a default bandwidth to the alignment code.

#include "MainTools.h"
#include "Superb.h"
#include "paths/ReadLoc.h"
#include "system/ParsedArgs.h"

int main(int argc, char *argv[])
{
     RunTime( );
  
     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_String_OrDefault(SUBDIR, "test");
     CommandArgument_String_OrDefault(ASSEMBLY, "linear_scaffolds0.clean.patched");
     CommandArgument_String_OrDefault_Doc(TIGS, "", 
          "if unspecified, process all contigs;"
          " otherwise it is one of the following: \n"
          "(a) a list of contig ids (in ParseIntSet format) or \n"
          "(b) the letter s followed by a list of scaffolds or \n"
          "(c) s<scaffold id>.<list of indices of contigs in the scaffold");
     CommandArgument_Int_OrDefault_Doc(START, -1, 
          "starting point on contig, provided that only one contig specified");
     CommandArgument_Int_OrDefault_Doc(STOP, -1, 
          "stopping point on contig, provided that only one contig specified");
     CommandArgument_Bool_OrDefault(SHOW_LOCS, True);
     CommandArgument_Bool_OrDefault(SHOW_ALIGNS, False);
     CommandArgument_Bool_OrDefault(SHOW_PARTNER_ALIGNS, False);
     CommandArgument_Bool_OrDefault(SHOW_QUALS_IN_ALIGNS, False);
     CommandArgument_Bool_OrDefault(SHOW_PILEUPS, False);
     CommandArgument_Double_OrDefault(PILEUP_PURITY_FILTER, 1.0);
     CommandArgument_String_OrDefault_Doc(READLOCS_PREFIX, "",
          "if specified, file extension is .READLOCS_PREFIX.readlocs "
          "instead of .readlocs");
     CommandArgument_String_OrDefault_Doc(LIBS, "{frag,jump,long}",
          "use reads from these library types");
     EndCommandArguments;

     // Define directories.

     String data_dir = PRE + "/" + DATA, run_dir = PRE + "/" + DATA + "/" + RUN;
     String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR; 

     // Parse LIBS.

     vec<Bool> libs(3, False);
     vec<String> libsx;
     ParseStringSet(LIBS, libsx);
     for ( int j = 0; j < libsx.isize( ); j++ )
     {    if ( libsx[j] == "frag" ) libs[0] = True;
          else if ( libsx[j] == "jump" ) libs[1] = True;
          else if ( libsx[j] == "long" ) libs[2] = True;
          else
          {    cout << "Illegal LIBS argument" << endl;
               exit(1);    }    }

     // Parse TIGS.

     vec<int> tigs;
     {    if (TIGS == "") 
          {    int n_tigs = MastervecFileObjectCount( sub_dir + "/" + ASSEMBLY
                    + ".contigs.fastb" );
               for ( int j = 0; j < n_tigs; j++ )
                    tigs.push_back(j);    }
          else if (TIGS.Contains("s", 0)) 
          {    TIGS = TIGS.After("s");
               vec<superb> scaffolds;
               ReadSuperbs(sub_dir + "/" + ASSEMBLY + ".superb", scaffolds);
               if (TIGS.Contains(".")) 
               {    int scaffold = TIGS.Before(".").Int();
                    ForceAssertLt(scaffold, scaffolds.isize());
                    vec<int> spos;
                    ParseIntSet(TIGS.After("."), spos);
                    for (int j = 0; j < spos.isize(); j++)
                         tigs.push_back(scaffolds[scaffold].Tig(spos[j]));    }
               else 
               {    vec<int> s;
                    ParseIntSet(TIGS, s);
                    for (int i = 0; i < s.isize(); i++) 
                    {    int scaffold = s[i];
                         ForceAssertLt(scaffold, scaffolds.isize());
                         for (int j = 0; j < scaffolds[scaffold].Ntigs(); j++)
                              tigs.push_back(scaffolds[scaffold].Tig(j));
                              }    }    }
          else ParseIntSet(TIGS, tigs);    }
     if ( !tigs.solo( ) && ( START >= 0 || STOP >= 0 ) )
     {    cout << "If you specify START or STOP, you have to specify just "
               << "one contig" << endl;    }

     // Load scaffolds.

     int ntigs 
          = MastervecFileObjectCount( sub_dir + "/" + ASSEMBLY + ".contigs.fastb" );
     vec<superb> scaffolds;
     ReadSuperbs( sub_dir + "/" + ASSEMBLY + ".superb", scaffolds );
     vec<int> to_super( ntigs, -1 ), to_super_pos( ntigs, -1 );
     vec<int> to_super_posr( ntigs, - 1 );
     for ( int i = 0; i < scaffolds.isize( ); i++ )
     {    int n = scaffolds[i].Ntigs( );
          for ( int j = 0; j < n; j++ )
          {    to_super[ scaffolds[i].Tig(j) ] = i;
               to_super_pos[ scaffolds[i].Tig(j) ] = j;   
               to_super_posr[ scaffolds[i].Tig(j) ] = n - j - 1;    }    }

     // Process contigs.

     String head = sub_dir + "/" + ASSEMBLY;
     if ( READLOCS_PREFIX != "" ) head += "." + READLOCS_PREFIX;
     read_locs_on_disk locs_file( head, run_dir );
     for ( int it = 0; it < tigs.isize( ); it++ )
     {    int TIG = tigs[it];
          vec<read_loc> locs;
          locs_file.LoadContig( TIG, locs );
          vec<Bool> to_delete( locs.size( ), False );
          for ( int j = 0; j < locs.isize( ); j++ )
          {    if ( !libs[ locs[j].ReadClass( ) ] ) to_delete[j] = True;
               if ( START >= 0 && locs[j].Stop( ) <= START ) to_delete[j] = True;
               if ( STOP >= 0 && locs[j].Start( ) >= STOP ) to_delete[j] = True;    }
          EraseIf( locs, to_delete );
	  if ( SHOW_LOCS )
          {    PrintWithPairInfo( cout, locs, scaffolds, to_super, to_super_pos, 
	            to_super_posr, run_dir, sub_dir, ASSEMBLY, SHOW_ALIGNS, 
		    SHOW_PARTNER_ALIGNS, SHOW_QUALS_IN_ALIGNS );    }
          if (SHOW_PILEUPS)
          {    vecbasevector tig;
               tig.ReadOne( sub_dir + "/" + ASSEMBLY + ".contigs.fastb", TIG );
               vec<dumbcall> calls;
               Pileup( tig[0], locs, run_dir, calls );
               cout << "\nCOVERAGE OF CONTIG " << TIG
                    << " (" << tig[0].size( ) << " BASES)\n";
               int start = ( START >= 0 ? START : 0 );
               int stop = ( STOP >= 0 ? STOP : tig[0].size( ) );
               for ( int j = start; j < stop; j++ )
               {    int ncalls = 0, ncalls_agree = 0;
                    for ( int k = 0; k < 6; k++ )
                    {    ncalls += calls[j].base[k];
                         if ( k == tig[0][j] ) ncalls_agree += calls[j].base[k];    }
                    double agree = double(ncalls_agree)/double(ncalls);
                    if ( ( ncalls == 0 && PILEUP_PURITY_FILTER == 1.0 ) 
                         || agree <= PILEUP_PURITY_FILTER )
                    {    cout << j << " ";
                         for ( int k = 0; k < 6; k++ )
                         {    for ( int l = 0; l < calls[j].base[k]; l++ )
                              {    if ( k < 4 ) cout << as_base(k);
                                   if ( k == 4 ) cout << "D";
                                   if ( k == 5 ) cout << "I";    }    }
                         if ( ncalls > 0 )
                         {    cout << " (" << as_base( tig[0][j] ) << " "
                                   << PERCENT_RATIO(3, ncalls_agree, ncalls)
                                   << ")";    }
                         cout << "\n";    }    }    }    }    }
