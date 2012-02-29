/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// Happening: show what a running process is doing by causing it to generate
// backtraces.  This only works if the process calls RunTime( ) or the equivalent.
//
// Arguments:
//      PID = process to be traced
//      PROG = name of program to be traced (alternative to PID, Linux only)
//      COUNT = the number of backtraces which are performed (default 1)
//      SEP = the wait time in seconds between backtraces (default 20)
//      DUMP = show raw traceback too
//      TALLY = do all the backtraces and then unique sort with counts.

#include <strstream>

#include <signal.h>

#include "FastIfstream.h"
#include "MainTools.h"

vec<String> tracebacks;
Bool TALLYG;

void Tally( )
{    if (TALLYG)
     {    Sort(tracebacks);
          vec< pair<int,String> > tbcount;
          for ( int i = 0; i < tracebacks.isize( ); i++ )
          {    int j;
               for ( j = i + 1; j < tracebacks.isize( ); j++ )
                    if ( tracebacks[j] != tracebacks[i] ) break;
               tbcount.push_back( make_pair( j - i, tracebacks[i] ) );
               i = j - 1;    }
          ReverseSort(tbcount);
          for ( int i = 0; i < tbcount.isize( ); i++ )
          {    if ( i > 0 ) cout << "\n";
               cout << "===== " << tbcount[i].first << " occurrences =====\n";
               cout << tbcount[i].second;    }    }
     cout << endl;    }

void simple_signal_handler( int signal_number, siginfo_t* info, void* context )
{    if ( signal_number == SIGINT )
     {    Tally( );
          cout << "\nCTRL-C received, exiting" << endl;
          _exit(0);    }
     _exit(1);    }

int main( int argc, char *argv[] )
{
     // RunTime( 1, &simple_signal_handler );
     RunTimeNoTraceback( );

     BeginCommandArguments;
     CommandArgument_Int_OrDefault(PID, -1);
     CommandArgument_String_OrDefault(PROG, "");
     CommandArgument_Bool_OrDefault(DUMP, False);
     CommandArgument_Int_OrDefault(COUNT, 1);
     CommandArgument_UnsignedInt_OrDefault(SEP, 20);
     CommandArgument_Bool_OrDefault(TALLY, False);
     EndCommandArguments;

     /*
     ArachneSignalHandler *pHandlerFunc = &simple_signal_handler;
     ArachneInterruptHandler( pHandlerFunc );
     */

     TALLYG = TALLY;
     ForceAssert( ( PID >= 0 ) ^ ( PROG != "" ) );
     int pid = PID;
     if ( PROG != "" )
     {    temp_file tempfile( "/tmp/tmp_LineOfOutput_XXXXXXX" );
          if ( System( "/sbin/pidof " + PROG + " > " + tempfile ) != 0 )
          {    cout << "Failed: pidof " << PROG << "\n";
               cout << "This couldn't be because your system doesn't have pidof "
                    << "or because " << PROG << " is not running.\n";
               exit(1);    }
          istringstream iline( FirstLineOfFile(tempfile).c_str( ) );
          while(1)
          {    int n;
               iline >> n;
               if ( !iline ) break;
               if ( pid >= 0 ) 
               {    cout << "There is more than one " << PROG << " process.\n";
                    exit(1);    }
               pid = n;    }
          if ( pid < 0 )
          {    cout << "Couldn't find a " << PROG << " process.\n";
               exit(1);    }    }
     String tracefile = "/tmp/traceback_from_process_" + ToString(pid);
     int kill_status = 0;
     for ( int i = 0; i < COUNT; i++ )
     {    if ( i > 0 ) sleep(SEP);
          Remove(tracefile);
          int signal = ( DUMP ? SIGUSR2 : SIGUSR1 );
          kill_status = kill( pid, signal );
          if ( kill_status != 0 ) break;
          while( !IsRegularFile(tracefile) ) usleep(10000);
          String traceback, line;
          {    fast_ifstream in(tracefile);
               while(1)
               {    getline( in, line );
                    if ( in.fail( ) ) break;
                    traceback += line + "\n";    }     }
          Remove(tracefile);
          if ( !TALLY )
          {    if ( i > 0 ) cout << "\n";
               cout << traceback;
               flush(cout);    }
          else tracebacks.push_back(traceback);    }
     Tally( );
     cout << endl;
     if ( kill_status != 0 )
          cout << "Last kill failed.  Perhaps the process finished.\n\n";    }
