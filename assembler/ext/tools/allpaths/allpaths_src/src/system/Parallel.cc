/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include <sys/wait.h>

#include "CoreTools.h"
#include "Set.h"

void RunCommandsInParallel( const vec<String>& commands,
			    const int max_threads,
			    ostream *pOut )
{    int pcount = 0, pcount_active = 0, status;
     set<int> pids;
     ofstream devnull ( "/dev/null" );
     ostream &out = pOut ? *pOut : devnull;
     for ( int i = 0; i < commands.isize( ); i++ )
     {    pcount_active++;
          pcount++;
	  out << Date( ) << ": starting batch_" << i << endl;
          int pid = Fork( commands[i] );
          pids.insert(pid);
          if ( pcount >= max_threads )
          {    int pid2;
               while(1)
               {    pid2 = wait( &status );
                    if ( pid2 == -1 )
                    {    cout << "Wait failed." << endl;
                         cout << "Abort." << endl;
                         exit(1);    }
                    if ( Member( pids, pid2 ) ) {
		      set<int>::iterator it = find( pids.begin( ), pids.end( ),
						    pid2 );
		      if ( it != pids.end( ) )
			out << Date( ) << ": batch_"
			    << distance( pids.begin( ), it ) << " done" << endl;
		      break;
		    }
	       }
               pcount_active--;
                if ( status != 0 ) 
                {
                    String err_msg = String("Process exited with nonzero status: ") + commands[0];
                    FatalErr(err_msg.c_str());
                }
            }   
        }
     for ( int i = 0; i < pcount_active; i++ )
     {    int pid2;
          while(1)
          {    pid2 = wait( &status );
               if ( pid2 == -1 )
               {    cout << "Wait failed." << endl;
                    cout << "Abort." << endl;
                    exit(1);    }
               if ( Member( pids, pid2 ) ) {
		 set<int>::iterator it = find( pids.begin( ), pids.end( ),
					       pid2 );
		 if ( it != pids.end( ) )
		   out << Date( ) << ": batch_"
		       << distance( pids.begin( ), it ) << " done" << endl;
		 break;
	       }
	  }
          if ( status != 0 ) FatalErr( "process failed" );    }    }
