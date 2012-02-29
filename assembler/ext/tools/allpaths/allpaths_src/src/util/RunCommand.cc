/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "system/System.h"
#include "system/HostName.h"
#include "util/RunCommand.h"
#include <sys/wait.h>
#include <unistd.h>

/**
 * RunCommand
 */
int RunCommand( const String &the_command )
{
  int exit_code = System( the_command );
  if ( ( WIFEXITED( exit_code ) && WEXITSTATUS( exit_code ) != 0 ) ||
       WIFSIGNALED( exit_code ) )
    FatalErr( the_command << " failed." );
  
  return exit_code;
}

/**
 * RunCommandWithLog
 */
int RunCommandWithLog( const String &the_command,
		       const String &log,
		       const bool append )
{
  if ( IsRegularFile( log ) && ! append ) Remove( log );
  String command = the_command + " 2>> " + log + " 1>> " + log;
  return RunCommand( command );
}

/**
 * RunCommandBool
 */
bool RunCommandBool( const String &the_command )
{
  int exit_code = System( the_command );
  if ( ( WIFEXITED( exit_code ) && WEXITSTATUS( exit_code ) != 0 ) ||
       WIFSIGNALED( exit_code ) )
    return false;

  return true;
}

/**
 * LocateCommand
 */
int LocateCommand( const vec<String>& commands, String key )
{
  String key_orig = key;
  int key_instance = 1, si;
  for ( si = key.isize( ) - 1; si >= 1; si-- )
    if ( !isdigit( key[si] ) ) break;
  ++si;
  if ( si >= 1 && si != key.isize( ) && key[si-1] == '.' ) {
    String begin = key.substr( si, key.size( ) - si + 1 );
    key_instance = begin.Int( );
    key = key.substr( 0, si-1 );
  }
  for ( int i = 0; i < commands.isize( ); i++ )
    if ( commands[i].Contains( key, 0 ) && --key_instance == 0 )
      return i;
  FatalErr( "Can't find " << key_orig << ".\n" );
  return 0;
}

// Print out header for an Assembly run.
void AssemblyHeader( const String &name, int argc, char **argv, ostream &out )
{
  String bar;
  while ( bar.size( ) < 80 ) bar += "=";

  // Original command

  String origcommand;
  for (int ii=0; ii<argc; ii++)
    origcommand += String( " " ) + argv[ii];
  
  // user name

  String start_date = Date( );
    char* username_char = getenv( "LOGNAME" );
  String username( "unknown" );
  if ( username_char != 0 ) username = username_char;

  // Compute datasize limit.
  
  String datasize_limit_string;
  {    ostrstream dls;
  struct rlimit rl;
  if ( getrlimit( RLIMIT_DATA, &rl ) != 0 ) FatalErr( "getrlimit failed" );
  double data_limit = rl.rlim_cur;
  double million = 1000000.0, billion = 1000.0 * million;
  if ( data_limit < million )
    dls << setprecision(3) << data_limit/million << " Mb";
  else dls << setprecision(3) << data_limit/billion << " Gb";
  dls << ends;
  datasize_limit_string = dls.str( );    }
  
  // Generate starting message for assembler "name"

  String beginline = "Begin " + name + " at "
    + start_date + " by " + username + " in\n" + LineOfOutput( "pwd" ) 
    + " on " + getHostName() + ",\nwith datasize limit = "
    + datasize_limit_string + ".\n";
  String foldcommand = origcommand;
  SlashFold( "Invoked by " + origcommand, foldcommand );

  out << "\n"
      << bar << "\n"
      << beginline << "\n"
      << foldcommand
      << bar << endl;
}

/**
 * ShortCommands
 */
void ShortCommands( const vec<String>& commands, vec<String> &shorts )
{
  shorts.resize( commands.size( ), "" );

  vec<String> modules;
  for (int ii=0; ii<(int)commands.size( ); ii++) {
    modules.push_back( commands[ii] );
    if ( modules.back( ).Contains( " " ) )
      modules.back( ) = modules.back( ).Before( " " );
  }
  sort( modules.begin( ), modules.end( ) );
  modules.erase( unique( modules.begin( ), modules.end( ) ), modules.end( ) );
  vec<int> counts( modules.size( ), 0 );

  vec<String>::iterator it;
  for (int ii=0; ii<(int)commands.size( ); ii++) {
    String thecomm = commands[ii];
    while ( thecomm.Contains( " " ) ) thecomm = thecomm.Before( " " );
    it = find( modules.begin( ), modules.end( ), thecomm );
    ForceAssert( it != modules.end( ) );
    int module_id = distance( modules.begin( ), it );
    counts[module_id] += 1;
    shorts[ii] = modules[module_id] + "." + ToString( counts[module_id] );
  }  
}

/**
 * StartCommand
 */
int StartCommand( const vec<String>& commands, const String &start )
{
  if ( start == "" ) return 0;
  return LocateCommand( commands, start );
}

/**
 * StopCommand
 */
int StopCommand( const vec<String>& commands, const String &stop )
{
  if ( stop == "" ) return commands.isize( ) - 1;
  return LocateCommand( commands, stop );
}

/**
 * CheckFilesExist
 */
bool CheckFilesExist( const vec<String> &files, ostream *log )
{
  ofstream devnull ( "/dev/null" );
  ostream &out = log ? *log : devnull;

  vec<String> missing;
  for (int ii=0; ii<(int)files.size( ); ii++)
    if ( ! IsRegularFile( files[ii] ) ) missing.push_back( files[ii] );
  
  if ( missing.size( ) > 0 ) {
    out << "Fatal error, some essential file(s) are missing:\n\n";
    for (int ii=0; ii<missing.isize( ); ii++)
      out << " " << missing[ii] << "\n";
    out << endl;
  }

  return ( missing.size( ) == 0 );
}

