/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef RUN_COMMAND_H
#define RUN_COMMAND_H

#include "Vec.h"
#include "String.h"
#include <ostream>

// Wrapper around System, with FatalErr on failure.
int RunCommand( const String &the_command );

// Like RunCommand, plus append cout/cerr to given log file.
int RunCommandWithLog( const String &the_command, const String &log,
		       const bool append = false );

// Just like RunCommand, only returns false iff failed (instead of asserting).
bool RunCommandBool( const String &the_command );

// Locate command for standard implementation of START/STOP in assembly runs.
int LocateCommand( const vec<String>& commands, String key );

// Print out header for an Assembly run.
void AssemblyHeader( const String &name, int argc, char **argv, ostream &out );

// Generate printable (short) versions of the given commands.
void ShortCommands( const vec<String> &commands, vec<String> &shorts );

// Find START (start can be "").
int StartCommand( const vec<String> &commands, const String &start );

// Find STOP (stop can be "").
int StopCommand( const vec<String> &commands, const String &stop );

// Check if given files exist (print message on failure if log is given).
bool CheckFilesExist( const vec<String> &files, ostream *log = 0 );

#endif
