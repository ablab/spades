///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef SYSTEM_H
#define SYSTEM_H

#include <cstddef>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <vector>
#include "String.h"
#include "system/SysConf.h"
#include "system/SysIncludes.h"
#include "system/Types.h"
#include "system/Exit.h"

#ifndef InputErr
     #define InputErr(message)                                               \
     {    cout << "\nFatal error at " << Date( ) << ": " << message << endl; \
          cout << "\nInvalid input detected." << endl  << endl;		     \
          CRD::exit(1);    }
#endif

#ifndef FatalErr
     #define FatalErr(message)                                               \
     {    cout << "\nFatal error (pid=" << int(getpid()) << ") at "          \
               << Date( ) << ":\n" << message << endl << endl;               \
          CRD::exit(1);    }
#endif

// ==============================================================================
//
// temp_file: this class is a front end to mkstemp, which appears to be the
// preferred (safest) way to create temporary files.  It takes input a template for
// filename, which is then replaced by an actual filename.  The destructor
// deletes the file.
//
// Sample usage:
//
// temp_file f( "/tmp/bob_XXXXXXX" );
// (do stuff with the file named f)
//
// ==============================================================================

class temp_file : public String {
     public:
     temp_file( const String& template_name );
     ~temp_file( );
};

String GetTempDir(const String &templatename);

// ==============================================================================
//
// GetNestedDirsFromKey: returns a set of nested directories based on the value
// of the key. For example: 10 -> 1/0,  125 -> 1/2/5,  3214 -> 3/2/1/4
// Use to spread sets of files across multiple directories.
//
// ==============================================================================

String GetNestedDirsFromKey(const unsigned int key);


// ==============================================================================
//
// GetDirFromKey: returns a directory based on the value of the key and specified
// capacity. Directory = key / capacity
// Use to spread sets of files across multiple directories.
//
// ==============================================================================

String GetDirFromKey(const unsigned int key, const unsigned int capacity);


// ==============================================================================
//
// Basename: return the final component of a pathname.  Same as basename
// in libgen.h, but libgen.h is not supplied with all platforms.
//
// ==============================================================================

String Basename( const String& path );

// ==============================================================================
//
// FilenameExtension: return component of a file/path name after the final ".".
// If the argument does not contain a ".", the entire string is returned.
//
// ==============================================================================

String FilenameExtension( const String& path );

// ==============================================================================
//
// Dirname: return all but the final component of a pathname.  Same as dirname
// in libgen.h, but libgen.h is not supplied with all platforms.  Specifically,
// ignoring trailing slashes, Dirname() returns the part of the path before the
// last slash, with the following exceptions:
//   1. Dirname( "/" ) returns "/"
//   2. Dirname( "" ) returns "."
//   3. Dirname() on any string with no non-trailing slashes returns "."
//      (which is a generalized version of #2)
//
// ==============================================================================

String Dirname( const String& path );


// ==============================================================================
//
// RealPath: return a resolved pathname, e.g. all symbolic links
// are expanded, extra '/' characters are removed, '/./' and '/../'
// are resolved, etc.
//
// ==============================================================================

String RealPath( const String& path );

// ==============================================================================
//
// Functions to run arbitrary commands.  Use these sparingly: they are not
// robust.
//
// System -- execute a given command
// SystemSucceed -- same as System, but check return status and abort if failure
// SystemSucceedQuiet -- same, but don't show output of command unless it fails
//
// Fork -- returns pid, does not wait for process to finish.
//
// StringOfOutput -- runs a command and returns the nth string which
//                   appears in its output
//
// LineOfOutput -- runs a command and returns the first line of its output
//
// AllOfOutput -- runs a command and returns the output as a vector<String>
//
// AllOfOutput1 -- runs a command and returns the output as a String
//
// ==============================================================================

int System( String command );
/// Abort if command fails.
void SystemSucceed( const String& command );
/// Like SystemSucceed, but pipe stdout and stderr to temp_file
/// /tmp/SystemSucceedQuiet_XXXXXX
void SystemSucceedQuiet( const String& command );

vector<String> TokenizeCommand(const String &sentence);
pid_t Fork( const String& command, Bool ignoreErrors = False );

int Csh( String command );

String StringOfOutput( String command, int n = 1, bool force = false );

String LineOfOutput( String command, bool force = false, bool err_too = False );

vector<String> AllOfOutput(String command);

inline String AllOfOutput1(String command)
{    String all;
     vector<String> line = AllOfOutput(command);
     for ( size_t i = 0; i < line.size( ); i++ )
          all += line[i] + "\n";
     return all;    }

/// Cp copies one file to another, destroying the original contents of the
/// second file.  CpAppend copies one file to another (or a stream),
/// appending onto the
/// original contents of the second file (or stream).  Both Cp and the first
/// version of CpAppend will retry in the event of failure.  However, they both
/// call system, and thus could fail at that point in the even of intermittent
/// network failure.

void Cp( String file1, String file2, Bool append = False );
void CpAppend( String file1, String file2 );
void CpAppend( String file1, ostream& file2 );

/// The following versions of Cp and CpAppend do not make system calls, thereby
/// circumventing the problem that the operating system might not be able to
/// find the "cp" command at a given instant in time.  However, they also do not
/// attempt to retry in the event of failure.
///
/// These assume that file1 exists.
///
/// File2 may be a directory.
///
/// Delete all occurrences of characters in chars_to_delete.

void Cp2( String file1, String file2, Bool append = False,
     String chars_to_delete = String( ) );
void CpAppend2( String file1, String file2 );

/// CpIfNeIfExists: if file1 = file2 (as strings) or file1 does not exist, do
/// nothing.  Otherwise, call Cp2.

void CpIfNeIfExists( String file1, String file2 );

/// Concatenate files.

void Cat( String infile1, String infile2, String outfile );

/// Symlink creates a soft link.  It fails if name_of_symbolic_link
/// already exists.  SymlinkForce deletes name_of_symbolic_link first.

void Symlink( String existing_file, String name_of_symbolic_link );
void SymlinkForce( String existing_file, String name_of_symbolic_link );

void Mv( String file1, String file2 );

/// Open: open a file with specified flags.  Abort if the open fails.
/// The optional mode is only important if flags include O_CREAT and
/// the file is in fact created by the Open() call.  In that case the
/// default mode is read/write/execute permissions for all of u,g,o.

int Open( const String& filename, int flags, int mode = 0664);

inline int OpenForRead( const String& filename )
{    return Open( filename, O_RDONLY );    }

/// Open file for writing, creating if it does not exist.

inline int OpenForWrite( const String& filename, int mode = 0664)
{    return Open( filename, O_WRONLY | O_CREAT, mode );    }

void Close( int fd );

String FirstLineOfFile( String filename );

/// Get the nth string from a file.

String StringOfFile( String filename, int n );

longlong FileSize( String filename );

int LineCount( const String& filename );

bool IsDirectory( String fn );

/// Function: IsRegularFile
/// Returns True if the given file is either a regular file
/// or a symlink that eventually resolves to a regular file.
bool IsRegularFile( String fn );

bool IsSymbolicLink( String fn );

/// Reads the symbolic link.  It's a fatal error if the link cannot be read.
String ReadSymbolicLink( String const& fn );

Bool IsSomeSortOfFile( const String& fn );

/// Returns true if fn1 was last modified before fn2, i.e. fn1 is older than fn2.
bool IsOlder( String fn1, String fn2 );

/// Copy fn1 to fn2 if fn2 does not exist or is older than fn1
void CpIfNewer( String fn1, String fn2, Bool ignoreErrors = False );

/// Check for stale path.

Bool Stale( const String& path );

/// Are these the same file?  Note that symlinks may fool it.

bool AreSameFile( String fn1, String fn2 );

/// LastModified: last modification time of file, in seconds since beginning
/// of 1970.  Returns -1 if file doesn't exist.

int LastModified( const String& fn );

inline double AgeInDays( const String& fn )
{   return double( time(0) - LastModified(fn) ) / ( 24.0 * 3600.0 );    }

inline double AgeInYears( const String& fn )
{   return AgeInDays(fn)/365.0;    }

// Return current date and local time in the following formats:
// default:  Fri Jan 16 14:19:03 2009
// iso8601:   2009-01-16T14:19:03

String Date( bool iso8601 = false );

/// The following has one defect.  If the directory already exists but does not
/// have permissions 777, no action is taken.

void Mkdir777( String dn );

/// Generates all the directories necessary to generate a given path
void Mkpath(String dn);

/// Remove(f) removes the file f.  If f is a soft link, the link is removed, not
/// the file.

void Remove( String f );

void Rename( String from, String to );

void RequireDirectory( String fn );
void RequireRegularFile( String fn );

/// Make sure we can write to this file, abort otherwise.
/// If the mode does not include ios::app, truncate the file to length 0.

void RequireWritePermission( String fn, ios::openmode mode=ios::out );

void SetDatasizeLimitMb( int n );

String Getenv( const String& var );

// ==============================================================================
//
// Tools to return lists of filenames.
//
// ==============================================================================

// Return a list of all files in a given directory.

vector<String> AllFiles( String dirname );

vector<String> AllFilesByTime( String dirname );

vector<String> AllFilesWithPrefix( const String &prefix, const String &directory );

// Return the list of files that "echo" would produce.

vector<String> AllFilesInSource( String source );

// ==============================================================================
//
/// Not every string is an acceptable file name.  But warn the user if you do this!
//
// ==============================================================================

String FilenameSafeString( String wannabe_fn );

// ==============================================================================
//
// Terminal escape sequences.
//
// ==============================================================================

#define START_BOLD "[01m"
#define START_UNDERLINE "[04m"
#define START_BOLD_UNDERLINE "[01;04m"
#define START_YELLOW "[01;33m"
#define START_LIGHT_YELLOW "[00;33m"
#define START_BLUE "[01;34m"
#define START_LIGHT_BLUE "[00;34m"
#define START_GREEN "[01;32m"
#define START_LIGHT_GREEN "[00;32m"
#define START_MAGENTA "[01;35m"
#define START_LIGHT_MAGENTA "[00;35m"
#define START_RED "[01;31m"
#define START_LIGHT_RED "[00;31m"
#define START_CYAN "[01;36m"
#define START_LIGHT_CYAN "[00;36m"
#define START_WHITE "[01;37m"
#define START_LIGHT_WHITE "[00;37m"
#define START_DARK_GREY "[01;30m"
#define END_ESCAPE "[0m"

#define BinRead(FILE, DATA) FILE.read( (char*) &DATA, sizeof(DATA) )
#define BinWrite(FILE, DATA) FILE.write( (char*) &DATA, sizeof(DATA) )

// Ifstream, Ofstream, and OfstreamMode define and open input and output
// streams for a given
// file, performing checks to make sure that the operations will work.
//
// OpenIfstream and OpenOfstream are intended only as machinery to make
// Ifstream and Ofstream work.

void OpenIfstream( ifstream& i, String f );
void OpenOfstream( ofstream& o, String f, ios_base::openmode mode = ios::out);
void OpenOfstream( ofstream& o, String s, String f,
                   ios_base::openmode mode = ios::out);

#define Ifstream(STREAMNAME, FILENAME)                         \
     ifstream STREAMNAME;                                      \
     OpenIfstream( STREAMNAME, FILENAME );

/// Truncate FILENAME to length 0 and open it to STREAMNAME

#define Ofstream(STREAMNAME, FILENAME)                         \
     ofstream STREAMNAME;                                      \
     OpenOfstream( STREAMNAME, #STREAMNAME, FILENAME );

/// Open FILENAME  to STREAMNAME with MODE. Truncate if !(MODE & ios::app).

#define OfstreamMode(STREAMNAME, FILENAME, MODE)               \
     ofstream STREAMNAME;                                      \
     OpenOfstream( STREAMNAME, #STREAMNAME, FILENAME, MODE );

/// PipeOstream: create an ostream STREAMNAME for output which is to be gzipped
/// and deposited in FILENAME.gz.  Note that if you use _exit, then you should
/// do STREAMNAME ## _pipe.close( ) first.

#define PipeOstream(STREAMNAME, FILENAME)                                           \
     if ( IsRegularFile(FILENAME) ) Remove(FILENAME);                               \
     String STREAMNAME ## _name_of_pipe                                             \
          = String("gzip -1 > ") + FILENAME + ".gz";                                \
     procbuf STREAMNAME ## _pipe( STREAMNAME ## _name_of_pipe.c_str( ), ios::out ); \
     ostream STREAMNAME( &STREAMNAME ## _pipe );

/// PipeIstream: create an istream STREAMNAME for input from the
/// gzipped file FILENAME.gz.  Note that if you use _exit, then you
/// should do STREAMNAME ## _pipe.close( ) first.

#define PipeIstream(STREAMNAME, FILENAME)                                           \
     String STREAMNAME ## _name_of_pipe                                             \
          = String("gzip -dc ") + FILENAME + ".gz";                                \
     procbuf STREAMNAME ## _pipe( STREAMNAME ## _name_of_pipe.c_str( ), ios::in ); \
     istream STREAMNAME( &STREAMNAME ## _pipe );

// The following PRINT macros are intended primarily for debugging.

#define PRCORE(X) #X " = " << X
#define PRINT(X) cout << PRCORE(X) << endl;
#define PRINT2(X, Y) cout << PRCORE(X) << ", " << PRCORE(Y) << endl;
#define PRINT3(X, Y, Z) cout << PRCORE(X) << ", " << PRCORE(Y) << ", " \
     << PRCORE(Z) << endl;
#define PRINT4(X, Y, Z, W) cout << PRCORE(X) << ", " << PRCORE(Y) << ", " \
     << PRCORE(Z) << ", " << PRCORE(W) << endl;
#define PRINT5(X, Y, Z, W, T) cout << PRCORE(X) << ", " << PRCORE(Y) << ", " \
     << PRCORE(Z) << ", " << PRCORE(W) << ", " << PRCORE(T) << endl;
#define PRINT6(X, Y, Z, W, T, U) cout << PRCORE(X) << ", " << PRCORE(Y) << ", " \
     << PRCORE(Z) << ", " << PRCORE(W) << ", " << PRCORE(T) << ", "             \
     << PRCORE(U) << endl;
#define PRINT7(X, Y, Z, W, T, U, V) cout << PRCORE(X) << ", " << PRCORE(Y) << ", " \
     << PRCORE(Z) << ", " << PRCORE(W) << ", " << PRCORE(T) << ", "             \
     << PRCORE(U) << ", " << PRCORE(V) << endl;
#define PRINT8(X, Y, Z, W, T, U, V, A) cout << PRCORE(X) << ", " << PRCORE(Y)  \
     << ", " << PRCORE(Z) << ", " << PRCORE(W) << ", " << PRCORE(T) << ", "  \
     << PRCORE(U) << ", " << PRCORE(V) << ", " << PRCORE(A) << endl;
#define PRINT9(X, Y, Z, W, T, U, V, A, B) cout << PRCORE(X) << ", " << PRCORE(Y)  \
     << ", " << PRCORE(Z) << ", " << PRCORE(W) << ", " << PRCORE(T) << ", "     \
     << PRCORE(U) << ", " << PRCORE(V) << ", " << PRCORE(A) << ", "               \
     << PRCORE(B) << endl;

#define PRINT_TO(O, X) O << PRCORE(X) << endl;
#define PRINT2_TO(O, X, Y) O << PRCORE(X) << ", " << PRCORE(Y) << endl;
#define PRINT3_TO(O, X, Y, Z) O << PRCORE(X) << ", " << PRCORE(Y) << ", " \
     << PRCORE(Z) << endl;
#define PRINT4_TO(O, X, Y, Z, W) O << PRCORE(X) << ", " << PRCORE(Y) << ", " \
     << PRCORE(Z) << ", " << PRCORE(W) << endl;
#define PRINT5_TO(O, X, Y, Z, W, T) O << PRCORE(X) << ", " << PRCORE(Y) << ", " \
     << PRCORE(Z) << ", " << PRCORE(W) << ", " << PRCORE(T) << endl;
#define PRINT6_TO(O, X, Y, Z, W, T, U) O << PRCORE(X) << ", " << PRCORE(Y) << ", " \
     << PRCORE(Z) << ", " << PRCORE(W) << ", " << PRCORE(T) << ", "                \
     << PRCORE(U) << endl;
#define DPRINT(X) cout << Date( ) << ": " << PRCORE(X) << endl;
#define DPRINT2(X, Y) cout << Date( ) << ": " << PRCORE(X) << ", " << PRCORE(Y) \
     << endl;
#define DPRINT3(X, Y, Z) cout << Date( ) << ": " << PRCORE(X) << ", " << PRCORE(Y) \
     << ", " << PRCORE(Z) << endl;
#define DPRINT4(X, Y, Z, W) cout << Date( ) << ": " << PRCORE(X) << ", " \
     << PRCORE(Y) << ", " << PRCORE(Z) << ", " << PRCORE(W) << endl;
#define DPRINT5(X, Y, Z, W, T) cout << Date( ) << ": " << PRCORE(X) << ", " \
     << PRCORE(Y) << ", " << PRCORE(Z) << ", " << PRCORE(W) << ", " << PRCORE(T) \
     << endl;
#define DPRINT6(X, Y, Z, W, T, U) cout << Date( ) << ": " << PRCORE(X) << ", " \
     << PRCORE(Y) << ", " << PRCORE(Z) << ", " << PRCORE(W) << ", " << PRCORE(T) \
     << ", " << PRCORE(U)  \
     << endl;
#define DPRINT7(X, Y, Z, W, T, U, V) cout << Date( ) << ": " << PRCORE(X) << ", " \
     << PRCORE(Y) << ", " << PRCORE(Z) << ", " << PRCORE(W) << ", " << PRCORE(T) \
     << ", " << PRCORE(U) << ", " << PRCORE(V)				\
     << endl;
#define DPRINT_TO(O, X) O << Date( ) << ": " << PRCORE(X) << endl;
#define DPRINT2_TO(O, X, Y) O << Date( ) << ": " << PRCORE(X) << ", " << PRCORE(Y) \
     << endl;
#define DPRINT3_TO(O, X, Y, Z) O << Date( ) << ": " << PRCORE(X) << ", " \
     << PRCORE(Y) << ", " << PRCORE(Z) << endl;

// PRF: print a float, together with its byte representation (shown as an
// unsigned int).  PRD is similar.

#define PRF(X)                                     \
     cout << #X " = " << X;                        \
     {    unsigned int *Xi;                        \
          float Xf = X;                            \
          Xi = (unsigned int*) &Xf;                \
          cout << " (" << *Xi << ")\n";    }
#define PRD(X)                                                  \
     cout << #X " = " << X;                                     \
     {    unsigned int *Xi;                                     \
          double Xf = X;                                        \
          Xi = (unsigned int*) &Xf;                             \
          cout << " (" << Xi[0] << " " << Xi[1] << ")\n";    }

// ==============================================================================
//
///Dot, illustrated by example.  If you call Dot(log, pass) with successive
///values pass = 0, 1,..., 99, you get
///
///........10 ........20 ........30 ........40 ........50 ........60 ........70
///........80 ........90 .......100
///
/// as output.
//
// ==============================================================================

void Dot( ostream& log, unsigned int pass );

/// Print a dot only every mod counts of pass.
inline void DotMod (ostream & log, unsigned int pass, unsigned int mod) {
  if (0 == pass % mod) Dot(log, pass / mod);
}

/// Macro: DotPerc
/// Print  percentage of work done in a loop.
///
/// Parameters:
///    pass - the current pass of the loop
///    total - the total number of passes of the loop

#define DotPerc(pass,total) do { \
    static int lastPerc = 0;  \
    int curPerc = int(float(pass) / float(total) * 100.0); \
    if (lastPerc > curPerc) \
      lastPerc = 0; \
    if (int(curPerc) - int(lastPerc) > 10) \
      cout << "..." << curPerc; \
    lastPerc = curPerc; \
  } while(0)

/*  dots_pct:
.....10%.....20%.....30%.....40%.....50%.....60%.....70%.....80%.....90%....100%
*/

inline 
void dots_pct(const size_t i, const size_t n) 
{
  unsigned u0 =  i      * 80 / n;
  unsigned u1 = (i + 1) * 80 / n;
  String s = "";
  for (unsigned u = u0 + 1; u <= u1; u++) {
    if      (u % 8 ==  0) s += "%";
    else if (u % 8 ==  7) s += "0";
    else if (u % 8 ==  6) s += ToString((u / 8 + 1) % 10);
    else if (u     == 77) s += "1";
    else                  s += ".";
    if (u == 80)          s += "\n";
  }
  cout << s << flush;
}

/// Number of bytes of memory in use.
int64_t MemUsageBytes();

/// MemUsage returns the memory usage in kB.
inline int64_t MemUsage() { return MemUsageBytes()/1024ul; }

/// Set the maximum amount of memory you'd like to use.
/// This is not enforced in any manner.  It's just an advisory that's used
/// by the MemAvailable() function which follows.  Setting the advisory
/// to 0 turns it off.  (I.e., it says you're willing to use all of the
/// physical memory.
void SetMaxMemory( size_t maxMemory );

/// Retrieve the memory advisory you set with SetMaxMemory, or all of the
/// physicalMemory() if you haven't set an advisory.
size_t GetMaxMemory();

/// Amount of memory that you just might be able to allocate, or maybe not.
/// Or maybe more. It's hard to say.
size_t MemAvailable();

/// Print memory usage in kB to out
inline void PrintMemUsage( ostream &out = cout )
{    out << "Memory used so far: " << MemUsage( ) << "k." << endl;    }

inline void PrintMemUsage( String stage, ostream &out = cout )
{    out << "Memory used (" << stage << "): " << MemUsage( ) << "k." << endl;    }

/// Print a timestamp, and memory usage in MB, to out
inline void
PrintDateAndMemUsage( ostream &out = cout ) {
  out << Date( ) << ": Memory used: " << ((MemUsage( )/1024)+1) << "M." << endl;
}

/// Print a timestamp, and memory usage in MB, to out
inline void
PrintDateAndMemUsage( String stage, ostream &out = cout ) {
  out << Date( ) << ": Memory used (" << stage << "): " << ((MemUsage( )/1024)+1) << "M." << endl;
}


///MULTIWRITE( vector<ostream*>, stuff to be written )
///If a vector element is zero, nothing is written to the corresponding stream.

#define MULTIWRITE( OSTREAMS, STUFF )                                          \
{    for ( unsigned int multi_i = 0; multi_i < OSTREAMS.size( ); multi_i++ )   \
          if ( OSTREAMS[multi_i] != 0 ) *OSTREAMS[multi_i] << STUFF;    }

// ============================================================================
//
/// Safequotient: compute the quotient of two integers, unless the denominator
/// is zero, in which case abort.
//
// ============================================================================

float SafeQuotient( longlong numerator, longlong denominator );

/*****************************************************************
 *
 * Simple Fstream io for any class with operators << and >> defined
 *
 *****************************************************************/

template<class T> Bool FstreamRead(const String & fn, T & data)
{
  if (IsRegularFile(fn)) {
    ifstream in_st;
    OpenIfstream(in_st, fn);
    in_st >> data;
    return True;
  }
  if (IsRegularFile(fn + ".gz")) {
    String pipe_command = "gzip -dc " + fn + ".gz";
    procbuf inp(pipe_command.c_str(), ios::in);
    istream in_st(&inp);
    if (!in_st)
      FatalErr( "Problem opening " << fn << ".gz." );
    in_st >> data;
    inp.close();
    return True;
  }
  return False;
}

template<class T> void FstreamReadOrFatal(const String & fn, T & data)
{
  if (!FstreamRead(fn, data))
    FatalErr( "Neither " << fn << " nor " << fn << ".gz found.\n" );
}

template<class T> void FstreamReadOrContinue(const String & fn, T & data)
{
  if (!FstreamRead(fn, data)) {} // just ignore 'file not found'
}

template<class T> void FstreamWrite(const String & fn, const T & data)
{
  Remove(fn + ".gz");  // make sure that there's no '<fn>.gz' to conflict with 'fn'
  ofstream out_st;
  OpenOfstream(out_st, fn);
  out_st << data;
}

template<class T> void FstreamWriteGZ(const String & fn, const T & data)
{
  Remove(fn);  // make sure that there's no 'fn' to conflict with '<fn>.gz' 
  String pipe_command = "gzip -1 > " + fn + ".gz";
  procbuf outp(pipe_command.c_str(), ios::out);
  ostream out_st(&outp);
  out_st << data;
  outp.close();
}


// ==============================================================================
//
// READ(FILE, TYPE, DATA)
//
// Example: READ( woof, vec<int>, x )
//
// declares vec<int> x, then reads it in from the file woof.gz (if it exists)
// or else from the file woof.  One or the other must exist.
//
// READN: identical to READ, except that if neither file exists, then no action is
// taken.
//
// READX: identical to READ, except that x is not declared.
//
// Note that an istream operator for TYPE must exist!
//
// ==============================================================================


#define READ( FILE, TYPE, DATA )                                       \
     TYPE DATA;                                                        \
     {    String f(FILE);                                              \
          if ( IsRegularFile(f) )                                      \
          {    ifstream temporary_read_stream( f.c_str( ) );           \
               if ( !temporary_read_stream )                           \
                    FatalErr( "Problem reading " << f << "." );        \
               temporary_read_stream >> DATA;    }                     \
          else if ( IsRegularFile( f + ".gz" ) )                       \
          {    String pipe_command = "gzip -dc " + f + ".gz";          \
               procbuf inp( pipe_command.c_str( ), ios::in );          \
               istream temporary_read_stream( &inp );                  \
               if ( !temporary_read_stream )                           \
                    FatalErr( "Problem reading " << f << ".gz." );     \
               temporary_read_stream >> DATA;                          \
               inp.close( );    }                                      \
          else FatalErr( "Neither " << f << " nor "                    \
               << f << ".gz found.\n" );    }

#define READN( FILE, TYPE, DATA )                                      \
     TYPE DATA;                                                        \
     {    String f(FILE);                                              \
          if ( IsRegularFile(f) )                                      \
          {    ifstream temporary_read_stream( f.c_str( ) );           \
               if ( !temporary_read_stream )                           \
                    FatalErr( "Problem reading " << f << "." );        \
               temporary_read_stream >> DATA;    }                     \
          else if ( IsRegularFile( f + ".gz" ) )                       \
          {    String pipe_command = "gzip -dc " + f + ".gz";          \
               procbuf inp( pipe_command.c_str( ), ios::in );          \
               istream temporary_read_stream( &inp );                  \
               if ( !temporary_read_stream )                           \
                    FatalErr( "Problem reading " << f << ".gz." );     \
               temporary_read_stream >> DATA;                          \
               inp.close( );    }    }

#define READX( FILE, DATA )                                            \
     {    String f(FILE);                                              \
          if ( IsRegularFile(f) )                                      \
          {    ifstream temporary_read_stream( f.c_str( ) );           \
               if ( !temporary_read_stream )                           \
                    FatalErr( "Problem reading " << f << "." );        \
               temporary_read_stream >> DATA;    }                     \
          else if ( IsRegularFile( f + ".gz" ) )                       \
          {    String pipe_command = "gzip -dc " + f + ".gz";          \
               procbuf inp( pipe_command.c_str( ), ios::in );          \
               istream temporary_read_stream( &inp );                  \
               if ( !temporary_read_stream )                           \
                    FatalErr( "Problem reading " << f << ".gz." );     \
               temporary_read_stream >> DATA;                          \
               inp.close( );    }                                      \
          else FatalErr( "Neither " << f << " nor "                    \
               << f << ".gz found.\n" );    }

#define WRITE( FILE, DATA )                            \
{    Remove( FILE + ".gz" );                           \
     Ofstream( temporary_write_stream, FILE );         \
     temporary_write_stream << DATA;    }

#define WRITEG( FILE, DATA )                                      \
{    Remove(FILE);                                                \
     String pipe_command = "gzip -1 > " + String(FILE) + ".gz";   \
     procbuf outp( pipe_command.c_str( ), ios::out );             \
     ostream temporary_write_stream( &outp );                     \
     temporary_write_stream << DATA;                              \
     outp.close( );    }

/// PERCENT_RATIO: As part of a <<-delimited output sequence, output
/// NUMERATOR/DENOMINATOR as a percentage.  Stream state is changed by
/// setting precision to PRECISION and old precision is lost.
///
/// If DENOMINATOR is 0, output "inf".

#define PERCENT_RATIO( PRECISION, NUMERATOR, DENOMINATOR )                    \
   setprecision((NUMERATOR==DENOMINATOR && DENOMINATOR != 0) ? 3 : PRECISION) \
   << ( (0 == DENOMINATOR)                                                    \
          ? numeric_limits<double>::infinity()                                \
          :  (NUMERATOR == DENOMINATOR)                                       \
               ? 100                                                          \
               : 100 * SafeQuotient(NUMERATOR, DENOMINATOR)) << "%"

/// PERCENT_RATIOB: Same as PERCENT_RATIO, but fix WIDTH and number of digits
/// to right of decimal point (PRECISION).  This allows one to get things to
/// line up nicely in a table.

#define PERCENT_RATIOB( WIDTH, PRECISION, NUMERATOR, DENOMINATOR )            \
   setiosflags(ios::fixed) << setprecision(PRECISION) << setw(WIDTH)          \
   << ( (0 == DENOMINATOR) ? numeric_limits<double>::infinity()               \
          : 100 * SafeQuotient(NUMERATOR, DENOMINATOR))                       \
   << resetiosflags(ios::fixed) << "%"

/// VALUE_AND_RATIO: As part of a <<-delimited output sequence, output
/// NUMERATOR and then, inside parentheses, PERCENT_RATIO (q.v.).
#define VALUE_AND_RATIO(PRECISION, NUMERATOR, DENOMINATOR) \
   (NUMERATOR) << " (" \
   << PERCENT_RATIO((PRECISION), (NUMERATOR), (DENOMINATOR)) << ")"

/// Echo(s, filename): write string s (and a newline) to filename.

inline void Echo( String s, String filename )
{    {    ofstream temporary_write_stream( filename.c_str( ), ios::app );
          temporary_write_stream << s << "\n";    }    }

/// ReadBytes and WriteBytes are like the system calls "read" and "write", but abort
/// in the event of failure.  Also they read and write in chunks of at most
/// 2^31 - 2^12 bytes to avoid failures on some operating systems.
///
/// SafeWrite is like WriteBytes, but does not correctly handle large nbytes values,
/// but does return an answer.

void ReadBytes( int filedes, const void* buffer, longlong nbytes );
void WriteBytes( int filedes, const void* buffer, longlong nbytes );
ssize_t SafeWrite( int filedes, const void *buffer, size_t nbytes );

///Copy bytes from source to target in BLOCKSIZE units.
/// Calls ReadBytes and WriteBytes, aborts on failure.
void CopyBytes(int fdsource, int fdtarget, longlong nbytes,
	       const int BLOCKSIZE = 1024);

/// PlainFold: fold a given string. If start and stop are given, it
/// prints the characters in [start, stop].

void PlainFold( const String &any_word, ostream& out, int start=0, int stop=-1 );

/// WhitespaceFold: Starting on a new line, write string to out,
/// folding it to fit in given width, breaking at whitespace. Both
/// tabs and spaces are replaced by a single space in output. Newlines
/// are preserved. If a
/// whitespace-free chunk exceeds width, print it on its own line.
void WhitespaceFold( const String &text, ostream& out, int width=70, String indent="" );

/// SlashFold: fold a given string, and indenting after the first line, and
/// putting a backslash in for each continuation line.

void SlashFold( String s, String& out, const int pagew = 80 );

/// Return the wall clock time, in seconds, as a double.

double WallClockTime( );

/// TimeSince(start): report wall clock time since "start", returning a string with
/// appropriate units (seconds, minutes, hours, or days).

String TimeSince( double start );

// START_TIMER, STOP_TIMER usage, by example:
//
//      START_TIMER( widget_clock, 100 );
//      widget( );
//      STOP_TIMER( widget_clock );
//
// will compute the wall clock time usage by each 100 calls to widget, and print it
// out, every 100 calls, so long as "USE_TIMERS" is defined.

#ifdef USE_TIMERS
     #define START_TIMER( timer_name, freq )                              \
          static int timer_name ## _call_count(0);                        \
          ++timer_name ## _call_count;                                    \
          static double timer_name(0);                                    \
          if ( timer_name ## _call_count > 0                              \
               && timer_name ## _call_count % freq == 0 )                 \
          {    cout << endl << "===== " << setprecision(3) << timer_name  \
                    << " seconds used by last " << freq << " operations"  \
                    << " on timer " << #timer_name << " =====" << endl;   \
               timer_name = 0;    }                                       \
          timer_name -= WallClockTime( );
     #define STOP_TIMER( timer_name )                                     \
          timer_name += WallClockTime( );
#else
     #define START_TIMER( timer_name, freq )
     #define STOP_TIMER( timer_name )
#endif

// EXIT_MAIN_NORMALLY: this is intended as a universal way to flush buffers
// and exit a main program.  It's reason for existence is that
// under gcc, it can use exit(0), which causes termination without
// destroying automatic variables.  This makes termination faster.
// However, under cxx, and perhaps under other compilers, you can't use
// exit(0), because it will terminate without flushing buffers.

#define EXIT_MAIN_NORMALLY { return(0); }


// CheckForCommand: test if the shell command can be found in the current path

inline bool IsCommandInPath(const String &cmd) {
  return (System("which " + cmd + " > /dev/null") == 0);
}

// ===============================================================================
//
// search_path_for_command: if a file named "cmd" is in one of the directories in
// the PATH environment variable, return the path to that file; otherwise return
// an empty string
//
// ===============================================================================

String search_path_for_command( const String &cmd );

// ===============================================================================
//
// command_name_of_process: determine the command name of a given process
//
// ===============================================================================

String command_name_of_process( int pid );

// ===============================================================================
//
// SafeMemcpy: invoke memcpy on chunks of size at most two billion, thereby
// circumventing a bug in gcc's implementation of memcpy.  gcc will replace
// calls to memcpy() with inline code it generates, which will generally be
// faster.  Unfortunately, in 2.95.2, it happens to contain a bug.  cxx uses
// the same strategy unless the -nointrinsics flag is used, but fortunately
// does not have the same bug.
//
// ===============================================================================

#ifdef __GNUC__
inline void SafeMemcpy( void* to, void* from, size_t nbytes )
{    const longlong two_billion = (longlong) 2000000 * (longlong) 1000;
     while(1)
     {    if ( nbytes <= (size_t) two_billion )
          {    memcpy( to, from, nbytes );
               break;    }
          else
          {    memcpy( to, from, two_billion );
               to = ((char*) to) + two_billion;
               from = ((char*) from) + two_billion;
               nbytes -= two_billion;    }    }    }
#else
inline void SafeMemcpy( void* to, void* from, size_t nbytes )
{
  memcpy( to, from, nbytes );
}
#endif

/// Rmdir: remove directory or abort.

void Rmdir( const String& dir );

/// BinaryOverwrite: remove file fn and write t into it.
/// To use this, there has to be a
/// preexisting function BinaryWrite( int fd, const T& t ).

template<class T> void BinaryOverwrite( const String& fn, const T& t ) {
  Remove(fn);
  BinaryWrite( fn, t );
}

/// BinaryWrite: write to fname, delegating to BinaryWrite(int fd,T& t ).

template<class T> void BinaryWrite( const String& fname, const T& t ) {
  int fd = OpenForWrite(fname);
  BinaryWrite( fd, t );
  Close(fd);
}

/// BinaryRead: read fname into t, delegating to BinaryRead(int fd,T& t ).

template<class T> void BinaryRead( const String& fname, T& t ) {
  int fd = OpenForRead(fname);
  BinaryRead( fd, t );
  Close(fd);
}

#define DEFINE_BINIO_FOR(T) \
 inline void BinaryWrite( int fd, const T& x ) \
 {    WriteBytes( fd, &x, sizeof(T) );    } \
 \
 inline void BinaryRead( int fd, T& x ) \
 {    ReadBytes( fd, &x, sizeof(T) );    } \
  typedef T __ ## T ## __binaryIoEatSemicolon__

DEFINE_BINIO_FOR(int);  // the rest are defined in StdMethods.h

// Output d digits to the right of decimal point, or if the number is less
// than 1 in absolute value, output d digits after leading zeros.

void RightPrecisionOut( ostream& out, const double x, int d );

/// Prints value <val> into output stream <out> using <sep> to separate
/// every third position (e.g. 9123456789 -> 9,123,456,789)

void PrintWithSep(ostream & out, unsigned int val, char sep=',');

// Utilities to test if an executable exists.

inline void TestExecutableByRunningIt( String executable, String options )
{    if ( System( executable + " " + options + " > /dev/null 2>&1" ) != 0 )
          FatalErr( "\nIt appears that " << executable
               << " is not properly installed on this system.\n"
               << "The reason why I think this is that when I ran "
               << "\"" << executable << " " << options << "\", it failed.\n"
               << "(Note that one possibility is that you have " << executable
	       << ", but your path does not include it.)\n" );    }

inline void TestExecutableByWhich( String executable )
{    if ( System( "which " + executable + " > /dev/null 2>&1" ) != 0 )
	  FatalErr( "\nIt appears that " << executable 
               << " is not properly installed on this system.\n"
               << "The reason why I think this is that when I ran "
               << "\"which " << executable << "\", it failed.\n"
               << "(Note that one possibility is that you have " << executable
               << ", but your path does not include it.)\n" );    }

// ARG macro.  This is a tool to simplify system calls to executables that are
// invoked with arguments in name=value form.  For example, instead of writing e.g.
//     + " LEN=" + ToString(len) + " OUT=" + OUTFILE
// you can write
//     + ARG(LEN, len) + ARG(OUT, OUTFILE)
// Values that are ints are automatically converted via ToString, and Bools are 
// also converted appropriately.  Note the possibility that a value could be
// converted inappropriately.  See also ParsedArgs.h.

inline String NameValueArg( char const* name, String const& value )
{    return String(" ") + name + "=" + value;    }
inline String NameValueArg( char const* name, int value )
{    return NameValueArg( name, ToString(value) );    }
inline String NameValueArg( char const* name, unsigned value )
{    return NameValueArg( name, ToString(value) );    }
inline String NameValueArg( char const* name, longlong value )
{    return NameValueArg( name, ToString(value) );    }
inline String NameValueArg( char const* name, double value )
{    return NameValueArg( name, ToString(value) );    }
inline String NameValueArg( char const* name, Bool value )
{    return NameValueArg( name, ToStringBool(value) ); }

#define ARG( NAME, VALUE ) NameValueArg( #NAME, VALUE )
#define ARGC(NAME) NameValueArg( #NAME, NAME )

#endif
