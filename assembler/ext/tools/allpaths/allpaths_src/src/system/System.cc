///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Note that System.h is included first because it contains some defines which
// affect the other includes.

#include "system/System.h"
#include "system/ErrNo.h"
#include "system/file/FileReader.h"
#include "TokenizeString.h"

#include <cerrno>
#include <locale.h>
#include <time.h>
#include <unistd.h>
#include <signal.h>
#include <ctype.h>

#include <sys/types.h>
#include <sys/wait.h>

#include <algorithm>
#include <fstream>
#include <strstream>
#include <vector>

#include <dirent.h>

int Open( const String& filename, int flags, int mode )
{    int ntries = 0, fd;
     retry: fd = open( filename.c_str( ), flags, mode );
     ++ntries;
     if ( fd < 0 )
     { ErrNo err;
       String msg = "Attempt to open file " + filename + " using flags = "
               + ToString(flags) + " failed" + err.text();

          // Handle case of interrupted function call.

          if ( errno == EINTR && ntries <= 25 )
          {    cout << msg << endl << "Retrying." << endl;
               if ( ntries <= 5 )
               {    cout << "retrying in 10 seconds" << endl;
                    sleep(10);    }
               else if ( ntries <= 15 )
               {    cout << "retrying in 5 minutes" << endl;
                    sleep(300);    }
               else if ( ntries <= 25 )
               {    cout << "retrying in 10 minutes" << endl;
                    sleep(600);    }
               goto retry;    }

          // Handle case where there are too many open files.

          if ( errno == EMFILE )
          {    cout << msg << endl;
               close(50); // Have to close *some* file or rest will fail!
               int pid = getpid( );
               cout << "\nHere is a list of all but one open file:\n" << endl;
               System( "lsof -p " + ToString(pid) );
               cout << "\nThis error is so bad that we have to abort." << endl;
               cout << "Bye." << endl;
               CRD::exit(1);    }

          // Bail.

          FatalErr(msg);    }
     return fd;    }

int SystemInternal( String command, const char *shell )
{
  flush(cout);

  // Theoretically, one could call system(...).  However, we don't, because
  // system(...) can fail if most of a machine's memory is in use.  Instead
  // we directly fork a process.

  int pid, status;

  // It is very important to use vfork instead of fork.  Otherwise, if you're
  // using most of the system memory, the fork attempt may fail.

  pid = vfork();
  if (pid == -1) {
    ErrNo err;
    cout << endl;
    cout << "Fatal error at " << Date( ) << ": can't fork process." << endl;
    cout << "Errno is" << err << endl;
    cout << "Aborting." << endl << endl;
    CRD::exit(1);
  }
  if (pid == 0) { // Child process
    // Start new process running <shell> -c <command>
    // Note path to shell is provided as arg0, by convention.
    // Tell the child process the process id of its actual parent
    // (as opposed to this transient fork process which the child will replace).
    // This will let the child log its parent correctly in <parsed_args::LogTheCommand()>.
    // Only do this if the command is of the form "programname param1 param2..."
    // and not anything more complex, such as "(programname ) &.

    execl( shell, shell, "-c", command.c_str(), NULL);
    // If the execl succeeds we never get here.  If it fails, the only
    // safe way to terminate the child process is with _exit, per
    // vfork() documentation.
    _exit(127);
  }
  while(1) {
    if ( waitpid(pid, &status, 0 ) == -1) {
      if (errno != EINTR) return -1;
    } else {
      if (WIFEXITED(status)) {
        // If program terminated by calling exit() or returning from main(),
        // return the exit status
        return WEXITSTATUS(status);
      } else {
        // program terminated via signal
        return -1;
      }
    }
  }
}

int System( String command )
{
  return SystemInternal(command, "/bin/sh");
}

void SystemSucceed( const String& command )
{    int status = System(command);
     if ( status != 0 )
     {    bool wasSignaled = WIFSIGNALED( status );
          cout << "\n\n";
          if ( wasSignaled ) {
            int signalNum = WTERMSIG( status );
            cout << "Signal " << signalNum << " received while executing command:\n";
          }
          else {
            int exitStatus = WEXITSTATUS( status );
            cout << "Attempt to execute command failed with exit status " << exitStatus << ":\n";
          }
          cout << command << endl;
          FatalErr( "SystemSucceed failed." );    }    }

void SystemSucceedQuiet( const String& command )
{    temp_file out( "/tmp/SystemSucceedQuiet_XXXXXXX" );
     int status = System( command + " > " + out + " 2>&1" );
     if ( status != 0 )
     {    cout << "\n\nAttempt to execute command failed with status " << status << ":\n";
          cout << command << endl << endl;
          Ifstream( in, out );
          String line;
          while(1)
          {    getline( in, line );
               if ( !in ) break;
               cout << line << endl;    }
          FatalErr( "SystemSucceed failed." );    }    }

// Csh is identical to System, except that it uses csh to execute the
// sh command, instead of sh.  At least at some point, doing this
// caused interrupts (CTRL-C) to be handled better.

int Csh( String command )
{
  return SystemInternal(command, "/bin/csh");
}

vector<String> TokenizeCommand(const String &sentence)
{
    istringstream iss(sentence);
    istream_iterator<String> itr = istream_iterator<String>(iss);
    istream_iterator<String> itr_end;
    char quote_flag = 0;
    vector<String> output;
    while(itr != itr_end)
    {
        if ((itr->at(0) == '"') || (itr->at(0) == '`'))
        {
            quote_flag = (*itr)[0];
            output.push_back(itr->substr(1, itr->length() - 1));
        } else if ((itr->at(itr->length() - 1) == quote_flag))
        {
            output.back() += " " + itr->substr(0, itr->length() - 1);
            quote_flag = 0;
        } else if (quote_flag)
        {
            output.back() += " " + *itr;
        } else
        {
            output.push_back(*itr);
        }
        ++itr;
    }
    return output;
}

pid_t Fork(const String& command, Bool ignoreErrors)
{    
    pid_t pid = fork();
    if (pid == -1)
    {
        if (ignoreErrors) { return -1; }
        ErrNo err;
        cout << "Can't fork process.  Exiting..." << endl;
        cout << "Error is" << err << endl;
        CRD::exit(-1);
    }
    if (pid == 0)
    {    
        // check for the command
        vector<String> args = TokenizeCommand(command);
        if (System("which " + args[0] + " > /dev/null") != 0)
        {    
            cout << "WARNING: Could not find the command: \"" << args[0] << "\" in your $PATH" << endl;
            cout << "         Please verify \"" << args[0] << "\" is available in your $PATH." << endl;
            cout << "         I will still attempt to run the command, but this will likely fail." << endl;
        }

        // The call to strdup might look like a memory leak, but it's not!
        // 1) We are about to call execvp, which will collapse our heap.
        // 2) If the execvp fails, we bail, which collapses our heap.
        char *argv[4];
        argv[0] = strdup("sh");
        argv[1] = strdup("-c");
        argv[2] = strdup(command.c_str());
        argv[3] = NULL;
        execv("/bin/sh", argv);
        ErrNo err;
        cout << "ERROR: execv() call to \"" << args[0] << "\" failed!" << endl;
        cout << "Error is" << err << endl;
        CRD::exit(127);
    }
    return pid;
}

// Cp copies one file to another.  To allow for intermittent network
// failures, it will retry many times if the copy fails.  Note that there
// are two possibilities here: first, the "cp" command itself may not be found,
// owing to network failure; second, one or both of the files may be
// inaccessible.  However, it is also possible that cp failure has nothing to do
// with network failure, but happens because either file1 is not there or there
// is a non-intermittent problem writing to file2.  Then what we do here doesn't
// make much sense.

void Cp( String file1, String file2, Bool append )
{    String command;
     if ( !append ) command = "cp " + file1 + " " + file2;
     else command = "cat " + file1 + " >> " + file2;
     for ( int n = 0; n < 10; n++ )
     {    if ( system( command.c_str( ) ) == 0 ) return;
          cout << "copy failed, retrying in 10 seconds" << endl;
          sleep(10);    }
     for ( int n = 0; n < 10; n++ )
     {    if ( system( command.c_str( ) ) == 0 ) return;
          cout << "copy failed, retrying in 5 minutes" << endl;
          sleep(300);    }
     for ( int n = 0; n < 20; n++ )
     {    if ( system( command.c_str( ) ) == 0 ) return;
          cout << "copy failed, retrying in 10 minutes" << endl;
          sleep(600);    }
     if ( !append )
     {    FatalErr( "Error in copying " << file1 << " to " << file2 << "." );    }
     else FatalErr( "Error in copying " << file1 << " >> " << file2 << "." );    }

void CpAppend( String file1, String file2 )
{    Cp( file1, file2, True );    }

void CpAppend( String file1, ostream& file2 )
{    Ifstream( in, file1 );
     while(1)
     {    char c;
          in.get(c);
          if ( !in ) break;
          file2 << c;    }    }

// Cp2 allows file2 to be a directory.

void Cp2( String file1, String file2, Bool append, String chars_to_delete )
{
     if ( !IsRegularFile(file1) )
          FatalErr( "Cp2: file1 = \"" << file1 << "\" does not exist." );
     ForceAssert( file1.size( ) > 0 );

     // If file2 is a directory, replace file2 by file2/last_part_of_file1.

     if ( IsDirectory(file2) )
     {    int slash_pos;
          for ( slash_pos = (int) file1.size( ) - 1; slash_pos >= 0; slash_pos-- )
               if ( file1[slash_pos] == '/' ) break;
          file2 += '/';
          for ( int i = slash_pos + 1; i < (int) file1.size( ); i++ )
               file2 += file1[i];    }

     // Test for identity.

     if ( file1 == file2 )
          FatalErr( "Cp2: Attempt to copy " << file1 << " to itself." );

     // Try to open file1 for read access.
     FileReader fr(file1.c_str());

     // If file2 doesn't exist, we attempt to set its mode to that of file1.
     mode_t mode1 = fr.getStat().st_mode;

     int file2_desc;
     if ( !append )
          file2_desc = open( file2.c_str( ), O_WRONLY | O_CREAT | O_TRUNC, mode1 );
     else file2_desc = open( file2.c_str( ), O_WRONLY | O_CREAT | O_APPEND, mode1 );
     if ( file2_desc < 0 )
     { ErrNo err;
       FatalErr( "Cp2 from " << file1 << " to " << file2 << " failed "
               "while opening " << file2 << err );    }

     // Set up the buffer.  Bigger seems to be better.

     int buf_size = 8192;
     char* buf = new char[ buf_size + sizeof(int) ];

     while(1)
     {    int bytes_read = fr.readSome( buf, buf_size );
          if ( bytes_read == 0 ) break;
          if ( chars_to_delete.size( ) > 0 )
          {    int count = 0;
               for ( int i = 0; i < bytes_read; i++ )
               {    int j;
                    for ( j = 0; j < (int) chars_to_delete.size( ); j++ )
                         if ( buf[i] == chars_to_delete[j] ) break;
                    if ( j == (int) chars_to_delete.size( ) )
                    {    if ( i != count ) buf[count] = buf[i];
                         ++count;    }    }
               bytes_read = count;    }
          if ( write( file2_desc, buf, bytes_read ) < 0 )
          { ErrNo err;
            FatalErr( "Cp2 from " << file1 << " to " << file2 << " failed "
                         "while writing " << file2 << err );    } }

     if ( close(file2_desc) < 0 )
     { ErrNo err;
       FatalErr( "Cp2 from " << file1 << " to " << file2 << " failed "
               "while closing " << file2 << err ); }

     delete [ ] buf;    }

void CpIfNeIfExists( String file1, String file2 )
{    if ( file1 == file2 ) return;
     if ( !IsRegularFile(file1) ) return;
     Cp2( file1, file2 );    }

void CpAppend2( String file1, String file2 )
{    Cp2( file1, file2, True );    }

// The following implementation of Cat should be replaced by one that does not
// do a system call.

void Cat( String infile1, String infile2, String outfile )
{    System( "cat " + infile1 + " " + infile2 + " > " + outfile );    }

void Mv( String file1, String file2 )
{    rename( file1.c_str( ), file2.c_str( ) );    }

void Symlink( String existing_file, String name_of_symbolic_link )
{    int answer = symlink( existing_file.c_str( ), name_of_symbolic_link.c_str( ) );
     if ( answer != 0 )
     { ErrNo err;
       FatalErr( "Symlink(" << existing_file << ", "
               << name_of_symbolic_link << ") failed" << err);    }    }

void SymlinkForce( String existing_file, String name_of_symbolic_link )
{    Remove(name_of_symbolic_link);
     Symlink( existing_file, name_of_symbolic_link );    }

String FirstLineOfFile( String filename )
{    Ifstream( in, filename );
     String line;
     getline( in, line );
     return line;    }

String StringOfFile( String filename, int n )
{    Ifstream( in, filename );
     String s;
     for ( int i = 0; i < n; i++ )
          in >> s;
     return s;    }

String StringOfOutput( String command, int n, bool force )
{    temp_file tempfile( "/tmp/tmp_StringOfOutput_XXXXXXX" );
     if ( System( command + " > " + tempfile ) != 0 && !force )
          FatalErr( "The command \"" + command + "\" failed." );
     return StringOfFile(tempfile, n);    }

String LineOfOutput( String command, bool force, bool err_too )
{    temp_file tempfile( "/tmp/tmp_LineOfOutput_XXXXXXX" );
     String C = command + " > " + tempfile;
     if (err_too) C += " 2>&1";
     if ( System(C) != 0 )
     {    if (force) return FirstLineOfFile(tempfile);
          else FatalErr( "The command \"" + command + "\" failed." );    }
     String s = FirstLineOfFile(tempfile);
     return s;    }

vector<String> AllOfOutput(String command) {
    vec<String> lines;
    String line;
    char value;

    FILE *pfile = popen(command.c_str(), "r");
    while ((value = fgetc(pfile)) != EOF) {
        if (value == '\n') {
            lines.push(line);
            line = "";
        } else {
            line += value;
        }
    }
    pclose(pfile);

    return lines;
}

// ==============================================================================
//
// Not every string is an acceptable file name.  But warn the user if you do this!
//
// ==============================================================================

String FilenameSafeString( String bad_fn ) {
  String fn;
  for(unsigned int i=0; i<bad_fn.size(); i++)
    if ( isalnum(bad_fn[i]) || bad_fn[i] == '_'
	 || bad_fn[i] == '-' || bad_fn[i] == '.' )
      fn += bad_fn[i];
    else
      fn += '_';
  return(fn);
}

Bool Stale( const String& path )
{    struct stat buffer;
     int answer = stat( path.c_str( ), &buffer );
     return answer != 0 && errno == ESTALE;    }

// Stat: this is stat, but is robust to certain intermittent filesystem problems.
// More specifically, if it encounters a "no such device" error (ENODEV), it retries
// a number of times.  This is there because of a problem which we experienced
// once.  Perhaps after further experience we can refine the behavior of this
// function.  Added: ETIMEDOUT and ESTALE treated like ENODEV.

int Stat( const char *path, struct stat *buffer )
{    int ntries = 0;
     while(1)
     {    int answer = stat( path, buffer );
          if ( answer != 0
               && ( errno == ENODEV || errno == ETIMEDOUT || errno == ESTALE )
               && ntries <= 25 )
          {    ++ntries;
               String error_message;
               if ( errno == ENODEV )
                    error_message = "errno = ENODEV (no such device); ";
               if ( errno == ETIMEDOUT )
                    error_message = "errno = ETIMEDOUT (connection timed out); ";
               if ( errno == ESTALE )
                    error_message = "errno = ESTALE (stale NFS file handle); ";
               if ( ntries <= 5 )
               {    cout << "stat of " << path << " failed with "
                         << error_message << "retrying in 10 seconds" << endl;
                    sleep(10);    }
               else if ( ntries <= 15 )
               {    cout << "stat of " << path << " failed with "
                         << error_message << "retrying in 5 minutes" << endl;
                    sleep(300);    }
               else if ( ntries <= 25 )
               {    cout << "stat of " << path << " failed with "
                         << error_message << "retrying in 10 minutes" << endl;
                    sleep(600);    }    }
          else return answer;    }    }

bool IsDirectory( String fn )
{    struct stat buf;
     if ( Stat( fn.c_str( ), &buf ) != 0 )
     {    if ( errno != ENOENT && errno != ENOTDIR )
          { ErrNo err;
            FatalErr( "IsDirectory(" << fn << ") failed" << err );    }
          return false;    }
     return S_ISDIR(buf.st_mode);    }

bool IsRegularFile( String fn )
{    struct stat buf;
     if ( Stat( fn.c_str( ), &buf ) != 0 )
     {    if ( errno != ENOENT && errno != ENOTDIR && errno != ELOOP )
          { ErrNo err;
            FatalErr( "IsRegularFile(" << fn << ") failed" << err );    }
          return false;    }
     return S_ISREG(buf.st_mode);    }

bool IsSymbolicLink( String fn )
{    struct stat buf;
     if ( lstat( fn.c_str( ), &buf ) != 0 )
     {    if ( errno != ENOENT && errno != ENOTDIR )
          { ErrNo err;
            FatalErr( "IsSymbolicLink(" << fn << ") failed" << err );    }
          return false;    }
     return S_ISLNK(buf.st_mode);    }

String ReadSymbolicLink( String const& fn )
{
    char buf[PATH_MAX+1];
    int len = readlink(fn.c_str(),buf,PATH_MAX);
    if ( len < 0 )
    { ErrNo err;
      FatalErr( "ReadSymbolicLink(" << fn << ") failed" << err );   }
    buf[len] = 0;
    return buf;
}

Bool IsSomeSortOfFile( const String& fn )
{    struct stat buf;
     if ( Stat( fn.c_str( ), &buf ) != 0 )
     {    if ( errno != ENOENT && errno != ENOTDIR )
          { ErrNo err;
            FatalErr( "IsSomeSortOfFile(" << fn << ") failed" << err );    }
          return False;    }
     return True;    }

bool AreSameFile( String fn1, String fn2 ) {
  struct stat buf1, buf2;
  return( Stat( fn1.c_str(), &buf1 ) == 0 &&
	  Stat( fn2.c_str(), &buf2 ) == 0 &&
	  buf1.st_dev == buf2.st_dev &&
	  buf1.st_ino == buf2.st_ino );
}

String Date(bool iso8601 )
{    char nowstr[80];
     time_t nowbin;
     const struct tm *nowstruct;
     static bool locale_has_been_set = false;
     if( ! locale_has_been_set ) {
       (void)setlocale(LC_ALL, "");
       locale_has_been_set = true;
     }
     if (time(&nowbin) == (time_t) - 1) return "(date unavailable - time failed)";
     nowstruct = localtime(&nowbin);
     bool failed;
     if (iso8601)
          failed = (strftime(nowstr, 80, "%Y-%m-%dT%H:%M:%S", nowstruct) == (size_t) 0);
     else
          failed = (strftime(nowstr, 80, "%a %b %d %H:%M:%S %Y", nowstruct) == (size_t) 0);
     if (failed)
          return "(date unavailable - strftime failed)";
     return String(nowstr);    }

void Mkdir777( String dn )
{    if ( IsRegularFile(dn) )
     {    FatalErr( "I can't make a directory named " << dn <<
               " because it is already a regular file." );    }
     if ( !IsDirectory(dn) )
     {    if ( mkdir( dn.c_str( ), 0777 ) != 0 )
          { ErrNo err;
            System( "/bin/ls -ld " + dn );
            FatalErr( "Mkdir777 failed on " << dn << err );    }    }    }

void Mkpath(String dn) {
        ForceAssert( dn != "" );
	vec<String> tokens;
	Tokenize(dn, '/', tokens);

	String origin = (dn[0] == '/') ? "/" : "./";

	for (unsigned int i = 0; i < tokens.size(); i++) {
		String dir = origin;
		for (unsigned int j = 0; j <= i; j++) {
			dir += tokens[j] + "/";
		}
		Mkdir777(dir);
	}
}

void Remove( String f )
{    remove( f.c_str( ) );    }

void Rename( String from, String to )
{    rename( from.c_str( ), to.c_str( ) );    }

void RequireDirectory( String fn )
{    if ( !IsDirectory(fn) )
          FatalErr( "I was hoping that " << fn << " was a directory, but "
               << "apparently it isn't.  Perhaps it doesn't exist at all." );    }

void RequireRegularFile( String fn )
{    if ( !IsRegularFile(fn) )
          FatalErr( "I was hoping that " << fn << " was a regular file, but "
               << "apparently it isn't.  Perhaps it doesn't exist at all." );    }

void RequireWritePermission( String fn, ios::openmode mode ) {
  if ( fn != "/dev/null" ) {
    if (!IsRegularFile(fn) || ! (mode & ios::app) ) {//create an empty file.
      close( creat( fn.c_str( ), 0664 ) );
    }
    chmod( fn.c_str( ), 0664 );
    if ( access( fn.c_str( ), W_OK ) != 0 )
      FatalErr( "I need write permission on " << fn <<
		" but don't have it." );
  }
}

void SetDatasizeLimitMb( int n )
{    rlimit data;
     longlong top = (longlong)(n) * 1000000;
     data.rlim_cur = top;
     data.rlim_max = top;
     if ( setrlimit( RLIMIT_DATA, &data ) != 0 )
          FatalErr( "Attempt to set datasize to " << n << " failed." );    }

vector<String> AllFiles( String dirname )
{    vector<String> allfiles;
     temp_file tempfile( "/tmp/tmp_AllFiles_XXXXXXX" );
     System( "/bin/ls -1 " + dirname + " > " + String(tempfile) );
     Ifstream( in, tempfile.c_str( ) );
     String line;
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) break;
          allfiles.push_back(line);    }
     return allfiles;    }

vector<String> AllFilesByTime( String dirname )
{    vector<String> allfiles;
     temp_file tempfile( "/tmp/tmp_AllFilesByTime_XXXXXXX" );
     System( "/bin/ls -1t " + dirname + " > " + String(tempfile) );
     ifstream temp( tempfile.c_str( ) );
     while (temp)
     {    String f;
          temp >> f;
          if ( f != "" ) allfiles.push_back(f);
          else break;    }
     return allfiles;    }

// Find all the files with prefix in directory
vector<String> AllFilesWithPrefix( const String &prefix, const String &dirname )
{
    vector<String> allfiles;
    DIR* dir_ptr = opendir( dirname.c_str( ) );

    if ( !dir_ptr )
    {
      switch ( errno )
      {
        case ENOENT:
          cout << "Error: " << dirname << " does not exist." << endl;
          break;
        case ENOTDIR:
          cout << "Error: " << dirname << " is not a directory." << endl;
          break;
        default:
        {
            ErrNo err;
            cout << "Error opening director" << err << endl;
            break;
        }
      }

      TracebackThisProcess();
    }

    // Walk over the entries of this dir.  If the entry has the prefix we're
    // interested in, save its complete path.

    for ( struct dirent* dirent_ptr = readdir(dir_ptr); dirent_ptr != 0;
          dirent_ptr = readdir(dir_ptr) )
    {
      if ( !strncmp( dirent_ptr->d_name, prefix.c_str(), prefix.size() ) )
        allfiles.push_back( dirname + "/" + dirent_ptr->d_name );
    }

    closedir( dir_ptr );
    return allfiles;
}

void OpenIfstream( ifstream& i, String f )
{    RequireRegularFile(f);
     i.open( f.c_str( ) );
     if ( !i ) FatalErr( "Problem opening " << f << " as an ifstream." );    }

void OpenOfstream( ofstream& o, String f, ios_base::openmode mode)
{    RequireWritePermission(f, mode);
     o.open( f.c_str( ), ios::app | ios::out );
     if ( !o ) FatalErr( "Problem opening " << f << " as an ofstream." );
     chmod( f.c_str( ), 0664 );    }

void OpenOfstream( ofstream& o, String s, String f ,
		   ios_base::openmode mode)
{    RequireWritePermission(f, mode);
     o.open( f.c_str( ), ios::app | ios::out );
     if ( !o )
          FatalErr( "My attempt to open " << s << " as an ofstream for "
               << "file " << f << " has failed." );
     chmod( f.c_str( ), 0664 );    }

float SafeQuotient( longlong numerator, longlong denominator )
{
  ForceAssert( denominator != 0 );
  // Use doubles internally so we don't lose too much precision in this computation
  return double(numerator) / double(denominator);
}

// FileSize: return the size of a file, which must exist.

longlong FileSize( String filename )
{
  return FileReader(filename.c_str()).getSize();
}

// LineCount: return the number of lines in the given file, which must
// be readable.

int LineCount( const String& filename )
{
  // Use old-school file descriptors for an edge on speed.  Since
  // we're doing pretty brain-dead I/O, the clumsiness of this
  // interface is less galling than in most situations.
  FileReader fr(filename.c_str());

  int line_count = 0;

  // read in 16k at a time
  const int buffer_size = 16*1024;
  char buffer[buffer_size+1];
  char* buffer_end = buffer;

  int bytes_read;

  while ( ( bytes_read = fr.readSome(buffer, buffer_size) ) )
  {
    // calculate where valid data ends
    buffer_end = buffer + bytes_read;

    // walk through the buffer, counting newlines using memchr()
    char* newline_ptr = buffer;
    while ( ( newline_ptr = (char*) memchr( newline_ptr, '\n',
					    buffer_end - newline_ptr ) ) )
    {
      ++newline_ptr;
      ++line_count;
    }
  }

  return line_count;
}

void Dot( ostream& log, unsigned int pass )
{    unsigned int passp = (pass+9)/10 * 10;
     String s = ToString(passp);
     if ( (pass % 10) < 10 - s.size( ) ) log << ".";
     else log << s[ (pass % 10) - (10 - s.size( )) ];
     if ( (pass % 10) == 9 )
     {    if ( passp/10 % 7 == 0 ) log << "\n";
          else log << " ";    }
     flush(log);    }

//

const longlong max_read_or_write = 2147479552; // 2^31 - 2^12

void ReadBytes( int filedes, const void* buffer0, longlong nbytes )
{    char* buffer = (char*) buffer0;
     longlong orig_nbytes = nbytes;
     ssize_t answer;
     while( nbytes > max_read_or_write )
     {    answer = read( filedes, buffer, max_read_or_write );
          if ( answer != max_read_or_write ) goto fail;
          buffer += max_read_or_write;
          nbytes -= max_read_or_write;    }
     answer = read( filedes, buffer, nbytes );
     if ( answer == nbytes ) return;
fail:
     ErrNo err;
     FatalErr( "Read of " << orig_nbytes << " bytes using file descriptor "
          << filedes << " failed" << err );    }

void WriteBytes( int filedes, const void* buffer0, longlong nbytes )
{    char* buffer = (char*) buffer0;
     longlong orig_nbytes = nbytes;
     ssize_t answer;
     while( nbytes > max_read_or_write )
     {    answer = write( filedes, buffer, max_read_or_write );
          if ( answer < 0 ) goto fail;
          buffer += max_read_or_write;
          nbytes -= max_read_or_write;    }
     answer = write( filedes, buffer, nbytes );
     if ( answer >= 0 ) return;
fail:
     ErrNo err;
     FatalErr( "Write of " << orig_nbytes << " bytes using file descriptor "
          << filedes << " failed" << err );    }

void CopyBytes(int fdsource, int fdtarget, longlong nbytes,
	       const int BLOCKSIZE) {
  char * c = new char[BLOCKSIZE];

  while (nbytes > BLOCKSIZE) {
    ReadBytes(fdsource, c, BLOCKSIZE);
    WriteBytes(fdtarget, c, BLOCKSIZE);
    nbytes -= BLOCKSIZE;
  }
  ReadBytes(fdsource, c, nbytes);
  WriteBytes(fdtarget, c, nbytes);
  delete[] c;
}


ssize_t SafeWrite( int filedes, const void *buffer, size_t nbytes )
{    ssize_t answer = write( filedes, buffer, nbytes );
     if ( answer < 0 )
     { ErrNo err;
       FatalErr( "Write of " << nbytes << " bytes using file descriptor "
               << filedes << " failed" << err );    }
     return answer;    }

void PlainFold( const String &any_word, ostream& out, int start, int stop )
{
  int page_width = 80;
  if ( -1 == stop )
    stop = int(any_word.size()) - 1;
  else
    stop = ( stop < int(any_word.size()) ) ? stop : int (any_word.size()) - 1;

  int col = 1;
  for (int ii=start; ii<stop+1; ii++) {
    if ( col > page_width ) {
      out << "\n";
      col = 1;
    }
    out << any_word[ii];
    col++;
  }
}

void WhitespaceFold( const String &text, ostream& out, int width, String indent )
{
  vec<String> lines;
  Tokenize(text,lines,"\n");
  for (int l=0; l != lines.isize(); ++l) {
    vec<String> words;
    Tokenize(lines[l], words);
    int col=0;
    out << "\n" << indent; // start on new line
    for (int i=0; i<words.isize(); ++i) {
      col += 1+words[i].size();
      if (col>width) {
	out << "\n" << indent; // line separator
	col = words[i].size() + indent.size();
      } else if (i>0) {
	out << " "; // word separator
      }
      out << words[i];
      if (col>width) { // word length > width
	out << "\n" << indent; // another line separator
	col = indent.size();
      }
    }
  }
}

void SlashFold( String s, String& out, const int pagew )
{    out.resize(0);
     String leftmarg = "";
     int n = s.Before( " " ).size( );
     for ( int i = 0; i < (int) s.size( ); i++ )
     {    out += leftmarg;
          int char_count = leftmarg.size( );
          int last, right = i + pagew - leftmarg.size( ) - 2;
          leftmarg = "";
          for ( int j = 0; j <= n; j++ )
               leftmarg += " ";
          if ( right >= (int) s.size( ) ) last = s.size( );
          else
          {    for ( last = right-1; last >= i; last-- )
                    if ( s[last] == ' ' ) break;
               if ( last <= i ) last = right;    }
          for ( int k = i; k < last; k++ )
          {    out += s[k];
               ++char_count;    }
          if ( last < (int) s.size( ) && s[last] == ' ' ) ++last;
          if ( last < (int) s.size( ) )
          {    for ( int j = char_count; j < pagew - 1; j++ )
                    out += " ";
               out += "\\";    }
          out += "\n";
          i = last - 1;    }    }

double WallClockTime( )
{    struct timeval tp;
     struct timezone tzp;
     int ftime_ret = gettimeofday( &tp, &tzp ) ;
     ForceAssert( ftime_ret == 0 );
     return (double) tp.tv_sec + (double) tp.tv_usec / 1000000.0;    }


// Note that we accidentally have two versions of Basename.  It's not
// clear that we need both.

String Basename( const String& path )
{
  int end = path.size()-1;
  int i;
  //Remove trailing '/'
  while (path[end] == '/' && end > 1) --end;

  for (i = end; i >= 0; i--)
    if ( path[i] == '/' )
      break;

  return path.substr( i+1, end - i);
}

String FilenameExtension( const String& path )
{
  int end = path.size()-1;
  int i;

  for (i = end; i>=0; i--)
    if ( path[i] == '.' )
      break;

  return path.substr( i+1, end - i);
}

String Dirname( const String& path )
{
  int i = path.size()-1;
  if ( i < 0 )
    return ".";

  // Skip trailing '/'
  while (path[i] == '/' && i > 1) --i;

  for ( ; i >= 0; i-- )
    if ( path[i] == '/' )
      break;

  if ( i > 0 )
    return path.substr( 0, i );
  else if ( i == 0 )
    return "/";
  else
    return ".";
}

String RealPath( const String& path )
{
  char real[ PATH_MAX ];
  if ( 0 == realpath( path.c_str(), real ) ) {
    ErrNo err;
    FatalErr( "RealPath(" << path <<") failed" << err );
  }
  return String(real);
}

int LastModified( const String& fn )
{    if ( !IsRegularFile(fn) && !IsDirectory(fn) ) return -1;
     struct stat buf;
     if ( Stat( fn.c_str( ), &buf ) != 0 )
     { ErrNo err;
       FatalErr( "LastModified(" << fn << ") failed" << err );    }
     return buf.st_mtime;    }

bool IsOlder( String fn1, String fn2 )
{    struct stat buf1, buf2;
     if ( Stat( fn1.c_str( ), &buf1 ) != 0 || Stat( fn2.c_str( ), &buf2 ) != 0 )
     { ErrNo err;
       FatalErr( "IsOlder(" << fn1 << "," << fn2 << ") failed" << err );    }
     return buf1.st_mtime < buf2.st_mtime;    }

/// Copy fn1 to fn2 if fn2 does not exist or is older than fn1.
/// If fn1 does not exist, prints a warning but does nothing.
void CpIfNewer( String fn1, String fn2, Bool ignoreErrors ) {

  // If fn2 is a directory, replace fn2 by fn2/last_part_of_fn1.

  if ( IsDirectory(fn2) ) {
    int slash_pos;
    for ( slash_pos = (int) fn1.size( ) - 1; slash_pos >= 0; slash_pos-- )
      if ( fn1[slash_pos] == '/' ) break;
    ForceAssert( slash_pos >= 0 );
    fn2 += '/';
    for ( int i = slash_pos + 1; i < (int) fn1.size( ); i++ )
      fn2 += fn1[i];
  }

  if ( !IsRegularFile( fn1 ) ) {
      cout << " Warning: could not find " << fn1
	   << " so not updating " << fn2 << endl;
      ForceAssert( ignoreErrors );
  } else {
    if ( !IsRegularFile( fn2 )  ||
	 IsOlder( fn2, fn1 ) )
      Cp( fn1, fn2 );
  }
}


temp_file::temp_file( const String& template_name )
{    String& fn = *this;
     fn = template_name;
     int fd = mkstemp( const_cast<char*>(fn.c_str()) );
     if ( fd == -1 )
     { ErrNo err;
       if ( errno == EMFILE )
        {    cout << "TempFile failed in creating " << fn << " from "
                    << template_name << err << endl;
             cout << "Aborting immediately." << endl;
             CRD::exit(-1);    }
          FatalErr( "TempFile failed in creating " << fn << " from "
               << template_name << " " << err );    }
     close( fd );
     chmod( fn.c_str( ), 0664 );
}

temp_file::~temp_file( )
{    String& fn = *this;
     Remove(fn);    }

String GetTempDir(const String &template_name)
{ String result(template_name);
  mkdtemp( const_cast<char*>(result.c_str()) );
  return result; }

String GetNestedDirsFromKey(const unsigned int key) 
{ String value = ToString(key);
  String dirs; dirs.reserve(2*value.size()-1);
  dirs.push_back(value[0]);
  for (size_t i = 1; i < value.size(); ++i) 
  { dirs.push_back('/');
    dirs.push_back(value[i]); }
  return dirs; }

String GetDirFromKey(const unsigned int key, const unsigned int capacity) 
{ String dir = ToString( key / capacity );
  return dir; 
}

String search_path_for_command( const String &cmd )
{
  String paths( getenv( "PATH" ) );
  String path;
  while( ! paths.empty() )
  {
    int pos = paths.Position(':');
    if ( pos == -1 )
    {
      path = paths;
      paths = "";
    }
    else
    {
      path = paths.substr( 0, pos );
      paths = paths.substr( pos+1, paths.size() - (pos+1) );
    }

    if ( IsRegularFile( path + "/" + cmd ) )
    {
      return path + "/" + cmd;
      break;
    }
  }

  return "";
}

// Note that command_name_of_process should not call (directly or indirectly)
// any code which could invoke FatalErr.  Otherwise, one can get into a
// circle.  That's one reason the code is so long: it is unwound from other
// functions.

String command_name_of_process( int pid )
{    String pscom = "ps -p " + ToString(pid) + " -o command";
     temp_file tempfile( "/tmp/tmp_process_name_XXXXXXX" );
     if ( System( pscom + " > " + tempfile ) != 0 )
     {    cout << "Bummer.  I can't run ps to determine the command name." << endl;
          cout << "failed: " << pscom << endl;
          struct stat buf;
          if ( Stat( tempfile.c_str( ), &buf ) != 0 )
          { ErrNo err;
            cout << "Stat of " << tempfile << " failed" << err << endl;
            CRD::exit(-1);    }
          if ( !S_ISREG(buf.st_mode) )
          {    cout << tempfile << "does not exist" << endl;
               CRD::exit(-1);     }
          {    ifstream in;
               in.open( tempfile.c_str( ) );
               if ( !in )
               {    cout << tempfile << " exists but I can't open it." << endl;
                    remove( tempfile.c_str( ) );
                    CRD::exit(-1);    }
               cout << "Here's what happened:" << endl;
               String line;
               while(1)
               {    getline( in, line );
                    if ( !in ) break;
                    cout << line << endl;    }
               cout << "Exiting." << endl;    }
          remove( tempfile.c_str( ) );
          CRD::exit(-1);    }
     String cmd = StringOfFile(tempfile, 2);
     if ( cmd.Contains( "[", 0 ) && cmd.Contains( "]", -1 ) )
          cmd = cmd.After( "[" ).Before( "]" );
     // if the path isn't fully specified, use path to track down the filename
     if ( cmd[0] != '-' && !cmd.Contains( "/" ) )
     {    String path = search_path_for_command( cmd );
          if ( ! path.empty() ) cmd = path;    }
     if ( cmd.Contains( "./", 0 ) ) cmd = cmd.After( "./" );
     return cmd;    }

vector<String> AllFilesInSource( String source )
{    vector<String> allfiles;
     {    temp_file tempfile( "/tmp/tmp_AllFilesInSource_XXXXXXX" );
          String tmp(tempfile), f;
          if ( System( "sh -c \"echo " + source + "\" > " + tmp + " 2>&1" ) == 0 )
          {    Ifstream( iall, tmp );
               while(iall)
               {    iall >> f;
                    if ( f.empty( ) ) break;
                    allfiles.push_back(f);    }    }
          if ( allfiles.size( ) == 0 )
          {    cout << "Your SOURCE expands to zero files.  Abort.\n";
               CRD::exit(1);    }    }
     return allfiles;    }

void Rmdir( const String& dir )
{    int status = rmdir( dir.c_str( ) );
     if ( status != 0 ) FatalErr( "Rmdir on " << dir << " failed." );    }

String TimeSince( double start )
{    double elapsed = WallClockTime( ) - start;
     ostrstream elapsed_str;
     if ( elapsed < 60.0 )
          elapsed_str << setprecision(3) << elapsed << " seconds";
     else if ( elapsed < 3600.0 )
          elapsed_str << setprecision(3) << elapsed/60.0 << " minutes";
     else if ( elapsed < 86400.0 )
          elapsed_str << setprecision(3) << elapsed/3600.0 << " hours";
     else elapsed_str << setprecision(3) << elapsed/86400.0 << " days";
     elapsed_str << ends;
     return elapsed_str.str( );    }

int64_t MemUsageBytes()
{
#ifdef __linux
    return StringOfFile("/proc/self/statm",2 ).Int()*pageSize();
#else
    struct rusage my_rusage;
    getrusage(RUSAGE_SELF, &my_rusage);
    return my_rusage.ru_maxrss;
#endif
}

namespace
{
    size_t gMaxMemory;
}

void SetMaxMemory( size_t maxMemory )
{
    gMaxMemory = std::min(maxMemory,physicalMemory());
}

size_t GetMaxMemory()
{
    return gMaxMemory ? gMaxMemory : physicalMemory();
}

size_t MemAvailable()
{
    size_t mem = GetMaxMemory();
    size_t used = MemUsageBytes();
    return mem > used ? mem-used : 0ul;
}

void Close( int fd )
{    int status = close(fd);
     if ( status != 0 )
     { ErrNo err;
       FatalErr( "Attempt to close file descriptor " << fd
               << " failed" << err );    }    }

void RightPrecisionOut( ostream& out, const double x, int d )
{    if ( x > -1.0 && x < 1.0 ) out << setprecision(d) << x;
     else
     {    out << setiosflags(ios::fixed) << setprecision(d)
               << x << resetiosflags(ios::fixed);    }    }

String Getenv( const String& var )
{    char* valp = getenv( var.c_str( ) );
     if ( valp == 0 )
     {    cout << "Getenv: environment variable " << var << " is undefined.\n";
          cout << "Abort.\n";
          TracebackThisProcess( );    }
     return valp;    }

void PrintWithSep(ostream & out, unsigned int val, char sep) {
  if ( val == 0 ) {
    out << val;
    return;
  }

  unsigned int blocks[5];
  int i = 0;

  while ( val > 0 ) {
    blocks[i++] = val%1000;
    val/=1000;
  }

  out << blocks[--i];

  while ( i > 0 ) {
    out << sep;
    i--;

    if ( blocks[i] < 10 ) {
      out << "00" << blocks[i];
    } else if ( blocks[i] < 100 ) {
      out << '0' << blocks[i];
    } else out << blocks[i];
  }

}

