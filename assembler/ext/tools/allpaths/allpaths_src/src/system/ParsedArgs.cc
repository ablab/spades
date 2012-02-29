///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <unistd.h>
#include <setjmp.h>
#include <signal.h>
#include <fstream>
#include <strstream>
#include <iostream>
#include <iomanip>
#include <sys/time.h>
#include <sys/times.h>
#include <sys/resource.h>
#include <sys/file.h> //for flock

#include "LinkTime.h"
#include "MemberOf.h"
#include "system/ParsedArgs.h"
#include "system/System.h"
#include "system/HostName.h"
#include "system/UseGDB.h"

#ifndef FatalErr
     #define FatalErr(message) { cout << message << endl << endl; exit(-1); }
#endif

#define CleanUpAndExit(EXITCODE) { RemoveTheTee(); exit(EXITCODE); }

const Bool         parsed_args::INVALID_BOOL = 2;
const int          parsed_args::INVALID_INT = numeric_limits<int>::min() + 1;
const unsigned int parsed_args::INVALID_UINT = (unsigned int) -1 ;
const longlong     parsed_args::INVALID_LONGLONG = numeric_limits<longlong>::min() + 1;
const double       parsed_args::INVALID_DOUBLE = -1.1234;
const String       parsed_args::INVALID_STRING = "<required>";
const parsed_args *parsed_args::ARGS = NULL;

Bool parsed_args::pretty_help;
String parsed_args::usage_format;  // usage line format
String parsed_args::param_format;  // parameter name format
String parsed_args::header_format; // header line format
String parsed_args::info_format;   // paramter information format
String parsed_args::doc_format;    // paramter documentation format

void parsed_args::PrintTheCommandPretty( ostream& out, const String& prefix )
{
  String bar( "-----------------------------------------------------------"
	      "---------------------" );
  
  std::string hostname = getHostName();
  hostname = hostname.substr(0,hostname.find('.'));
  
  // Get the process ID and pad it to 5 characters.
  String pid_str = ToString( getpid( ) );
  while ( pid_str.size( ) < 5 ) pid_str += " ";
  
  // Write the first line of the command header.
  out << prefix << "\n" << prefix << bar << "\n" << prefix
      << Date( ) << " run on " << hostname << ", pid=" << pid_str
      << " [" << LINK_TIMESTAMP << " R" << SVN_REVISION << " ]\n";
  
  if ( orig_command_ != "" ) {
    String oc;
    SlashFold( orig_command_, oc );
    out << prefix << oc;
    out << "........................................................."
	<< ".......................\n";
  }
  
  String tc;
  SlashFold( TheCommand( ), tc );
  if ( prefix != "" )  {
    ForceAssert( !prefix.Contains("\n") );
    
    for ( int i = 0; i+1 < tc.isize(); i++ ) // tc changes size in the loop, which is ok
      if ( tc[i] == '\n' )
	tc.replace( i+1, 0, prefix );
  }
  out << prefix << tc << prefix << bar << "\n" << prefix << endl;
}

parsed_args::parsed_args( int argc, char *argv[] )
  : command_( Basename(argv[0]) ), doc_(0), outputRedirectionPipe( NULL )
{
     ARGS = this;
     char* pretty_help_ptr = getenv( "ARACHNE_PRETTY_HELP" );
     if (pretty_help_ptr != 0 && String(pretty_help_ptr) == "Color") {
       usage_format = START_BOLD;
       param_format = START_MAGENTA;
       header_format = START_BOLD;
       info_format = START_LIGHT_BLUE;
       doc_format = "";
       pretty_help = True;
     } else if (pretty_help_ptr != 0 && String(pretty_help_ptr) == "Bold") {
       usage_format = START_BOLD;
       param_format = START_BOLD;
       header_format = START_BOLD;
       info_format = "";
       doc_format = "";
       pretty_help = True;
     } else {
       pretty_help = False;
     }
     name_.reserve( argc - 1 );
     def_.reserve( argc - 1 );
     used_.reserve( argc - 1 );
     if ( argc == 1 || ( argc > 1 && (strcmp( argv[1], "-h" ) == 0  || strcmp( argv[1], "--help" ) == 0 ) ) ) {
       get_help_ = True;
       if (argc == 3)
	 get_help_command_ = String(argv[2]);
     } else if ( argc == 2 && (strcmp( argv[1], "-v" ) == 0  || strcmp( argv[1], "--version" ) == 0 ) ) {
       PrintVersion();
       exit(0);
     } else
       {    static Bool included_args;  // static because of longjmp below
       included_args = False;
          get_help_ = False;
          for ( int i = 1; i < argc; i++ )
          {    if ( String(argv[i]).Contains( "@@", 0 ) )
               {    included_args = True;
                    String source = String(argv[i]).After( "@@" );
                    if ( !IsRegularFile(source) )
                    {     PRINT2(i, argv[i]);
                          InputErr( "Failed to locate " << source << "." );     }
                    Ifstream( in, source );
                    String line, x, inargs;
                    while(in)
                    {    getline( in, line );
                         if ( !in ) break;
                         if ( line.Contains( "#", 0 ) ) continue;
                         istrstream ini( line.c_str( ) );
                         while(ini)
                         {    ini >> x;
                              if ( x == "" ) break;
                              if ( !ParseString(x) )
                              {    PRINT2(i, argv[i]);
	                           InputErr( "Could not parse argument " << x 
                                        << " from file " << source 
                                        << "." );    }    }    }    }
               else if ( ! this->ParseString( argc, argv, i) )
               {    PRINT2(i, argv[i]);
	            InputErr( "You did not invoke this program with " <<
		          "white-space-free arguments of the form PARAMETER=VALUE." 
		          << "\n\nTo see the syntax for this command, type \"" 
		          << command_ << " -h\"." );    }    }

	  MemMonitor();   // Launch MemMonitor if required
	  ParseTheTee();	  // Create TEE if required

          // Print the command name, etc.
     
          if ( TheCommand( ).Contains( "Assemble ", 0 ) ) return;
          if ( TheCommand( ).Contains( "Rerun ", 0 ) ) return;
          if ( TheCommand( ).Contains( "PackSource ", 0 ) ) return;
          if ( TheCommand( ).Contains( " NO_HEADER=True" ) ) return;
          if ( TheCommand( ).Contains( " NH=True" ) ) return;

	  // Add any command arguments passed via a file (using @@filename) to
	  // the original command line string for pretty printing. Ignore TEE.

          if (included_args)
          {
	    Bool first = True;
	      for ( int i = 0; i < argc; i++ )
               {    if ( String(orig_command_).substr(0, 4) == "TEE=" )
		       continue;
	            if ( !first ) orig_command_ += " ";
		    first = False;
                    orig_command_ += argv[i];    }    }

          PrintTheCommandPretty( cout );
	  PrintTheTee();    }    }

parsed_args::~parsed_args() {
  ARGS = NULL;
  RemoveTheTee();
}

bool parsed_args::ParseString( int argc, char **argv, int & i) {
  int len = strlen(argv[i]);
  if (i == argc-1 
      || argv[i][len-1] != '=' //e.g. PARAM=blah
      || String(argv[i+1]).Contains("=") ) { 
    return ParseString(argv[i]);
  }
  else { //e.g. PARAM= blah
    String temp(argv[i]);
    temp += argv[i+1];
    ++i; //consume the next argument!
    return ParseString(temp);
  }
}

void parsed_args::FetchArgsFromFile( const String& filename )
{    if ( IsRegularFile(filename) )
     {    Ifstream( in, filename );
          String line, a;
          while(1)
          {    getline( in, line );
               if ( !in ) break;
               if ( line.Contains( "#", 0 ) ) continue;
               if ( line.size( ) == 0 ) continue;
               istrstream iline( line.c_str( ) );
               while(iline)
               {    iline >> a;
                    if ( a.empty( ) ) break;
                    if ( !this->ParseString(a) )
	                 InputErr( "\nThe file " << filename << "\ndoes not contain "
		              << "white-space-free arguments of the form "
		              << "PARAMETER=VALUE." );    }    }    }    }

bool parsed_args::ParseString( const String& string )
{
  int eq = string.Position( "=" );
  if ( eq <= 0 ) 
    return false;

  String name = string.substr( 0, eq );
  if ( !Member( name_, name ) )
  {
    def_.push_back( string.substr( eq+1, string.size( ) ) );
    name_.push_back( name );
    used_.push_back(False);
  }

  return true;
}
  
void parsed_args::ConversionFailure(int i) {
  cout << "I couldn't convert the command line value \"" << def_[i] <<
    "\" for argument \"" << name_[i] << "\" to the required type.\n";
  System( command_ + " -h " + name_[i] );
  cout << "Aborting.\n\n";
  CleanUpAndExit(1);
}

void parsed_args::ValidationFailure(int i, const String& valid ) {
  cout << "Invalid command line value \"" << def_[i]
       << "\" for argument \"" << name_[i] << "\".\n" 
       << "Valid value set/range is " << valid << "\n";
  System( command_ + " -h " + name_[i] );
  cout << "Aborting.\n\n";
  CleanUpAndExit(1);
}

void parsed_args::ValidationFailureDefault(const String& name,
					   const String& deflt,
					   const String& valid ) {
  cout << "Invalid default command line value \"" << deflt
       << "\" for argument \"" << name << "\".\n" 
       << "Valid value set/range is " << valid << "\n";
  System( command_ + " -h " + name );
  cout << "Aborting.\n\n";
  CleanUpAndExit(1);
}

void parsed_args::NoValue(String n)
{    cout << "You must specify a command-line value for " << n << ".\n";
     System( command_ + " -h " + n );
     cout << "Aborting.\n\n";
     CleanUpAndExit(1);    }

String parsed_args::CheckForArgAndAbbreviation( String n, String abbr )
{
  if ( abbr == "" )
    return n;

  bool gotName(False), gotAbbr(False);
  for ( unsigned int i = 0; i < name_.size( ); i++ )
     {
       if ( name_[i] == n )
	 {
	   gotName = True;
	   continue;
	 }
       if ( name_[i] == abbr )
	 {
	   gotAbbr = True;
	 }
     }

  if ( gotName && gotAbbr )
    {
      cout << "You must not specify both an argument-" << n 
	   << " and its abbreviation-" << abbr << ".\n"
	   << "\nTo see the syntax for this command, type \"" 
	   << command_ << " -h\".\n\n";
     CleanUpAndExit(1);    
    }

  if ( gotName )
    return n;
  
  return abbr;

}   

void parsed_args::CheckEnv( const String& n, const String& abbr )
{    for ( int i = 0; i < name_.isize( ); i++ )
          if ( name_[i] == n || name_[i] == abbr ) return;
     int npasses = ( abbr == "" ? 1 : 2 );
     for ( int pass = 1; pass <= npasses; pass++ )
     {    String env_var = "ARACHNE_" + ( pass == 1 ? n : abbr );
          char* env = getenv( env_var.c_str( ) );
          if ( env != 0 )
          {    def_.push_back( String(env) );
               name_.push_back(n);
               used_.push_back(False);    }    }    }

template<>
void parsed_args::GetValue<String>
( String n, String & value, const String & abbr, const String & deflt,
  const String & valid) {
  value = GetStringValue(n, abbr, deflt, valid); 
}
    
template<>
void parsed_args::GetValue<vec<String> >
( String n, vec<String> & value, const String & abbr, const String & deflt,
  const String & valid) {
    ParseStringSet(GetStringValue(n, abbr,deflt), value, false); 
}
    
template<>
void parsed_args::GetValue<vec<int> >
( String n, vec<int> & value, const String & abbr, const String & deflt,
  const String & valid) {
    ParseIntSet(GetStringValue(n, abbr,deflt), value, false); 
}
template<>
void parsed_args::GetValue<vec<longlong> >
( String n, vec<longlong> & value, const String & abbr, const String & deflt,
  const String & valid) {
    ParseLongLongSet(GetStringValue(n, abbr,deflt), value, false); 
}
    
template<>
void parsed_args::GetValue<vec<double> >
( String n, vec<double> & value, const String & abbr, const String & deflt,
  const String & valid) {
    ParseDoubleSet(GetStringValue(n, abbr,deflt), value, false); 
}

template<>
void parsed_args::GetValue<unsigned int>
( String n, unsigned int & value, const String & abbr, const String & deflt,
  const String & valid) {
  value = deflt.IsInt() 
    ? GetUnsignedIntValue(n, abbr, deflt.Int(), valid)
    : GetUnsignedIntValue(n, abbr, INVALID_UINT, valid);
  return;
}

template<>
void parsed_args::GetValue<int>
( String n, int & value, const String & abbr, const String & deflt, 
  const String & valid) {
  value = deflt.IsInt() 
    ? GetIntValue(n, abbr, deflt.Int(), valid)
    : GetIntValue(n, abbr, INVALID_INT, valid);
  return;
}

template<>
void parsed_args::GetValue<longlong>
( String n, longlong & value, const String & abbr, const String & deflt, 
  const String & valid) {
  value = deflt.IsInt() // only checks format, not range
    ? GetLongLongValue(n, abbr, deflt.Int(), valid)
    : GetLongLongValue(n, abbr, INVALID_LONGLONG, valid);
  return;
}


template<>
void parsed_args::GetValue<double>
( String n, double & value, const String & abbr, const String & deflt,
  const String & valid) { 
  value = deflt.IsDouble() 
    ? GetDoubleValue(n, abbr, deflt.Double(), valid)
    : GetDoubleValue(n, abbr, INVALID_DOUBLE, valid);
  return;
}


String parsed_args::GetStringValue( String n, String abbr,
				    String deflt, String valid ) { 

  if (deflt != INVALID_STRING && valid != "" && !IsMemberOf(valid, deflt))
    ValidationFailureDefault(n, deflt, valid);

  CheckEnv( n, abbr );
  String toFind = CheckForArgAndAbbreviation( n, abbr );
  for ( unsigned int i = 0; i < name_.size( ); i++ ) {
    if ( name_[i] == toFind ) {
      used_[i] = True;
      if ( def_[i].size( ) > 0 && def_[i][0] == '"' 
	   && def_[i][ def_.size( ) - 1 ] == '"' ) {
	for ( unsigned int j = 1; j < def_.size( ) - 1; j++ )
	  def_[j-1] = def_[j];
	def_.resize( def_.size( ) - 2 );
      }
      if (valid == "" || IsMemberOf(valid, def_[i]))
	return def_[i];
      else 
	ValidationFailure(i, valid);
    }
  }
  if ( deflt != INVALID_STRING )
    return deflt;

  NoValue(n);
  return GARBAGE;    
}


Bool parsed_args::GetBoolValue( String n, String abbr, Bool deflt )
{ 
  CheckEnv( n, abbr );
  String toFind = CheckForArgAndAbbreviation( n, abbr );
  for ( unsigned int i = 0; i < name_.size( ); i++ )
          if ( name_[i] == toFind ) 
          {    used_[i] = True;
              if ( ToLower(def_[i]) == "true" ) return True;
              else if ( ToLower(def_[i]) == "false" ) return False;
               else FatalErr( "The value of " << n 
                    << " must be either True or False." );    }
     if ( deflt != INVALID_BOOL ) return deflt;
     NoValue(n);
     return GARBAGE;    
}

template<>
void parsed_args::GetValue<Bool>
( String n, Bool & value, const String & abbr, const String & deflt,
  const String & valid) {
    value = (!deflt.IsBool()) ? GetBoolValue(n, abbr) 
    : ("True" == deflt) 
      ?  GetBoolValue(n, abbr, True)
      :  GetBoolValue(n, abbr, False);
    return;
  }


unsigned int parsed_args::GetUnsignedIntValue( String n, String abbr,
					       unsigned int deflt, String valid ) {
  
  if (deflt != INVALID_UINT && valid != "" && !IsMemberOf(valid, deflt))
    ValidationFailureDefault(n, ToString(deflt), valid);
  
  CheckEnv( n, abbr );
  String toFind = CheckForArgAndAbbreviation( n, abbr );
  for ( int i = 0; i < name_.isize( ); i++ ) {
    if ( name_[i] == toFind ) {
      if ( def_[i].size( ) > 0 && def_[i].IsInt( ) ) {
	longlong answer = def_[i].Int( );
	static longlong answer_bound((longlong) 4 *
				     (longlong) (1024 * 1024 * 1024) );
	if ( answer < answer_bound ) {
	  used_[i] = True;
	  unsigned int answerUInt = static_cast<unsigned int>(answer);
	  if (valid == "" || IsMemberOf(valid, answerUInt))
	    return answerUInt;
	  else
	    ValidationFailure(i, valid);
	}
      }
      ConversionFailure(i);
    }
  }
  if ( deflt != INVALID_UINT )
    return deflt;
  
  NoValue(n);
  return 0;
}


int parsed_args::GetIntValue( String n, String abbr, int deflt, String valid ) {

  if (deflt != INVALID_INT && valid != "" && !IsMemberOf(valid, deflt))
    ValidationFailureDefault(n, ToString(deflt), valid);

  CheckEnv( n, abbr );
  String toFind = CheckForArgAndAbbreviation( n, abbr );
  for ( int i = 0; i < name_.isize( ); i++ ) {
    if ( name_[i] == toFind ) {
      if ( def_[i].size( ) > 0 && def_[i].IsInt( ) ) {
	longlong answer = def_[i].Int( );
	static longlong answer_top(INT_MAX);
	static longlong answer_bottom(INT_MIN);
	if ( answer <= answer_top && answer >= answer_bottom ) {
	  used_[i] = True;
	  int answerInt = static_cast<int>(answer);
	  if (valid == "" || IsMemberOf(valid, answerInt))
	    return answerInt;
	  else
	    ValidationFailure(i, valid);
	}
      }
      ConversionFailure(i);
    }
  }
  if ( deflt != INVALID_INT )
    return deflt;
  
  NoValue(n);
  return 0;
}


longlong parsed_args::GetLongLongValue( String n, String abbr, longlong deflt, String valid ) {

  if (deflt != INVALID_LONGLONG && valid != "" && !IsMemberOf(valid, deflt))
    ValidationFailureDefault(n, ToString(deflt), valid);

  CheckEnv( n, abbr );
  String toFind = CheckForArgAndAbbreviation( n, abbr );
  for ( int i = 0; i < name_.isize( ); i++ ) {
    if ( name_[i] == toFind ) {
      if ( def_[i].size( ) > 0 && def_[i].IsInt( ) ) {  // isInt() doesn't range check
	longlong answer = def_[i].Int( );  // Int() returns a longlong, no range checking
	used_[i] = True;
	if (valid == "" || IsMemberOf(valid, answer))
	  return answer;
	else
	  ValidationFailure(i, valid);
      }
      ConversionFailure(i);
    }
  }
  if ( deflt != INVALID_LONGLONG )
    return deflt;
  
  NoValue(n);
  return 0;
}


double parsed_args::GetDoubleValue( String n, String abbr,
				    double deflt, String valid ) {  

  if (deflt != INVALID_INT && valid != "" && !IsMemberOf(valid, deflt))
    ValidationFailureDefault(n, ToString(deflt), valid);

  CheckEnv( n, abbr );
  String toFind = CheckForArgAndAbbreviation( n, abbr );
  for ( unsigned int i = 0; i < name_.size( ); i++ ) {
    if ( name_[i] == toFind ) {
      if ( def_[i].size( ) > 0 && def_[i].IsDouble( ) ) {
	double answer = def_[i].Double( );
	used_[i] = True;
	if (valid == "" || IsMemberOf(valid, answer))
	  return answer;
	else
	  ValidationFailure(i, valid);
      }
      ConversionFailure(i);
    }
  }
  if ( deflt != INVALID_DOUBLE )
    return deflt;

  NoValue(n);
  return GARBAGE;
}


void parsed_args::MemMonitor() {
  int pos = Position(name_, String("MM"));
  if (pos == -1)  pos = Position(name_, String("MEM_MONITOR"));
  if (pos != -1) {
    if ( def_[pos] == "True" ) {
      int pid = getpid();
      String cmd = "MemMonitor";
      if (IsCommandInPath(cmd) ) {
	String mmcmd = cmd + " NH=True PID=" + ToString(pid);
	for (unsigned int mmindex = 0; mmindex < used_.size( ); mmindex++) {
	  if (name_[mmindex].Contains("_MM_")) {
	    mmcmd += " " + name_[mmindex].After("_MM_") + "=" + def_[mmindex];
	  }
	}
	Fork(mmcmd, true);
      }
      used_[pos] = True;
    }
  }    
}

void parsed_args::ParseTheTee()
{
  int pos = Position(name_, String("TEE"));
  if (pos != -1) {
    MakeTheTee(def_[pos]);
    used_[pos] = True;
  }
}

void parsed_args::MakeTheTee(const String filenames)
{
  if (filenames != "") {
    outputRedirectedTo = filenames;
    String pipe_command = "tee " + outputRedirectedTo;
    outputRedirectionPipe = new procbuf( pipe_command.c_str( ), ios::app );
    ostream *pipeout = new ostream( outputRedirectionPipe );
    std::cout.rdbuf( pipeout->rdbuf( ) );
  }
}


void parsed_args::RemoveTheTee()
{
  if ( outputRedirectionPipe && outputRedirectionPipe->is_open() )
    outputRedirectionPipe->close();
}

void parsed_args::PrintTheTee()
{
  if ( outputRedirectedTo != "" ) {
    cout << "\nRedirecting standard output to the following files:\n";
    String files = outputRedirectedTo;
    DeleteLeadingWhiteSpace(files);
    DeleteTrailingWhiteSpace(files);
    files.GlobalReplaceBy(" ", "\n");
    cout << files << "\n\n";
  }
}


void parsed_args::CheckForExtraArgs( bool abort_if_found )
{    Bool junk_found = False;
     for ( unsigned int i = 0; i < used_.size( ); i++ )
          if ( used_[i] == False )
          {
	       if ( name_[i] == "GDB" )
               {    if ( def_[i] == "True" ) use_gdb_for_tracebacks = True;
                    else if ( def_[i] != "False" )
                    {    cout << "illegal value for parameter GDB" << endl;
                         CleanUpAndExit(1);    }
                    continue;    }
               if ( name_[i] == "NO_HEADER" )
               {    if ( def_[i] != "True" && def_[i] != "False" )
                    {    cout << "illegal value for parameter NO_HEADER" << endl;
                         CleanUpAndExit(1);    }
                    continue;    }
               if ( name_[i] == "NH" )
               {    if ( def_[i] != "True" && def_[i] != "False" )
                    {    cout << "illegal value for parameter NH" << endl;
                         CleanUpAndExit(1);    }
                    continue;    }
               if ( name_[i] == "MEM_MONITOR" || name_[i] == "MM" )
               {    if ( def_[i] != "True" && def_[i] != "False" )
		    {    cout << "illegal value for parameter " << def_[i] << endl;
		         CleanUpAndExit(1);    }
		    continue;    }
               if ( name_[i].Contains("_MM_") ) continue;
               if ( name_[i] == "TV" )
		 {    TraceValCommon::SetStopVal( (TraceValCommon::timestamp_t)atoll(def_[i].c_str()));
                    continue;    }
               cout << "You should not have specified a value for " <<
                    name_[i] << ".\n";
               junk_found = True;    }
     if (junk_found) 
       if ( abort_if_found )
       {
	 cout << "\nTo see the syntax for this command, type \"" << command_
	      << " -h\". \n" << endl;
	 CleanUpAndExit(-1);    
       }
       else {
         cout << "Ignoring extra arguments.\n" << endl;
       }
}

/// The special command-line argument TEE="file1 file2..."
/// Redirects cout to the specified files.  This method tells the caller the
/// names of the files (if given) to which we are redirecting output, or the
/// empty string if not redirecting. Allows a module to turn off its own
/// redirection as needed.
String parsed_args::GetOutputRedirection() const {
  return outputRedirectedTo;
}

void parsed_args::SetOutputRedirection(const String filenames) {
  if (outputRedirectedTo != "") {
    cout << "Unable to redirect output to: " << filenames << endl;
    cout << "Output already redirected via TEE parameter." << endl;
  } else {
    MakeTheTee(filenames);
  }
}


void parsed_args::PrintVersion()
{
  cout << command_ << " r" << SVN_REVISION << endl;
  return;
}


void parsed_args::AddArgHelp( const String& arg, 
			      const String& type, 
			      const String& abbreviation,
			      const String& default_value,
			      const String& valid_values,
			      const String& documentation)
{
  arg_help_.push_back(parsed_arg_help
		      (arg, type, abbreviation, default_value, valid_values, documentation));
}

void parsed_args::PrintArgHelp()
{

  bool print_special = (get_help_command_ == "special");    
  if (get_help_command_ != "" && print_special != true) {
    for ( unsigned int i = 0; i < arg_help_.size( ); i++ )
      if ( arg_help_[i].GetName() == get_help_command_) {
	cout << "\n" << arg_help_[i] << "\n\n"
	     << "For more help type: "<< command_ << " -h\n" << endl;
	return;
      }
    cout << "\nUnknown argument: " << get_help_command_ 
	 << "  - No help available.\n\nFor a list of arguments type: "
	 << command_ << " -h\n" << endl;
    return;
  }

  if (parsed_args::pretty_help) cout << usage_format;
  cout << "\nUsage: " << command_ << " arg1=value1 arg2=value2 ..." << endl;
  if (parsed_args::pretty_help) cout << END_ESCAPE;
  if (doc_) {
    WhitespaceFold(String(doc_), cout);
    cout << "\n";
  }
  vec<parsed_arg_help>::iterator firstOpt =
    stable_partition(arg_help_.begin(), arg_help_.end(),
		     mem_fun_ref(&parsed_arg_help::isRequired));
  if (firstOpt != arg_help_.begin()) {
    if (parsed_args::pretty_help) cout << header_format;
    cout << "\nRequired arguments:\n\n";
    if (parsed_args::pretty_help) cout << END_ESCAPE;
    copy( arg_help_.begin(), firstOpt, 
	  ostream_iterator<parsed_arg_help>(cout, "\n") );
  }
  if (firstOpt != arg_help_.end() ) {
    if (parsed_args::pretty_help) cout << header_format;
    cout << "\nOptional arguments:\n\n";
    if (parsed_args::pretty_help) cout << END_ESCAPE;
    copy( firstOpt, arg_help_.end(), 
	  ostream_iterator<parsed_arg_help>(cout, "\n") );
  }
    
  if (print_special) {
    if (parsed_args::pretty_help) cout << header_format;
    cout << "\nSpecial arguments:\n\n";
    if (parsed_args::pretty_help) cout << END_ESCAPE;

    cout << parsed_arg_help( "GDB", "Bool", "", "False", "", "Whether to use GDB for tracebacks." )
	 << "\n"
	 << parsed_arg_help( "NO_HEADER", "Bool", "NH", "False", "",
			     "Whether to suppress the normal command-line header block." )
	 << "\n"
	 << parsed_arg_help( "MEM_MONITOR", "Bool", "MM", "False", "",
			     "Whether or not to fork an instance of MemMonitor. "
			     "All arguments specifiable to MemMonitor can be supplied here by prefixing them with '_MM_'." )
	 << "\n"
	 << parsed_arg_help( "TEE", "String", "", "", "",
			     "Redirect standard out to the supplied space separated list of files." )
	 << "\n";
  } else {
    cout << endl << "To see additional special arguments, type: " << command_ << " --help special" << endl;
  }
  cout << endl;
}

ostream& operator<< (ostream& out, const parsed_arg_help& an_arg)
{
  // Output first line describing the argument briefly.
  if (parsed_args::pretty_help) out << parsed_args::param_format;
  out << an_arg.name_;
  if (parsed_args::pretty_help) out << END_ESCAPE;
  if ( !an_arg.abbreviation_.empty() )
  {    out << ", or ";
       if (parsed_args::pretty_help) out <<  parsed_args::param_format;
       out << an_arg.abbreviation_;
       if (parsed_args::pretty_help) out << END_ESCAPE;    }
  if (parsed_args::pretty_help) out <<  parsed_args::info_format;
  out << " (" << an_arg.type_;
  if ( an_arg.valid_.empty() )
    out << ") ";
  else
    out << ": " << an_arg.valid_ << ") ";
  if ( ! (an_arg.isRequired() || an_arg.default_.empty() ) )
    out << "default: " << an_arg.default_ << " "; 
  if (parsed_args::pretty_help) out << END_ESCAPE;
  
  // Now fold and display the documentation line, if any.
  if (!an_arg.doc_.empty()) {
    if (parsed_args::pretty_help) out <<  parsed_args::doc_format;
    WhitespaceFold(an_arg.doc_, out, 70, "  ");
    if (parsed_args::pretty_help) out << END_ESCAPE;
  }
  return out;
}

String parsed_args::TheCommand( ) const
{    String com = command_;
     for ( unsigned int i = 0; i < name_.size( ); i++ )
     {    Bool need_quotes = False;
          const String& d = def_[i];
	  if ( name_[i] == "TEE" ) continue;
          for ( unsigned int j = 0; j < d.size( )  &&  !need_quotes; j++ )
          {    if ( isspace( d[j] ) ) need_quotes = True;
               if ( d[j] == '{' || d[j] == '}' ) need_quotes = True;
               if ( d[j] == '[' || d[j] == ']' ) need_quotes = True;
               if ( d[j] == '(' || d[j] == ')' ) need_quotes = True;
               if ( d[j] == '%' || d[j] == ',' || d[j] == '=' ) need_quotes = True;
               if ( d[j] == '*' ) need_quotes = True;    }
          if ( !need_quotes ) com += " " + name_[i] + "=" + d;
          else com += " " + name_[i] + "=\"" + d + "\"";    }
     return com;    }

bool parsed_args::RemoveArg(String n, String abbr)
{  
  String toFind = CheckForArgAndAbbreviation( n, abbr );
  for ( int i = 0; i < name_.isize( ); i++ )
    if ( name_[i] == toFind ) {
      name_.erase(name_.begin() + i, name_.begin() + i+1);
      def_.erase(def_.begin() + i,  def_.begin() + i+1);
      used_.erase(used_.begin() + i,  used_.begin() + i+1);
      return true;
    }
  return false;
}
