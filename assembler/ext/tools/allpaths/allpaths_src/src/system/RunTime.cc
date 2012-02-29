///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include <unistd.h>

#include "system/RunTime.h"

#include "String.h"
#include "system/Assert.h"
#include "system/Exit.h"
#include "system/file/FileReader.h"
#include "system/MemTracker.h"
#include "system/System.h"
#include "system/SysIncludes.h"
#include "system/Types.h"
#include "system/UseGDB.h"

#include <stdlib.h>
#include <unistd.h>
#include <string>
#include <strstream>


#ifdef __alpha
     #include <excpt.h>
     #include <setjmp.h>
#endif

#ifdef __i386__
#ifndef __solaris
     #include <execinfo.h>
#endif
#endif

#if __ia64__ || __x86_64__
     #include <unwind.h>
#endif

/// ===========================================================================
///
/// ReturnAddress(i), where 0 <= i <= 100: get the return address.  The
/// implementation given here only works for g++, and is only known to work on
/// Intel boxes.  (It does not work on an alpha.)
///
/// ===========================================================================

inline void* ReturnAddress(int i)
{
     #ifdef __GNUC__


          #define RETURN_ADDR(I) if ( i == I ) return __builtin_return_address(I);

          RETURN_ADDR(0)
          RETURN_ADDR(1)
          RETURN_ADDR(2)
          RETURN_ADDR(3)
          RETURN_ADDR(4)
          RETURN_ADDR(5)
          RETURN_ADDR(6)
          RETURN_ADDR(7)
          RETURN_ADDR(8)
          RETURN_ADDR(9)
          RETURN_ADDR(10)
          RETURN_ADDR(11)
          RETURN_ADDR(12)
          RETURN_ADDR(13)
          RETURN_ADDR(14)
          RETURN_ADDR(15)
          RETURN_ADDR(16)
          RETURN_ADDR(17)
          RETURN_ADDR(18)
          RETURN_ADDR(19)
          RETURN_ADDR(20)
          RETURN_ADDR(21)
          RETURN_ADDR(22)
          RETURN_ADDR(23)
          RETURN_ADDR(24)
          RETURN_ADDR(25)
          RETURN_ADDR(26)
          RETURN_ADDR(27)
          RETURN_ADDR(28)
          RETURN_ADDR(29)
          RETURN_ADDR(30)
          RETURN_ADDR(31)
          RETURN_ADDR(32)
          RETURN_ADDR(33)
          RETURN_ADDR(34)
          RETURN_ADDR(35)
          RETURN_ADDR(36)
          RETURN_ADDR(37)
          RETURN_ADDR(38)
          RETURN_ADDR(39)
          RETURN_ADDR(40)
          RETURN_ADDR(41)
          RETURN_ADDR(42)
          RETURN_ADDR(43)
          RETURN_ADDR(44)
          RETURN_ADDR(45)
          RETURN_ADDR(46)
          RETURN_ADDR(47)
          RETURN_ADDR(48)
          RETURN_ADDR(49)
          RETURN_ADDR(50)
          RETURN_ADDR(51)
          RETURN_ADDR(52)
          RETURN_ADDR(53)
          RETURN_ADDR(54)
          RETURN_ADDR(55)
          RETURN_ADDR(56)
          RETURN_ADDR(57)
          RETURN_ADDR(58)
          RETURN_ADDR(59)
          RETURN_ADDR(60)
          RETURN_ADDR(61)
          RETURN_ADDR(62)
          RETURN_ADDR(63)
          RETURN_ADDR(64)
          RETURN_ADDR(65)
          RETURN_ADDR(66)
          RETURN_ADDR(67)
          RETURN_ADDR(68)
          RETURN_ADDR(69)
          RETURN_ADDR(70)
          RETURN_ADDR(71)
          RETURN_ADDR(72)
          RETURN_ADDR(73)
          RETURN_ADDR(74)
          RETURN_ADDR(75)
          RETURN_ADDR(76)
          RETURN_ADDR(77)
          RETURN_ADDR(78)
          RETURN_ADDR(79)
          RETURN_ADDR(80)
          RETURN_ADDR(81)
          RETURN_ADDR(82)
          RETURN_ADDR(83)
          RETURN_ADDR(84)
          RETURN_ADDR(85)
          RETURN_ADDR(86)
          RETURN_ADDR(87)
          RETURN_ADDR(88)
          RETURN_ADDR(89)
          RETURN_ADDR(90)
          RETURN_ADDR(91)
          RETURN_ADDR(92)
          RETURN_ADDR(93)
          RETURN_ADDR(94)
          RETURN_ADDR(95)
          RETURN_ADDR(96)
          RETURN_ADDR(97)
          RETURN_ADDR(98)
          RETURN_ADDR(99)
          RETURN_ADDR(100)
          cout << "Stack too deep.\n";
          exit(1);

     #else

          cout << "ReturnAddress only works if you compile with g++." << endl;
          exit(1);
          return( cout << endl ) ; // returns (void*)

     #endif

           }


/// ===========================================================================
///
/// unwindOneLevel: I don't know what this does.  Only defined for alpha.
///
/// ===========================================================================

#ifdef __alpha
     void unwindOneLevel(
          struct sigcontext *unwindSignalContext,
          void *             runtimeProcedurePtr)
     {    unwind(unwindSignalContext, (pdsc_crd *)runtimeProcedurePtr);    }
#endif


// ============================================================================
///
/// SimplifyFunctionName: Replace its argument, say GerbilSpit( int, int ), with
/// GerbilSplit( ... ).  Some exceptions are dealt with.
//
// ============================================================================

String SimplifyFunctionName( String f )
{
     // If the function name is short, we may as well print it in its entirety:

     if ( f.size( ) < 60 ) return f;

     // Don't mess with operators, as their full name may contain useful
     // information:

     if ( f.Contains( "operator" ) ) return f;

     // I'm not sure what this case would be:

     if ( !f.Contains( "(" ) ) return f;

     // The main cases:

     if ( f.Contains( ")", -1 ) ) return f.Before( "(" ) + "( ... )";
     if ( f.Contains( ") const", -1 ) ) return f.Before( "(" ) + "( ... )";

     // Fallen through, not sure why:

     return f;    }


// ===========================================================================
///
/// IgnoreFunction: Ignore stack dump entries for certain functions.
///
// ===========================================================================

Bool IgnoreFunction( String function )
{    Bool first_print(True);
     if ( first_print && function == "??" ) return True;
     if ( function.Contains( "TracebackThisProcess", 0 ) ) return True;
     if ( function.Contains( "Assert", 0 ) ) return True;
     if ( function.Contains( "arachne_signal_handler", 0 ) ) return True;
     if ( function.Contains( "dump_stack", 0 ) ) return True;
     first_print = False;
     return False;    }


// ==========================================================================
///
/// PrintFrame: print a given return address, prettified.
///
// ==========================================================================

void PrintFrame( void* returnAddress, ostream& out )
{    static int count(0);
     if ( count == 0 ) out << "[raw stack entries =" << flush;
     if ( count > 0 && count % 3 == 0 ) out << "\n                     ";
     else out << " ";

     // The following foolishness has the affect of getting returnAddress
     // printed in a 14-character wide field.  The direct approach didn't work.

     strstream hex_stream;
     hex_stream << hex << returnAddress << endl;
     String sp_addr;
     hex_stream >> sp_addr;
     out << setw(18) << sp_addr << flush;
     ++count;    }

Bool interrupt_detected(False);

// ===========================================================================
///
/// PrintStack: given as input the stack addresses, print a stack dump, exit.
/// This only works with g++.
//
// ===========================================================================

void PrintStack( const vector<void*>& sp_addresses, String command,
     ostream& out, Bool exit_when_done = True, Bool minimal = False )
{    if ( !minimal) out << "]\n" << endl; // closure for calls to PrintFrame
     temp_file tempfile( "/tmp/temp_PrintStack_XXXXXXX" );
     signal( SIGINT, SIG_DFL );
     char* backtracer = getenv( "BACKTRACER" );
     // Below is deprecated and gives a warning in gcc 4.2
     // if ( ! backtracer ) backtracer = "addr2line";
     if ( ! backtracer ) {
        // This means backtracer is null, so allocate space then set
        backtracer = new char[10];
        strcpy(backtracer, "addr2line");
     }

     String addr2line_command( backtracer );
     addr2line_command += " -i -e " + command + " -f -s -C ";

     for ( unsigned int i = 0; i < sp_addresses.size(); ++i )
     {    strstream hex_stream;
          #if __ia64 || __x86_64
               hex_stream << hex << (void *)(((char *)sp_addresses[i])-1) << endl;
          #else
               hex_stream << hex << sp_addresses[i] << endl;
          #endif
          String sp_addr;
          hex_stream >> sp_addr;
          addr2line_command += sp_addr + " ";    }
     addr2line_command += String("> ") + tempfile;
     if ( System( addr2line_command ) != 0 )
     {    out << "Call to " << backtracer << " failed." << endl;
          out << "Perhaps an interrupt was received, or perhaps "
               << backtracer << " could not be found." << endl << endl;
          _exit(1);    }
     ifstream tempstream( tempfile.c_str( ) );
     int depth = 0;
     String function, where;
     for ( unsigned int i = 0; i < sp_addresses.size(); ++i )
     {    if ( !tempstream ) break;
          getline( tempstream, function );
          getline( tempstream, where );
          if ( !IgnoreFunction(function) )
          {    out << depth++ << ". ";
               if ( function == "??" ) out << "??" << endl ;
               else out << SimplifyFunctionName(function) << ", in " << where
                    << endl ;    }
          if ( function == "main" ) break;    }
     remove( tempfile.c_str( ) );
     if ( !minimal ) out << endl;
     if (exit_when_done)
     {    if (interrupt_detected) _exit(1);
          else CRD::exit(1);    }

     delete [] backtracer;
}

// =========================================================================
///
/// dump_stack: generate a human-readable printout of the stack, exit.
/// This is only known to work on Alpha and Intel platforms, and only works under
/// g++.
//
// =========================================================================

  #if __ia64__ || __x86_64__

struct ia64_unwind_struct
{
  static const int max_frames = 101;
  int num_frames;
  void * stacklist[max_frames];
};

_Unwind_Reason_Code ia64_unwind_callback(struct _Unwind_Context *info, void *arg)
{
  unsigned long ip;
  struct ia64_unwind_struct *unwind_data = (struct ia64_unwind_struct *)arg;

  if ( unwind_data->num_frames == unwind_data->max_frames )
    return _URC_END_OF_STACK;

  ip = _Unwind_GetIP( info );
  unwind_data->stacklist[unwind_data->num_frames++] = (void *)ip;

  return _URC_NO_REASON;
}

  #endif

void dump_stack( String command, ostream& out, Bool exit_when_done = True,
     Bool minimal = False )
{    if ( !minimal ) out << "\nDump of stack:" << endl << endl;
     vector<void*> sp_addresses;

#ifdef __solaris
     cout << "Sorry, stack dump does not work on Solaris " << endl;

     exit(1);
#else //__solaris
     #ifndef __GNUC__
          cout << "Sorry, stack dump only works if code was "
               << "compiled with g++." << endl;
          exit(1);
     #endif

     #ifdef __i386__

          void* return_addresses[101];
          int nback = backtrace( return_addresses, 101 );
          for ( int j = 0; j < nback - 1; j++ )
          {    if ( !minimal ) PrintFrame( ReturnAddress(j), out );
               sp_addresses.push_back( ReturnAddress(j) );   }

     #endif

  #if __ia64__ || __x86_64__

          struct ia64_unwind_struct unwind_data;
          unwind_data.num_frames = 0;

          _Unwind_Backtrace( &ia64_unwind_callback, &unwind_data );

          for ( int depth = 0; depth < unwind_data.num_frames; ++depth )
          {
            if ( !minimal ) PrintFrame( unwind_data.stacklist[ depth ], out );
            sp_addresses.push_back( unwind_data.stacklist[ depth ] );
          }

  #endif

     #ifdef __alpha

          #define RETURNADDRREG (26)
          #define FAULTING_ADDRESS sc_traparg_a0

          // Get current execution context.

          jmp_buf context;
          setjmp(context);

          // Set the initial context for the unwind.

          struct sigcontext* signalContextPtr = (struct sigcontext *)context;
          struct sigcontext unwindSignalContext = *signalContextPtr;

          // Discard the frame for dump_stack() and TracebackThisProcess().

          int numLevelsToDiscard = 2;
          for ( int level = 0; level < numLevelsToDiscard; level++ )
          {    unsigned long int programCounter = unwindSignalContext.sc_pc;
               void* runTimeProcedurePtr = find_rpd(programCounter);
               unwindOneLevel(&unwindSignalContext, runTimeProcedurePtr);    }

          // Pick out the return address and program counter.

          unsigned long returnAddress  = unwindSignalContext.sc_regs[RETURNADDRREG];
          unsigned long programCounter = unwindSignalContext.sc_pc;

          // This is the address that caused the fault when we tried to access

          unsigned long faultingAddress = signalContextPtr->FAULTING_ADDRESS;

          // Special cases for bogus program counter values. If the program
          // counter is zero or the fault occurred when we were trying to
          // fetch an instruction (because the program counter itself was bad)
          // then we cannot unwind the stack.

          if (programCounter == 0)
               out << "\nPC is zero - stack trace not available.\n";
          else if (programCounter == faultingAddress)
               out << "\nbad PC (" << faultingAddress
                    << ") - stack trace not available.\n";

          else
          {    unsigned int sameSpCount = 0;

               // Loop through all the stack frames.

               unsigned long stackpointer = 0;
               void *runTimeProcedurePtr;

               while  ((returnAddress != 0) && (programCounter != 0))
               {
                    // Get the run time procedure descriptor for this frame.

                    runTimeProcedurePtr = find_rpd(programCounter);

                    if ( !minimal ) PrintFrame( (void*) returnAddress, out );
                    sp_addresses.push_back( (void*) returnAddress );

                    // Unwind one level.

                    unwindOneLevel(&unwindSignalContext, runTimeProcedurePtr);
                    returnAddress = unwindSignalContext.sc_regs[RETURNADDRREG];
                    programCounter = unwindSignalContext.sc_pc;

                    if ((unsigned int) unwindSignalContext.sc_sp <= stackpointer)
                    {    if ( ++sameSpCount == 10 ) break;    }
                    else
                    {    sameSpCount  = 0;
                         stackpointer = unwindSignalContext.sc_sp;    }    }    }
     #endif

     PrintStack( sp_addresses, command, out, exit_when_done, minimal );

#endif //__solaris
}

// ===========================================================================
/// \fn TracebackThisProcess
/// TracebackThisProcess: generate a human-readable stack dump for this process.
/// This only works if g++ was used, and is only known to work on Alpha and
/// Intel boxes.
///
/// Sometimes gdb will produce a more informative stack dump.  Somehow it does
/// a better job than addr2line does of converting addresses into
/// function names/filenames/line numbers.  If you want to see the gdb stack
/// dump
/// instead of what we would otherwise give, set the global variable
/// use_gdb_for_tracebacks to True.  This is done by ParsedArgs.cc, when the
/// command-line option GDB=True is given.
///
// ==========================================================================

void TracebackThisProcess( ostream& out, Bool exit_when_done, Bool minimal )
{
     // If called by gdb or emacs, crash.

     String mom = command_name_of_process( getppid( ) );
     if ( mom == "gdb" || mom == "emacs" )
     {    out << "called by gdb or emacs; crashing...\n" << flush;
          abort( );    }

     // Get the id of this process.

     int pid = getpid( );

     // Find the executable that was used, or a link to it.  On Linux systems,
     // we check to make sure that it was not deleted or overwritten, and if it
     // was moved, we give a link that points to the moved file.  Otherwise,
     // we don't do this.

     String exe;

     #ifdef __linux

          exe = "/proc/" + ToString(pid) + "/exe";
          char* buf = new char[500];
          if ( readlink( exe.c_str( ), buf, 500 ) < 0 )
          {    out << "Attempt to link to executable failed.  Weird." << endl;
               CRD::exit(1);    }
          String buf_string(buf);
          if ( buf_string.Contains( "(deleted)" ) )
          {    out << "Traceback is impossible because the original executable "
                    << "no longer exists.  Bummer." << endl << endl;
               CRD::exit(1);    }
     #else

          exe = command_name_of_process(pid);

     #endif

     // On Alpha and Intel boxes, we use a built-in stack dump.  Otherwise,
     // call gdb.  None of this will work if you did not compile with g++.

     #if __alpha || __i386__ || __ia64__ || __x86_64__
          if ( !use_gdb_for_tracebacks)
               dump_stack( exe, out, exit_when_done, minimal );
          if ( !exit_when_done ) return;
     #endif

     out << "\nInvoking gdb to backtrack...\n";
     out << "(If process stops, you may have to foreground (fg) it.)" << endl;
     Echo( "attach " + ToString(pid), "bt_file" );
     Echo( "bt",  "bt_file" );
     Echo( "quit",  "bt_file" );
     if ( System( "gdb -batch -x bt_file -q " + exe ) != 0 )
          out << "Bummer.  My call to gdb failed." << endl;
     Remove( "bt_file" );
     if (exit_when_done) CRD::exit(-1);    }


// SA_NOMASK is a linux-ism, so if it's not on your system, define it
#ifndef SA_NOMASK
#define SA_NOMASK SA_NODEFER
#endif

// Stock installs of Mac OS X have a fixed 64M stack size limit.
// Use that as an upper limit instead of recommended 100M.
// WARNING: 100M is the lowest stack size at which we trust an Arachne run.
#ifndef MACOSX
#define STACK_SIZE 100
#else
#define STACK_SIZE 64
#endif

void our_new_handler( )
{    cout << "\nCall to new failed, memory usage before call = "
          << MemUsage( ) << "k." << endl;
     cout << "\nHere is the output of top:\n" << endl;
#ifdef __alpha
     System( "top -b | head -15" );
#else
     System( "top -b -n 1 | tail -n +7 | sort -nr -k 10 | head -15" );
#endif
     cout << "Aborting." << endl;
     TracebackThisProcess( );    }

// ===============================================================================
//
//  SetLimits: Make sure that stacksize is always at least 100M.
//
//  On some systems, this may not have an effect.
//
// WARNING WARNING WARNING!!
// If stacksize is 8192 kbytes (or less), Arachne (or some other program)
// may crash in a totally inexplicable way, causing untold grief.  And I
// don't know the exact value which is needed.  So don't lower the limit set
// here unless you have a very compelling reason.
//
// ===============================================================================

void SetLimits()
{
    // Stock installs of Mac OS X have a fixed 64M stack size limit.
    // Use that as an upper limit instead of recommended 100M.
#ifndef MACOSX
    static unsigned long const STACK_SIZE_KB = 100000;
#else
    static unsigned long const STACK_SIZE_KB = 64000;
#endif
    static char envNameStr[] = "CRD_STKSIZ_REEXEC";
    static char envValStr[] = "CRD_STKSIZ_REEXEC=1";

    unsigned long stackBytes = STACK_SIZE_KB*1024ul;

    rlimit rlim;
    getrlimit( RLIMIT_STACK, &rlim );
    if ( getenv(envNameStr) ) // if we're waking up after a re-exec (see below)
    {
        if ( rlim.rlim_cur < stackBytes ) // and the stack is still too small
            FatalErr("This program requires " << STACK_SIZE_KB
                        << "KB of stack space.\n"
                        << "We adjusted it, and re-exec'd ourselves, but the "
                            "revised limit didn't stick, somehow.\n"
                        << "We don't know what else to do, so we're quitting.");

        if ( putenv(envNameStr) ) // clean up so a sub-program isn't confused
            FatalErr("Unable to remove re-exec marker in environment.");
    }

    // If the value is too small, or if the value is unlimited, reset it.
    // The reason that unlimited is no good is that it causes pthreads to use a
    // default-sized stack which is too small for our purposes, particularly
    // for omp threads over which we have no run-time control of the stack size.
    if ( rlim.rlim_cur < stackBytes || rlim.rlim_cur == RLIM_INFINITY )
    {
        if ( rlim.rlim_max < stackBytes ) // if it's hopeless to try to increase
            FatalErr("This program requires " << STACK_SIZE_KB
                        << "KB of stack space.\n"
                        << "Your system is configured to allow a stack size no "
                            "larger than " << rlim.rlim_max/1024ul << "KB.\n"
                        << "You'll need to ask a system administrator to "
                            "increase this limit for you.\n");

        rlim.rlim_cur = stackBytes;
        if ( setrlimit(RLIMIT_STACK,&rlim) != 0 )
            FatalErr("This program requires " << STACK_SIZE_KB
                        << "KB of stack space.\n"
                        << "It looked like getting that much was possible, but "
                            "when we asked, we were refused.\n"
                        << "We don't know what else to do, so we're quitting.");

        // We need to re-exec this program so that the pthread library and omp
        // will reinitialize using the new value for stack size.

        if ( putenv(envValStr) ) // note that we're attempting a re-exec
            FatalErr("Unable to mark re-exec in environment.");

        // get a buffer for command args
        size_t argMax = sysconf(_SC_ARG_MAX)+1;
        if ( !argMax )
            FatalErr("Can't get max args length from sysconf.");
        char* buf = new char[argMax];
        memset(buf,0,argMax);

        // read the args that were used to invoke this program
        ssize_t totLen = 0;
        if ( true )
        {
            FileReader fr("/proc/self/cmdline");
            totLen = fr.readSome(buf,argMax-1);
        }

        // chop up the whole mess into separate strings
        std::vector<char*> args;
        char* ppp = buf;
        while ( totLen > 0 )
        {
            args.push_back(ppp);
            size_t len = strlen(ppp) + 1;
            ppp += len;
            totLen -= len;
        }
        args.push_back(0);

        cout << "Performing re-exec to adjust stack size." << endl;
        execv("/proc/self/exe",&args[0]);
        FatalErr("Re-exec to adjust stack size failed.");
    }
}

// ===============================================================================
//
// NoDump: turn off core dumps.
//
// ===============================================================================

void NoDump( )
{    rlimit core;
     core.rlim_cur = 0;
     core.rlim_max = 0;
     setrlimit( RLIMIT_CORE, &core );    }

// ===============================================================================
//
// arachne_signal_handler: this is the code that is called when an interrupt
// is detected.
//
// ===============================================================================

static const char* SIGSEGV_msg = 
"Segmentation violation.  An attempt will be made to backtrace,\n"
"but the backtrace itself may fail.";
static const char* SIGFPE_msg = "Arithmetic exception.";
static const char* SIGINT_msg = "Interrupt received (perhaps a ctrl-c).  Stopping.";
static const char* SIGBUS_msg = "SIGBUS interrupt.";
static const char* SIGILL_msg = "Illegal instruction.  Something has gone badly wrong.  Stopping.";
static const char* SIGABRT_msg = "Abort.  Stopping.";
static const char* SIGTERM_msg = "Killed.  Stopping.";

void print_signal_handler_message( int signal_number, siginfo_t* info )
{  
  cout << "\n" << Date() << ".  ";
  switch ( signal_number ) {
    case SIGSEGV: cout << SIGSEGV_msg << endl; break;
    case SIGFPE:  cout << SIGFPE_msg << endl;  break;
    case SIGINT:  cout << SIGINT_msg << endl;  break;
    case SIGBUS:  cout << SIGBUS_msg << endl;  break;
    case SIGILL:  cout << SIGILL_msg << endl;  break;
    case SIGABRT: cout << SIGABRT_msg << endl; break;
    case SIGTERM: cout << SIGTERM_msg << endl; break;
    case SIGCHLD:
        if (info->si_code == CLD_EXITED)
        {
            cout << "Child process " << info->si_pid << " exited with value "
                << info->si_status << ". Stopping." << endl;
        }
        else if (info->si_code == CLD_KILLED)
        {
            cout << "Child process " << info->si_pid
                << " was killed. Stopping." << endl;
        }
        else if (info->si_code == CLD_DUMPED)
        {
            cout << "Child process " << info->si_pid
                << " terminated abnormally. Stopping." << endl;
        }
        else
        {
            cout << "Harmless SIGCHLD - we shouldn't be here." << endl;
        }
        break;        
    default:      
      cout << "Unrecognized signal (" << signal_number << ") received.  Stopping." << endl;
  }
}

extern Bool interrupt_detected;

void arachne_signal_handler( int signal_number, siginfo_t* info, void* context,
    Bool no_ctrlc, Bool nonzero_child )
{
  interrupt_detected = True;

  // ignore SIGCHLDs unless they indicate a "bad" child process termination
  if (signal_number == SIGCHLD &&
      ((info->si_code == CLD_EXITED &&
        (!nonzero_child || info->si_status == 0)) ||
      (info->si_code != CLD_EXITED && info->si_code != CLD_KILLED 
        && info->si_code != CLD_DUMPED)) )
  {
    return;
  }

  if ( signal_number != SIGUSR1 && signal_number != SIGUSR2 )
       print_signal_handler_message( signal_number, info );

  // We do not attempt to trace back on SIGTERM interrupts, which are the interrupts 
  // resulting when one kills a process (without specifying an interrupt type).
  // The reason for this is that our traceback process includes a stack trace,
  // and to generate the stack trace, we have to fork a process, which in turn
  // asks the operating system for memory.  Under certain circumstances,
  // this can wreak havoc.  Once such circumstance is the simultanous killing of 
  // multiple large memory processes.  We did this once and the operating system
  // did not cope well.
  //
  // However, if you do want to kill a process, and get a traceback, this is still
  // possible: use "kill -INT".

  if ( signal_number == SIGTERM ) _exit(1);
  if ( no_ctrlc && signal_number == SIGINT ) _exit(1);

  if ( signal_number == SIGUSR1 || signal_number == SIGUSR2 )
  {    String tracefile = "/tmp/traceback_from_process_" + ToString( getpid( ) );
       Ofstream( out, tracefile + "_initial" );
       TracebackThisProcess( out, False, signal_number == SIGUSR1 );
       out.close( );
       Rename( tracefile + "_initial", tracefile );    }
  else
  {    cout << "\nGenerating a backtrace..." << endl;
       TracebackThisProcess( cout, True, False );
       //  something not right, ugly fix attempt
       printf( "%s\n%s\n", "Back from backtrace.  This should not have occurred, "
            "but probably you don't care.", "Exiting." );
       _exit(1);    }

  interrupt_detected = False;
}

void arachne_signal_handler_standard( int signal_number, siginfo_t* info,
    void* context )
{    arachne_signal_handler( signal_number, info, context, False, False  );    }

void arachne_signal_handler_no_ctrlc_traceback( int signal_number,
    siginfo_t* info, void* context )
{    arachne_signal_handler( signal_number, info, context, True, False );    }

void arachne_signal_handler_nonzero_child(int signal_number, siginfo_t* info,
    void* context)
{
    arachne_signal_handler(signal_number, info, context, False, True);
}

// ===============================================================================
//
// ArachneInterruptHandler: set up run-time interrupt catching.
//
// ===============================================================================

void ArachneInterruptHandler(ArachneSignalHandler* pSigFunc)
{
  // On an Alpha, if you want to catch an unaligned access, you have to first tell
  // the operating system that unaligned accesses are to generate SIGBUG signals.

  #ifdef __alpha
  if ( System( "uac p sigbus" ) != 0 )
  { 
      cout << "Unable to catch unaligned accesses." << endl;
      exit(1);
  }
  #endif

  //  set up the sigaction function by defining members
  struct sigaction act, oact;

  act.sa_sigaction = pSigFunc;
  act.sa_flags = SA_SIGINFO;

  //  init mask to empty
  sigemptyset(&act.sa_mask);

  //  define signals to block
  sigaddset( &act.sa_mask, SIGSEGV );
  sigaddset( &act.sa_mask, SIGFPE  );
  sigaddset( &act.sa_mask, SIGINT  );
  sigaddset( &act.sa_mask, SIGBUS  );
  sigaddset( &act.sa_mask, SIGILL  );
  sigaddset( &act.sa_mask, SIGABRT  );
  sigaddset( &act.sa_mask, SIGTERM  );
  sigaddset( &act.sa_mask, SIGUSR1  );
  sigaddset( &act.sa_mask, SIGUSR2  );
  sigaddset( &act.sa_mask, SIGCHLD  );

  //  now send the signal to the handler
  if ( sigaction( SIGSEGV, &act, &oact ) < 0 )
  {
    cout << "Error catching SIGSEGV signal." << endl;
    exit(1);
  }
  if ( sigaction( SIGFPE, &act, &oact ) < 0 )
  {
    cout << "Error catching SIGFPE signal." << endl;
    exit(1);
  }
  if ( sigaction( SIGBUS, &act, &oact ) < 0 )
  {
    cout << "Error catching SIGBUS signal." << endl;
    exit(1);
  }
  if ( sigaction( SIGILL, &act, &oact ) < 0 )
  {
    cout << "Error catching SIGILL signal." << endl;
    exit(1);
  }
  if ( sigaction( SIGABRT, &act, &oact ) < 0 )
  {
    cout << "Error catching SIGABRT signal." << endl;
    exit(1);
  }
  if ( sigaction( SIGUSR1, &act, &oact ) < 0 )
  {
    cout << "Error catching SIGUSR1 signal." << endl;
    exit(1);
  }
  if ( sigaction( SIGUSR2, &act, &oact ) < 0 )
  {
    cout << "Error catching SIGUSR2 signal." << endl;
    exit(1);
  }
  if ( sigaction( SIGTERM, &act, &oact ) < 0 )
  {
    cout << "Error catching SIGTERM signal." << endl;
    exit(1);
  }
  if ( sigaction( SIGCHLD, &act, &oact ) < 0 )
  {
    cout << "Error catching SIGCHLD signal." << endl;
    exit(1);
  }
  act.sa_flags = SA_RESETHAND ^ SA_NOMASK;
  if ( sigaction( SIGINT, &act, &oact ) < 0 )
  {
    cout << "Error catching SIGINT signal." << endl;
    exit(1);
  }
}

// ===============================================================================
//
// RunTime: provide for various run-time services: turn off core dumps, make sure
// stack size is big enough, catch some interrupts.
//
// ===============================================================================

void RunTime( int no_dump, ArachneSignalHandler* pSigFunc )
{
  // Decouple C-style stdio and C++-style c{out,err,log} streams.

  ios::sync_with_stdio(false);

  // Turn off core dumps.

  if ( no_dump ) NoDump();

  // Set up to catch failures of new.

  std::set_new_handler(our_new_handler);

  // Make sure stack size is big enough.

  SetLimits(); // Don't even think about removing this.
               // Stack size too small will cause hideously untraceable errors.

  // Prepare to catch interrupts.

  if ( pSigFunc ) ArachneInterruptHandler(pSigFunc);

  // Ensure that shell is compliant with XPG4.  The only known affect of this
  // is to prevent procbuf closures from producing spurious blank lines.

  // Below is deprecated in gcc 4.2.  Must declare/alloc string 1st, then pass.
  // putenv( "BIN_SH=xpg4" );
  char tempenvstr[] = "BIN_SH=xpg4";
  putenv( tempenvstr );

  // Register the memory tracker's atexit function.  This function will have
  // no effect if no modules are doing memory tracking.

  atexit( &check_memory_tracker_records );
}
