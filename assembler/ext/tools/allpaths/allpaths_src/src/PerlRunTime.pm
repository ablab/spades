# PerlRunTime.pm
#
# A Perl module to provide automatic run-time signal handling, in much the same
# style that Arachne allows signal handling (implemented in system/RunTime.cc)
# to its C++ executables.
#
# Exported functions:
#
# run_or_die: A replacement for system(), with the following improvements:
# -- If the command returns nonzero, the parent script will abort.  This makes
#    run_or_die much better than system() for calling C++ modules within Perl
#    scripts; if you use system(), the user will be unable to stop the script
#    with Ctrl-C.
# -- You can call run_or_die with a second argument: a reference to a
#    "bail function", (e.g., a user-supplied exception handler) which will be
#    called in the event of command failure.
# -- You can call run_or_die in non-null context
#    (i.e., @output = run_or_die("command") instead of run_or_die("command"))
#    in which case it redirects stdout into the output variable you specify.
# -- Note that run_or_die is no more secure than system; you should not give it
#    any unchecked user data.
#
# Perl syntax examples:
# run_or_die ($command)
# $ouput = run_or_die($command)
# @output = run_or_die($command, \&bail_function)
#
#
# To activate this module's signal handling, just "use PerlRunTime" in your
# Perl script.  You may get an error ("Can't locate PerlRunTime.pm in @INC...")
# which you can fix by using the FindBin module as follows.  (You may need to
# adjust the second line, depending on the location of your script relative to
# PerlRunTime.pm.)
#
# use FindBin;
# use lib "$FindBin::Bin/../"
#
# Josh Burton
# April 2008 

package PerlRunTime;
use strict;
use Config; # gives us the localized conversion from signal numbers to names
use Carp; # allows 'confess', which is like 'die', but it prints a stack trace

use vars qw($VERSION @ISA @EXPORT @EXPORT_OK);
require Exporter;
@ISA = qw(Exporter AutoLoader);
@EXPORT = qw(run_or_die);


my @signals = split / /, $Config{'sig_name'}; # Convert from signal ID to name
my %message; # explanatory messages for each signal



# See top of module for documentation
sub run_or_die
{
    my ($command, $bail) = @_;
    my @output = ();
    
    
    # If run_or_die was called in null context, run the command with 'system'
    if (!defined wantarray) {
	system ($command);
    }
    # Otherwise, run it through 'open' to capture stdout and stderr
    else {
	open (OUTPUT, "$command 2>&1 |") or 
            confess("can't open pipe on command '$command'."); 
	@output = <OUTPUT>;
	close OUTPUT;
    }
    # The command's exit status is now in $?
    
    # build the return values, regardless of exit status
    my @ret = ();
    my $ret = "";
    if (defined wantarray) {
        if (wantarray) { @ret = @output; }          # expects a list
        else           { $ret = join '', @output; } # expects a scalar
    }

    # Run ok!  Return the data in the requested context
    return (@ret ? @ret : $ret) if ($? == 0);
    
    
    # Something went wrong!
    # There are 3 possible reasons for a process to have nonzero exit status:
    # an interrupt signal (e.g., a seg fault, or user hit Ctrl-C),
    # a return failure (e.g., an exit(-1)), or a core dump
    my $exit_value = $? >> 8;
    my $signal_num = $? & 127;
    my $core_dump  = $? & 128;
    
    $exit_value = -1 if ($exit_value == 255);
    my $signal = $signals[$signal_num];
    
    # See if something is actually wrong (some signals should be ignored)
    return (@ret ? @ret  : $ret) unless ($message{$signal} ||
                                         $exit_value ||
                                         $core_dump);
    
    
    # Call the bail function, if it was supplied
    if (ref $bail eq 'CODE') {
	&$bail;
	return (@ret ? @ret : $ret);
    }
    
    # Print error message
    my $time = localtime;
    print STDERR "\n$time - $0\n";
    print STDERR "A shell command ";
    
    if    ($signal_num) { print STDERR "detected signal SIG$signal ($message{$signal}).\n"; }
    elsif ($exit_value) { print STDERR "returned exit value $exit_value.\n"; }
    elsif ($core_dump)  { print STDERR "caused a core dump.\n"; }
    print STDERR "Exiting Perl.\n\n";
    
    print STDERR "Output was:\n\n@output\n\n" if (@output);
  
    print STDERR "Stack trace:\n";
    
    # 'confess' gives a thorough stack trace, then exits
    confess("\tShell command '$command'\n");
}






# Signal/interrupt handler subroutine
sub signal_handler ($ )
{
    my ($signal) = @_;
    
    # Print the current time and a useful message
    my $time = localtime;
    print STDERR "\n$time - $0\n";
    print STDERR "Detected signal SIG$signal ($message{$signal}).  Exiting Perl.\n\n";
    print STDERR "Stack trace:";
    
    # 'confess' gives a thorough stack trace, then exits
    confess("\n");
}








# Signal handling: this causes the subroutine "signal_handler" in this module
# to be called on each of these signals
# Signals not on this list will be ignored (as is standard for those signals)
$SIG{'HUP' } = \&PerlRunTime::signal_handler;
$SIG{'INT' } = \&PerlRunTime::signal_handler;
$SIG{'QUIT'} = \&PerlRunTime::signal_handler;
$SIG{'ILL' } = \&PerlRunTime::signal_handler;
$SIG{'TRAP'} = \&PerlRunTime::signal_handler;
$SIG{'ABRT'} = \&PerlRunTime::signal_handler;
$SIG{'BUS' } = \&PerlRunTime::signal_handler;
$SIG{'FPE' } = \&PerlRunTime::signal_handler;
$SIG{'SEGV'} = \&PerlRunTime::signal_handler;
$SIG{'TERM'} = \&PerlRunTime::signal_handler;
$SIG{'STOP'} = \&PerlRunTime::signal_handler;


# Explanatory messages for each handled signal (plus a few more)
%message =
    (
     'HUP'  => 'Terminal hangup',
     'INT'  => 'Interrupt received - perhaps a ctrl-c',
     'QUIT' => 'Quit - perhaps a ctrl-\\',
     'ILL'  => 'Illegal instruction',
     'TRAP' => 'Trace/breakpoint trap',
     'ABRT' => 'Abort',
     'BUS'  => 'Bus error',
     'FPE'  => 'Arithmetic or floating-point exception',
     'KILL' => 'Kill (unblockable)',
     'SEGV' => 'Segmentation fault',
     'PIPE' => 'Broken I/O pipe',
     'TERM' => 'Process terminated',
     'STKFLT' => 'Stack fault',
     'STOP' => 'Stopped signal',
     );



1;
