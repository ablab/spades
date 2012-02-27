############################################################################
# Title:          RunCmd.pm
# Author:         Mark Chaisson, Glenn Tesler
# Created:        2008
# Last modified:  02/23/2010
# 
# Copyright (c) 2008-2010 The Regents of the University of California
# All Rights Reserved
# See file LICENSE for details.
############################################################################

package RunCmd;

use POSIX;

sub RunCommand {
  my ($cmd) = @_;
	print "$cmd\n";
	print LOG "running $cmd\n";
  $output = `$cmd`;
  $status = $?;
  print LOG "$status $cmd\n";
  print LOG "output:\n";
  print LOG "$output";
  if ($status != 0) {
    $host = $ENV{"HOST"};
    $dir  = $workDir;
    print RESULT "$jobID FAILED, $status at $host, $dir, command $cmd : $argstr\n";
    close LOG;
    close RESULT;
    exit(0);
  }
}

###############################################################################
# $output = RunCommandLog($cmd, $desc, $verbose, $capture)
#
# $cmd: shell command to execute
# $desc: short description (couple of words)
# $verbose:
#    -1: just print commands; don't execute them.
#     0: don't display output from command
#     1: display output from command
# $capture:
#     0: don't return output from command (may run faster), just return undef
#     1: return output from command
#
# If $cmd fails, print an error message and die.
#
#
# Global variable (for now):
#   $RunCmd::reportfile --- if defined, failed $cmd will also print
#   message to this file
#
#
# Possible TODOs:
# * Existing RunCommand has a number of global vars as input; deal with them.  
# * Send output to a log file instead of STDOUT
# * Capture STDERR too
# * Verbosity levels
# * Severe level (for nonzero return status)
###############################################################################

# For STDOUT outpu:
$FormatKey = "\%-16s ";
$FormatLine = "\%-16s \%s\n";
$LineWidth = 79;

# For .report file:
$FormatLineReport = $FormatLine;

# $RunCmd::reportfile --- global variable for now

sub RunCommandLog {
		my ($cmd, $desc, $verbose, $capture) = @_;

		if ($verbose == -1) {
				print $cmd, "\n";
				return;
		}

		print '#' x $LineWidth;
		print "\n";
		printf $FormatLine, "Task:", $desc;
		printf $FormatLine, "Command:", $cmd;
		printf $FormatLine, "Start:", CurTimeString();

		my $possibly_nested = ($cmd =~ /\.pl/);

		if ($verbose) {
				print PadDashLine($possibly_nested ? "BEGIN OUTPUT [$desc]" : "BEGIN OUTPUT"), "\n";
		}

		my ($out,$errmsg) = RunCommandPlain($cmd, $verbose, $capture);
		
		if ($verbose) {
				print PadDashLine($possibly_nested ? "END OUTPUT [$desc]" : "END OUTPUT"), "\n";
		}

		my $endstr = $errmsg ? "Abort" : "End";
		if ($errmsg) {
				print "\n", AlertLine(), "\n\n";
		}
		if ($possibly_nested) {
				printf $FormatLine, "${endstr} [${desc}]:", CurTimeString();
		} else {
				printf $FormatLine, "${endstr}:", CurTimeString();
		}

		if ($errmsg) {
				printf $FormatLine, "Error:", $errmsg;
				printf $FormatLine, "", "$desc failed";
				printf $FormatLine, "Task:", $desc;
				printf $FormatLine, "Command:", $cmd;

				if ($RunCmd::reportfile
						&& open(REPORT,">>",$RunCmd::reportfile)) {
						print REPORT "\n", AlertLine(), "\n\n";
						printf REPORT $FormatLineReport, "Aborted:", CurTimeString();
						printf REPORT $FormatLineReport, "Command:", $cmd;
						printf REPORT $FormatLineReport, "Error:", $errmsg;
						print REPORT "\n", PadDashLine(), "\n\n";
						close(REPORT);
				}

#				exit(0);
		}
		return $out;
}

sub RunCommandPlain {
		my ($cmd, $verbose, $capture) = @_;

		my $out;       # string to return
		my $err;       # error flag to return
		my $status;    # system status

		if ($capture) {
				$out = `$cmd`;
				$status = $?;

				if ($verbose) {
						print $out;
						if (length($out) > 0
								&& substr($out,-1,1) ne "\n") {
								print "\n[Output did not end with a newline.]\n";
						}
				}
		} else {
				if ($verbose) {
						system($cmd);
						$status = $?;
				} else {
						system("( $cmd ) > /dev/null");
						$status = $?;
				}
		}

#		print "status = $status\n";

		my $errmsg = undef;
		if ($status == -1) {
				$errmsg = "failed to execute: $!";
				$err = 1;
		} elsif ($status & 127) {
				$errmsg =
						sprintf "process died with signal %d, %s coredump",
						($status & 127),  ($status & 128) ? 'with' : 'without';
				$err = 1;
		}	else {
				$err = $status >> 8;
				if ($err) {
						$errmsg =
								"process exited with nonzero status $err";
				}
		}
#		if ($errmsg) {
#				print $errmsg, "\n";
#		}

		return ($out,$errmsg);
}

# OnlyLogCommand($cmd, $desc, $reportSwitch, $stdoutSwitch)
#   Generate log entries w/o running command
#   Useful for commands that don't do log entries themselves
#      $reportSwitch true: log to .report file
#      $stdoutSwitch true: log to STDOUT
sub OnlyLogCommand {
		my ($cmd, $desc, $reportSwitch, $stdoutSwitch) = @_;

		if ($stdoutSwitch) {
				print '#' x $LineWidth;
				print "\n";
				printf $FormatLine, "Task:", $desc;
				printf $FormatLine, "Command:", $cmd;
				printf $FormatLineReport, "Directory:", getcwd();
				printf $FormatLine, "Start:", CurTimeString();
		}

		if ($reportSwitch && $RunCmd::reportfile
				&& open(REPORT,">>",$RunCmd::reportfile)) {
#				printf REPORT $FormatLineReport, "Running:", $0;
				printf REPORT $FormatLineReport, "Directory:", getcwd();
				printf REPORT $FormatLineReport, "Start:", CurTimeString();
				printf REPORT $FormatLineReport, "Command:", $cmd;
				print REPORT "\n", PadDashLine(), "\n\n";
				close(REPORT);
		}
}


sub CurTimeString {
		return strftime("%a %b %e %H:%M:%S %Z %Y",localtime());
}

sub AlertLine {
		my $rep = "* ";
		my $s = $rep x int($LineWidth/2);
		my $rem = $LineWidth - length($s);
		if ($rem > 0) {
				$s .= substr($rep, 0, $rem);
		}
		return $s;
}

sub PadDashLine {
		my ($s) = @_;
		$s =~ s/\s+$//;
		if (length($s) < $LineWidth) {
				if ($s ne '') {
						$s .= " ";
				}
				my $n = $LineWidth - length($s);
				if ($n > 0) {
						$s .= '-' x $n;
				}
		}
		return $s;
}


return 1;
