#!/usr/bin/env perl

use Fcntl ':flock';
use POSIX;

if (!exists($ENV{"HOST"})) {
  print "set environment variable HOST to the name of this machine\n";
  exit(0);
}
$host = $ENV{"HOST"};

$usage = "usage: $0 load_dir machines_file\n ".
   "load_dir:      directory that will have $host.avail and $host.jobs \n".
   "machines_file: that should contain a line with $host load, where load \n".
   "               is the maximum load for $host\n. ";


print "starting client on $host\n";
$loadDir        = shift @ARGV;
$machinesFile   = shift @ARGV or die "$usage\n";
$delay    = shift @ARGV or ($delay = 30);

$availFile      = $loadDir . "/$host.avail";
$jobsFile       = $loadDir . "/$host.jobs";
# create the daemon.



chdir '/'                 or die "Can't chdir to /: $!";
umask 0;
#open STDIN, '/dev/null'   or die "Can't read /dev/null: $!";
#open STDOUT, '>/dev/null' or die "Can't write to /dev/null: $!";
#open STDERR, '>/dev/null' or die "Can't write to /dev/null: $!";
defined(my $pid = fork)   or die "Can't fork: $!";
if ($pid) {
  print "main is exiting\n";
  exit;
}
setsid                    or die "Can't start a new session: $!";


while (1) {
  open (MACHINES, "$machinesFile") or die "cannot open $machinesFile\n";
  # share the machines file with other processes that are reading it.
  flock(MACHINES, LOCK_SH);
  $loadLimit = getLoadLimit(\*MACHINES, $host);
  close MACHINES; # unlocks the file 

  # read jobs that are given to this machine
  # consume a jobs file
  ( -e $jobsFile)  or `touch $jobsFile`;
  open(JOBS, "$jobsFile") or "die cannot open $jobsFile\n";
  @jobs = <JOBS>;
  chomp @jobs;
  # reset the jobs file
  close JOBS;
  `rm -f $jobsFile`;

  $numSubmitted = 0;
  while ($#jobs >= 0) {
    $job = shift @jobs;
    if ($job ne "") {
      print "$host started $job\n";
      system("$job &");
      $numSubmitted++;
    }
  }

  $loadAvg   = getLoadAvg();

  # produce availability file with the number of jobs this machine can take

  # total load = current load + numSubmitted (assume will generate a load)
  # jobs can take = total limit - total load
  $availJobs = POSIX::floor($loadLimit - $loadAvg - $numSubmitted);
  open (AVAIL, ">$availFile") or die "cannot open $availFile\n";
  flock(AVAIL, LOCK_EX);
  print AVAIL "$availJobs\n";
  close AVAIL;

  sleep ($delay);
}



sub getLoadLimit {
  my ($machFh, $host) = @_;
  $loadLimit = 0;
  while (<$machFh>) {
    $line = $_;
    if ($line =~ /$host\s([\d\.]+)/) {
      $loadLimit = $1;
    }
  }
  return $loadLimit;
}

sub getLoadAvg {
  $os = `uname`;
  chomp $os;
  if ($os eq "FreeBSD") {
    $la = `getloadavg.bsd`;
  }
  else {
    $la = `getloadavg`;
  }
  $la =~ /([\d\.]+)\s+.*/;
  return $1;
}
