#!/usr/bin/env perl

use Fcntl ':flock';
use POSIX;

exists($ENV{"HOST"}
  print "set environment variable HOST to the name of this machine\n";
  exit(0);
}
$host = $ENV{"HOST"}

$usage = "usage: $0 load_dir machines_file\n"
  " load_dir: directory that will have $HOST.avail and $HOST.jobs \n"
  " machines_file: that should contain a line with $HOST load, where load is the maximum load for $HOST\n";



$loadDir        = shift @ARGV;
$machinesFile   = shift @ARGV or die "$usage\n";


$availFile      = $loadDir . "/$HOST.avail";
$jobsFile       = $loadDir . "/$HOST.jobs";
# create the daemon.



chdir '/'                 or die "Can't chdir to /: $!";
umask 0;
open STDIN, '/dev/null'   or die "Can't read /dev/null: $!";
#open STDOUT, '>/dev/null' or die "Can't write to /dev/null: $!";
open STDERR, '>/dev/null' or die "Can't write to /dev/null: $!";
defined(my $pid = fork)   or die "Can't fork: $!";
exit if $pid;
setsid                    or die "Can't start a new session: $!";


while (1) {
  open (MACHINES, "$machinesFile") or die "cannot open $machinesFile\n";
  # share the machines file with other processes that are reading it.
  flock(MACHINES, LOCK_SH);
  $loadLimit = getLoadLimit(\*MACHINES, $host);
  close MACHINES; # unlocks the file 
  $loadAvg   = getLoadAvg();
  
  # post how many jobs this machine can take
  $availJobs = POSIX::floor($loadLimit - $loadAvg);
  open (AVAIL, ">$availFile") or die "cannot open $availFile\n";
  flock(AVAIL, LOCK_EX);

  # read jobs that are given to this machine
  open (JOBS, ">$jobsFile") or die "cannot open $jobsFile\n";
  flock(JOBS, LOCK_EX);
  @jobs = <JOBS>;
  # reset the jobs file
  print JOBS "";  
  close JOBS;
  flock(JOBS, LOCK_UN); # just in case removing the file didn't unlock it. 
  
  foreach $job (@jobs) {
    `job &`;
  }
  
  $availJobs = $availJobs - ($#jobs + 1);
  print "$host wants $availJobs\n";
  print AVAIL "$availJobs\n";
  close AVAIL;

  sleep ($delay);
}



sub getLoadLimit {
  my ($machFh, $host) = @_;
  $loadLimit = 0;
  while (<$machFh>) {
    $line = $_;
    if ($line =~ /$host\s(\d+)/) {
      $loadLimit = $1;
    }
  }
  return $loadLimit;
}



sub getLoadAvg {
  my $topstr = `top -n 1`;
  my @topstr = split /\n/, $topstr;
  my $s;
  foreach $s (@topstr) {
    if ($s =~ /load averages:  (\d+\.\d+)/) { 
      return $1;
    }
    elsif ($s =~ /load average: (\d+\.\d+)/) {
      return $1;
    }
  }
  return -1;
}
  
