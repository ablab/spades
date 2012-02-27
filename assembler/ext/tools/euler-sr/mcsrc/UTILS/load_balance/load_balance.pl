#!/usr/bin/env perl
use strict;
use Socket;
use IO::Handle;
use POSIX;

if ($#ARGV < 1) {
  print "usage: $0 machines_file jobs_file [delay] [out_file]\n";
  exit(0);
}

# spawn a child process.
defined(my $pid = fork)   or die "Can't fork: $!";
exit if $pid;
setsid;



my $delay = 60;
if ($#ARGV >= 2) {
  $delay = @ARGV[2];
}
if ($#ARGV == 3) {
  my $out = @ARGV[3];
  open STDOUT, ">$out" or die "cannot open $out\n";
  open STDERR, ">$out.err" or die "cannot open $out.err\n";
}

open STDIN, '/dev/null'   or die "Can't read /dev/null: $!";

my $machinesfile = @ARGV[0];

# initialize for running in the background


my @machines = ();
my %limits = ();
my %loads =    ();
my %fhrd = ();
my %fhwr = ();
my %pids  =  ();
my $pid   = 0;
my %online= ();
my $line;
my $chfh;
my $status;
my @comandQueue = ();
my @available = ();
my @busy = ();
my @unavailable = ();

getMachines($machinesfile, \@machines, \%limits);
my $machine;

if (@machines < 1) {
  print "must specify at least 1 machine!!!\n";
  exit(0);
}

my ($machine, $avg);
my ($load, $numjobs);
$pid = 1;
my @loadAvgs = ();
print "getting load averages\n";
@loadAvgs = getLoadAvgs(@machines);
print "done with load averages\n";
my $index = 0;
foreach $load (@loadAvgs) {
  $machine = @machines[$index];
  if ($load eq "") {
    $online{@machines[$index]} = 0;
  }
  elsif ($limits{@machines[$index]} > $load) {
    $numjobs = POSIX::floor($limits{@machines[$index]} - $load);
    my $i;
    for $i (1 .. $numjobs) {
#      print "$machine is available\n";
      push @available, $machine;
      $online{@machines[$index]} = 1;
    }
  }
  else {
    push @busy, @machines[$index];
    $online{@machines[$index]} = 1;
  }
  $index = $index + 1;
}


my $unavailable;

# remove unavailable
foreach $machine (@machines) {
  if ($online{$machine} == 0) {
      pop @machines;
    }
}

# read what jobs to submit
open(JOBS, "@ARGV[1]");
my @jobs = ();
@jobs = <JOBS>;
chomp @jobs;

# variables to keep track of what job is running, result of jobs, etc.
my $res;
my $countJobs =0;
my $job;
my $togo;

while ($#jobs >= 0) {
  print "beginning of loop\n";
  # submit as many jobs as possible to available processors
  while ($#available >= 0 && $#jobs >= 0) {
    $countJobs += 1;
    $job        = pop @jobs;
    $machine    = pop @available;
    $res = system("ssh  $machine -q -f -t \"$job &\" & ");
  }
  print "submitted jobs\n";
  # Wait for loads to be updated on computers.
  sleep($delay);     # will change to 60 later.

  # query to see which machines are available now, for all machines
  # that can fit another job, add them to a queue. 
  $index = 0;
  @loadAvgs = getLoadAvgs(@machines);
  foreach $load (@loadAvgs) {
    $machine = @machines[$index];
    $index = $index + 1;
    if ($limits{$machine} > $load) {
      # $machine can fit at least one job, add it to the queue.
      $numjobs = POSIX::floor($limits{$machine} - $load);
      my $i;
      for $i (1 .. $numjobs) {
	push @available, $machine;
      }
    }
  }
  print "got available: @available\n";
  $togo = $#jobs+1;
  print "$countJobs submitted, $togo to go\n";
}
#print "$countJobs submitted \n";

sub getLoadAvgs {
  my @machines = @_;
  my @loads = ();
  my ($machine, $loadavg);
  foreach $machine (@machines) {
    $loadavg = getLoadAvg($machine);
    if ($loadavg >= 0) {
      push @loads, $loadavg;
    }
    else {
      push @loads, "9999";
    }
  }
  return @loads;
}
    
sub getMachines {
  my($machinesfile, $machines, $limits) = @_;
  @{$machines} = ();
  %{$limits} = ();
  print "opening machines file $machinesfile\n";
  open (MACHINES, "$machinesfile") or die "cannot open machines file\n"; 
  while (<MACHINES>) {
    chomp($line = $_);
  $line =~ /(\w+)\s+([\d\.]+)/;
    push @{$machines}, $1;
    # start with all machines available, will make busy after checking.
    $$limits{$1} = $2;
  }
  print "done with machines file\n";
  close MACHINES;
}

sub getLoadAvg {
  my ($machine) = @_;
  my $topstr = `ssh -f -t -q $machine \"top -n 1 \"`;
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
  
