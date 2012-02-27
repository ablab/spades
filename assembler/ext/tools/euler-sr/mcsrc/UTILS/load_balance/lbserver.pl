#!/usr/bin/env perl
use Socket;
use IO::Handle;
use POSIX;
use Fcntl ':flock';

$usage =  "usage: $0 loads_dir machines_file jobs_file [delay] [out_file]\n";


$loadsDir = shift @ARGV;
$machinesFile = shift @ARGV;
$jobsFile = shift @ARGV or die "$usage\n";
$delay   = shift @ARGV or ($delay = 10);

# spawn a child process.
defined(my $pid = fork)   or die "Can't fork: $!";
exit if $pid;
setsid;



if ($#ARGV >= 0) {
  # redirect output if told to
  $out = shift @ARGV;
  open STDOUT, ">$out" or die "cannot open $out\n";
  open STDERR, ">$out.err" or die "cannot open $out.err\n";
}

open STDIN, '/dev/null'   or die "Can't read /dev/null: $!";

my @machines = ();
my %limits = ();

# flush the loads for each machine in case clients are not running.
getMachines($machinesFile, \@machines, \%limits);

foreach $host (@machines) {
  $availableFile = $loadsDir . "/$host.avail";
  `rm $availableFile`;
}


while (1) {
  # get the jobs;
  open (JOBS, "$jobsFile") or die "cannot open $jobsFile\n";
  flock(JOBS, LOCK_EX);
  @jobs = <JOBS>;
  getMachines($machinesFile, \@machines, \%limits);
  # cycle through machines submitting jobs
  $machIdx = 0;
  while ($machIdx <= $#machines and $#jobs >= 0) {
    # get the load for the machine.
    $availableFile = $loadsDir . "/@machines[$machIdx].avail";

    # consume availability  file
    (-e $availableFile) or `touch $availableFile`; 
    open(LOAD, $availableFile) or die "cannot open $load File\n";
    flock(LOAD, LOCK_EX);
    $limit = <LOAD>;
    if ($limit eq "") {
      $limit = 0;
    }
    close LOAD;
    `rm $availableFile`;

    if ($limit > 0) {
      # produce a jobs file
      # write as many jobs to the jobs file as the load permits
      (! -e $submitFile) and `touch $submitFile`;
      $l = 1;
      $submitFile  = $loadsDir . "/@machines[$machIdx].jobs";
      open (SUBMIT, ">>$submitFile") or die "cannot open $submitFile\n";
      #    flock(SUBMIT, LOCK_EX);
      while ($l <= $limit and $#jobs >= 0) {
	$jobline = shift @jobs;
	chomp $jobline;
	# strip the date
	$jobline =~ /[^\#]+\#(.*)/;
	$job = $1;
	print SUBMIT "$job\n";
	$l = $l + 1;
      }
      #    flock(SUBMIT, LOCK_UN);
      close SUBMIT;
      # write back jobs that haven't been processed yet.
    }
    $machIdx = $machIdx + 1;
  }
  close JOBS;
  # not thread safe, but oh well
  open (JOBS, ">$jobsFile") or die "cannot open $jobsFile\n";
  flock(JOBS, LOCK_EX);
  for $job (@jobs) {
    print JOBS "$job";
  }
  close JOBS; # release the lock on jobs
  sleep($delay);
}



sub getMachines {
  my($machinesfile, $machines, $limits) = @_;
  @{$machines} = ();
  %{$limits} = ();
  open (MACHINES, "$machinesfile") or die "cannot open machines file\n"; 
  while (<MACHINES>) {
    chomp($line = $_);
    $line =~ /(\S+)\s+([\d\.]+)/;
    push @{$machines}, $1;
    # start with all machines available, will make busy after checking.
    $$limits{$1} = $2;
  }
  close MACHINES;
}
