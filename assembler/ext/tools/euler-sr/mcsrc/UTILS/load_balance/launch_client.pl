#!/usr/bin/env perl

if ($#ARGV != 2) { 
  print "usage: $0 machines_file client_script load_balance_dir\n";
  exit(0);
}
$machinesFile = shift @ARGV;
$client = shift @ARGV;
$lbdir  = shift @ARGV;

my @machines = ();
my %limits = ();

# flush the loads for each machine in case clients are not running.
getMachines($machinesFile, \@machines, \%limits);


for $machine (@machines) {
  $command = "ssh -f $machine \"/home/mchaisso/UTILS/load_balance/lc.csh $client $lbdir $machinesFile\"";
  print "$command";
 system($command);
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
