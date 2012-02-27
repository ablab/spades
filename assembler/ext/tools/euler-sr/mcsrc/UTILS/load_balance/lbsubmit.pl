#!/usr/bin/env perl
use Fcntl ':flock';


$usage = "usage: $0 commands_file  sever_commands_file\n";


$commands = shift @ARGV;
if ($#ARGV >= 0) {
  $serverCommands = shift @ARGV or die "$usage\n";
}
else {
  if (exists($ENV{"LBCOMMANDS"})) {
    $serverCommands = $ENV{"LBCOMMANDS"};
  }
  else {
    print "must either specify server commands file, or set in environment LBCOMMANDS\n";
    exit;
  }
}

(-e $serverCommands ) or `touch $serverCommands`;

open (COMMANDS, "$commands") or die "Cannot open $commands\n";
@commands = <COMMANDS>;
chomp @commands;

open (SERVER_COMMANDS, ">>$serverCommands") or die "cannot open $serverCommands\n";
flock(SERVER_COMMANDS, LOCK_EX);
$d = `date`;
chomp $d;
foreach $c (@commands) {
  print SERVER_COMMANDS "$d#$c\n";
}

