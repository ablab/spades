#!/usr/bin/env perl
use Getopt::Long;

(exists $ENV{"HOME"}) or die "HOME environment variable not set\n";
$home = $ENV{"HOME"};
%options = ();
@emails = ();
@commands = ();
@dirs = ();
$priority = 0;

if (! -e "$home/.notify") {
  print "UPDATE: this now writes temporary files to a directory ~/.notify for sending e-mails.\n";
  print "Creating this directory now\n";
  system("mkdir $home/.notify");
}

if (GetOptions('directory=s'=> \@dirs, 
	       'email=s'=> \@emails, 
	       'commands=s'=>\@commands, 
	       'priority'=>\$priority, 
	       'startup=s'=>\@startups) == false || $#commands < 0) {
  print "usage: $0 -d dir [-d dir2...] -e e-mail [-e email2 ...] [-s startup e-mail]  -c \"command\" [-c \"command2\"...]\n";
  print "       $0, a script to run jobs and let you know when they are done.\n";
  print "       The first dir specified corresponds to the first command, etc.\n";
  print "       Priority is ignored.\n";
  print "       -e address: sends an e-mail to 'address' when the job is complete.\n";
  print "       -s address: sends and e-mail to 'address' when the job starts.\n";
  print "       -c \"command\":  a command to run (the quotes around the command are required).\n";
  print "       -d dir:     run the commands in directory 'dir'.\n";
  print " Note: make sure all paths are full (don't use ~)\n";
  exit(0);
}

if ($#startups >= 0) {
  $startAddresses = join(' ', @startups);
}

print "sending to @emails\n";

if ($#emails >= 0) {
  $from = $ENV{"USER"};
  $out = "/home/$from/.notify/notify.$$";
  open (OUT, ">$out");
  print OUT "~s notify.pl results\n";
  print OUT "Message from $from\n";
  print "got commands: @commands\n";
}
foreach $c (@commands) {
  $startTime = `date`;
  if ($#startups >= 0) {
    $startFile = "$home/.notify/start.$$\n";
    open (START, ">$startFile");
    print START "~s job started\n";
    print START "A message from notify.pl\n";
    print START "Started job: \"$c\"\n";
    print START "Start time: $startTime\n";
    close START;
    `mailx $startAddresses < $startFile`;
  }

  if ($#dirs >= 0) {
    chdir(shift @dirs);
  }
  $res = `$c`;
  $endTime  = `date`;
  $status = $?;
  if ($#emails >= 0) {
    print OUT "job   : $c\n";
    print OUT "status: $status\n"; 
    print OUT "start : $startTime"; 
    print OUT "end   : $endTime\n";
    print OUT "res:  $res\n";
  }
  if ($#startups > 0) {
    print "removing `$startFile' \n";
    `rm $startFile`;
  }
}


if ($#emails >= 0) {
  close(OUT);
  if ($#emails >= 0) {
    $addresses = join(' ', @emails);
    `mailx $addresses < $out`;
  }
  `rm $out`;
}
