#!/usr/bin/env perl

BEGIN {
		unshift(@INC, $ENV{"EUSRC"} . "/assembly/");
}

use RunCmd;

if ($#ARGV < 1) {
  print <<ENDusage;
usage: FixErrors.pl reads.fasta tupleSize
For additional options, use:
  Assemble.pl reads.fasta tupleSize -onlyFixErrors -script > my_script.sh
and edit the script.  Run
  fixErrors
with no parameters to get a list of options to customize in the script.
ENDusage
  exit(0);
}
$readsFile = shift @ARGV;
$tupleSize = shift @ARGV;



$EUSRC = $ENV{"EUSRC"};
$fixErrorsCmd = "$EUSRC/assembly/Assemble.pl $readsFile $tupleSize -onlyFixErrors";
RunCmd::RunCommand($fixErrorsCmd);
