package Simulate;


sub RunCmd {
  my ($cmd) = @_;
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


return 1;
