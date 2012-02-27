package FwgLib;

sub InitDir {
  my ($dir) = @_;
	if (! -e $dir) {
    `mkdir -p $dir`;
  }
}

sub GetBase {
		my ($file) = @_;
		if ($file =~ /\./) {
				$file =~ /^(.*)\.[^\.]+$/;
				return $1;
		}
		else {
				return $file;
		}
}

sub RemovePath {
  my ($file) = @_;
  if ($file =~ /\/([^\/]+)$/) {
    $base = $1;
    return $base;
  }
  else {
    return $file;
  }
}

sub CrucialGetEnv {
  my ($envName) = @_;
  if (!exists $ENV{$envName}) {
    print "ERROR! You must set the environment variable $envName to run $0\n";
    exit(1);
  }
  return $ENV{$envName};
}


sub SetupWorkingDir {

  my ($dir, $inputFile) = @_;
	if (! -e $inputFile) {
		print "ERROR, could not find input $inputFile\n";
		exit(1);
	}
	my $base = $inputFile;
  if (! -e $dir) {
    `mkdir $dir`;
  } 
	my $prevDir = `pwd`;
	chomp $prevDir;
	if (!chdir($dir)) {
    print "Could not create or go to direcory requested: $dir\n";
    exit(1);
  }
	my $base = RemovePath($inputFile);
	if (! -e "./$base") {
		`ln -sf $prevDir/$inputFile ./$base`;
	}
  return $base;
}


sub WaitOnJobIds {
  my @queryIds = @_;
  $delay = 20;
  $done = 0;
  #hashify queryIds
  while ($done == 0) {
    @curCmds = `qstat -u fwg.mchaisso`;
		# the first wo lines are junk
    shift @curCmds;
    shift @curCmds;
    # the rest of the lines have useful information
    @curCmdIds = ();
    grep(/(\d+)\s\S.*/ && push(@curCmdIds, $1), @curCmds);
#    print "got running ids: @curCmdIds\n";
#		print "waiting on @queryIds\n";
    %cidh = ();
    foreach $cid (@curCmdIds) {
      if ($cid ne "") {
				$cidh{$cid} = 1;
      }
    }
    $done = 1;
    foreach $qid (@queryIds) {
      if (exists $cidh{$qid}) {
				$done = 0;
      }
    }
    if ($done == 0) {
      sleep($delay);
    }
  }
}

sub SubmitCommandFile {
  my ($cmdFile, $dir, $njobs, $nproc) = @_;
  system("cd $dir; jobify $cmdFile $njobs -nproc 2");
  $cmdFile =~ /(.*)\.([^\.]*)$/;
  $base = $1;
  @commandIds = ();
  @commandIds = SubmitCommands($dir, $base, $njobs);
	
  return @commandIds;
}
  

sub SubmitCommands {
  my ($dir, $base) = @_;
  @commandList = `submit_jobs $dir/$base`;
  @jobIds = ();
  grep(/Your job (\d+)/ && push(@jobIds, $1), @commandList);
  return @jobIds;
}

sub SMPRunCommands {
  my ($dir, $nproc, $commandFile) = @_;
	$curdir = `pwd`;
	print "curdir: $curdir\n";
	print "smp running $dir $commandFile\n";
  $res = `cd $dir; jobify $commandFile 1 -nproc $nproc`;
	$EUSRC = $ENV{"EUSRC"};
	$base = GetBase($commandFile);
	$command = "$base.0.csh";
  `chmod +x $command`;
  # this will wait for the command file to finsih
  system("$dir/$command");
}

sub SubmitCommandListAndWait {
  my ($dir, $njobs, $nproc, $smp, @commands) = @_;
  $commandName = "cmd.". $$ . ".txt";
  open(CMDS, ">$commandName");
  chomp $commands;
  for ($c = 0; $c <= $#commands; $c++ ) {
    print CMDS "$commands[$c]\n";
  }
  if ($smp > 0) {
    SMPRunCommands($dir, $nproc, $commandName);
  }
  else {
    SubmitCommandsAndWait($commandName, $dir, $njobs, $nproc);
	}
	`cd $dir; rm cmd.$$.*.{pl,csh}`;
}
sub SubmitCommandsAndWait {
  my ($cmdFile, $dir, $njobs) = @_;
  @jobIds = SubmitCommandFile($cmdFile, $dir, $njobs, $nproc);
  WaitOnJobIds(@jobIds);
}
    
sub PrintLog {
    my ($logFile, $logStr) = $_;
		if ($logFile ne "") {
				print LOG "$logStr\n";
		}
}

sub RunLoggedCommand {
		my ($cmd, $logFile) = @_;
		system("echo $cmd >> $logFile");
		RunCommand($cmd);
}

sub RunCommand {
		my($cmd) = @_;
		print "$cmd\n";
		system($cmd);
		$res = $?;
		if ($res != 0) {
				print "FAILED: $cmd\n";
				exit(0);
		}
}

return 1;
