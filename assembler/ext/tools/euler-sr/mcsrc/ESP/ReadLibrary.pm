package ReadLibrary;

sub GetFastaTitle {
  my ($fastaHeader) = @_;
  if ($fastaHeader =~ />(\S+)/) {
    return $1;
  }
  else {
    return 0;
  }
}

sub IsABIQuality {
  my ($title) = @_;
  if ($title !~ /CHROMAT_FILE/ or $title !~ /PHD_FILE/) {
    return 1;
  }
  else {
    return 0;
  }
}

sub ParseABIQuality {
  my ($title) = @_;
  $title =~ />(\S+)/;
  $name = $1;
  return $name;
}

sub ParseABITitle {
  my ($title) = @_;
#  if (IsABIQuality($title)) {
#    print "parsing quality: $title\n";
#    $name = ParseABIQuality($title);
#    return (0,0,0,$name);
#  }
  my $base = "";
  my $dir = 0;
  my $type = -1;
  my $name ="";
  if ($title =~ />/) {
    $title =~ />(\S+)/;
    $name = $1;
  }
  else {
    $name = $title;
  }
  if ($title =~ /S401149/) {
 #   print "matching: $title\n";
 # example: >S401149.PSM-100-M13F.ab1
    $title =~ /S401149\.(.{6,16})([FR])/;
    $base = $1;
    $dir  = $2;
    $type = 0;
  }
  # 000001228602A01.F.ab1
  elsif ($title =~ /0+(.{8,12})\.([FR])\.ab1/) {
    $base = $1;
    $dir  = $2;
    $type = 1;
  }
   # example 000000597258A01@13750R1.guf
  elsif ($title =~ /@/){
    if ($title =~ /000000(.{4,10})@.{4,10}([FR]).\.[bg]uf/) {
      $base = $1;
      $dir = $2;
      $type = 3;
#      print "got $base $dir $type\n";
    }
  }
  else {
    print "failed to parse title: $title\n";
    exit(0);
  }
  return ($base, $dir, $type, $name);
}

@MATE_SPANS = (5000,60000,40000);
@MATE_MEAN  = (2660,36000,40000);
@MATE_STDEV = (830, 4000, 4400);

sub ParseCloneName {
  my ($cloneName) = $_;
  # parse a name of the format:
  # S401149.PSM-999-M13F.ab1
  $cloneName =~ /.*PSM\-(\d+)\-(.*)\.ab1/;
  my $number = $1;
  my $well   = $2;
  return ($number, $well);
}


sub ReadTitleFile {
  my ($titleFileName, $names) = @_;
  open (NF, $titleFileName);
  # read the names
  @$names = ();
  while(<NF>) {
    $title = $_;
    if ($title =~ />/) {
      chomp($title);
      my ($base, $dir, $type, $name) = ReadLibrary::ParseABITitle($title);
      push @$names, $name;
    }
  }
}

sub ReadMatePairFile {
  my ($mpFileName, $mphash) = @_;
  open(MPFILE, $mpFileName) or die "cannot open matepair file $mpFileName\n";
  while (<MPFILE>) {
    $line = $_;
    chomp $line;
    my @vals = split(/\s+/, $line);
    $$mphash{$vals[0]} = $vals[1];
    $$mphash{$vals[1]} = $vals[0];
  }
  close MPFILE;
}


sub ReadMatePairIndexFile {
  my ($mpFileName, $mpindex, $mphash) = @_;
  open(MPFILE, $mpFileName) or die "cannot open $mpFileName\n";
  @$mpindex = <MPFILE>;
  chomp @$mpindex;
}

# current format of the locations file:
#579	1	580	2199729	2200308	0.0	1	0	580	1
$MLength = 0;
$MRefStart = 1;
$MRefEnd   = 2;
$MQryStart = 3;
$MQryEnd   = 4;
$MScore    = 5;
$MStrand   = 6;
$MMaskStart= 7;
$MMaskEnd  = 8;
$MPctIdent = 9;


sub ReadLocationsFile {
  my ($locationsFileName, $locations) = @_;
  open(LOC, $locationsFileName) or die "cannot open $locationsFileName\n";
  %locations = ();
  my ($base, $dir, $type, $name);
  while (<LOC>) {
    $line = $_;
    chomp $line;
    if ($line =~ />/) {
      ($base, $dir, $type, $name) = ParseABITitle($line);
    }
    else {
      $match = new ReadMatch($line);
      if (! exists $$locations{$name}) {
	$$locations{$name} = ();
      }
      push @{$$locations{$name}}, $match;
    }
  }
}

sub FindTrimStart {
  my ($seq) = @_;
  my $startTrimSeq = "";
  my $endTrimSeq = "";
  $startTrimSeq = "";
  $endTrimSeq = "";
  $seqStarted = 0;
  while ($seqStarted == 0) {
    if ($seq =~ /^(N+)[^N].*/) {
      $match = $1;
      $startTrimSeq .= $match;
      $seq = substr($seq, length $match);
      if ($seq =~ /^([^N]{1,5}N+)/) {
	$match = $1;
	$startTrimSeq .= $match;
	$seq = substr($seq, length $match);
      }
    } else {
      $seqStarted = 1;
    }
  }
  $seqEnded = 0;
  while ($seqEnded == 0) {
    if ((length $seq) == 0) {
      $seqEnded = 1;
    } else {
      if ($seq =~ /[^N]*(N+)$/) {
	$match = $1;
	$l = length $match;
	$endTrimSeq = $match . $endTrimSeq;
	$trim = length($seq) - $l;
	$seq = substr($seq, 0, $trim);
	if ($seq =~ /(N+[^N]{1,5}$)/) {
	  $match = $1;
	  $endTrimSeq = $match . $endTrimSeq;
	  $l = length $match;
	  $trim = length($seq) - $l;
	  $seq = substr($seq, 0, $trim);
	}
      } else {
	$seqEnded = 1;
      }
    }
  }
  my $startTrim = length $startTrimSeq;
  my $endTrim = length $endTrimSeq;
  return ($startTrim, $endTrim);
}

sub ReadNamesFile {
  my($readNamesFile, $reads, $names) = @_;
  %$reads = ();
  open(RN, "$readNamesFile") or die "cannot open $readsNameFile\n";
  $pos = 0;
  while (<RN>) {
    if (/>/) {
      my ($base, $strand, $type, $name) = ReadLibrary::ParseABITitle($_);
      $$reads{$name} = $pos;
      push @$names, $name;
      $pos++;
    }
  }
  close RN;
}




return 1;


