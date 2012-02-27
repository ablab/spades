use DBI;
use Bio::Tools::Run::StandAloneBlast;
package common;


sub CreateTempFileName {
  my ($ext) = @_;
  $name = $$ . "ext";
  $test = $name . ".test";
  if (open(TMPF, ">/state/partition1/$test")) {
    close TMPF;
    system("rm /state/partition1/$test");
    return "/state/partition1/$name";
  }
  else {
    return $name;
  }
}

sub max {
  my ($a, $b) = @_;
  if ($a > $b) { return $a; }
  return $b;
}

sub min {
  my ($a, $b) = @_;
    if ($a < $b) { return $a;}
  return $b;
}

sub ReadLocusFile {
  my ($locusFileName, $loci) = @_;
  open(LF, $locusFileName) or die "cannot open $locusFileName\n";
  @locilines = <LF>;
  foreach $locusLine (@locilines) {
    push @$loci, [(split(/\s+/, $locusLine))];
  }
}

sub ReadIndexFile {
  my($indexFileName, $indices) = @_;
  open (IF, $indexFileName) or die "cannot open $indexFileName\n";
  while(<IF>) {
    $line = $_;
    $line =~ /(\S+)\s+(\S+)/;
    ${$indices}{$2} = $1;
    #print "($2) =  '$1'\n";
  }
}

sub MakeConnectString {
  $dbConnect = "dbi:mysql:$dbName";
  $importcommand = "";
  if ($scope eq "global") {
    if (exists $ENV{"GSCK"}) {
      $socket = $ENV{"GSCK"};
      print "using global socket: $socket\n";
      $dbConnect .= ";mysql_socket=$socket";
    }
  }
  else {
    if (exists $ENV{"LSCK"}) {
      $socket = $ENV{"LSCK"};
      $dbConnect .= ";mysql_socket=$socket";
    }
  }
  return $dbConnect;
}


sub ParseNCBIName {
  my ($name) = @_;
  @vals = split(/\|/, $name);
  return @vals;
}


sub IsRef {
  my ($name, $speciesList) = @_;
  my %foundSpecHash = ();
  $nspec = $#$speciesList;
  foreach $spec (@$speciesList) {
    if ($name =~/$spec/) {
      $foundSpecHash{$spec} = 1;
    }
  }

  my $specA = @foundSpec[0];
  my $specB = @foundSpec[1];

  if ($name =~ /$specB.*$specA/) {
    my $tmpSpec = $specA;
    $specA = $specB;
    $specB = $tmpSpec;
  }

  if ($name =~ /$specA\_$specB\_(\d+)_(\d+)/) {
    return 1;
  }
  else {
    $name =~ /$specA\_(\d+)\_(\d+)_$specB/;
    @res = ($specA, $specB, $1, $2);
    return 0;
  }
}

sub ParseName {
  my ($name, $speciesList) = @_;
  @res = ();
  if ($#$speciesList == 0) {
    if ($name =~ /(\w+)_(\d+)_(\d+)_(\w+)/) {
      print "simple ref parsed \n";
      @res = ($1,$4,$2,$3);
    }
    elsif ($name =~ /(\w+)_(\w+)_(\d+)_(\d+)/) {
      print "simple query parsed\n";
      @res = ($1,$2,$3,$4);
    }
  }
  else {
#    print "looking for species in pn in $#$speciesList\n";
    my %foundSpecHash = ();
    $nspec = $#$speciesList;
    foreach $spec (@$speciesList) {
#      print "checking $name for $spec\n";
      if ($name =~/\A$spec\_/ or $name=~/\d+\_$spec\Z/ or $name=~/[^\d]\_$spec\_\d/) {
	$foundSpecHash{$spec} = 1;
#	print "$found $spec\n";
      }
    }
    my @foundSpec = keys %foundSpecHash;
    my $s1, $s2;
    my $spec1, $spec2;
    $s1 = 0;
    while ($s1 < $#foundSpec) {
      $s2 = $s1 + 1;
      $spec1 = $foundSpec[$s1];
      while ($s2 <= $#foundSpec) {
	$spec2 = $foundSpec[$s2];
	if ($spec1 =~ /$spec2/) {
#	  print "$spec1 contians $spec2\n";
	  splice @foundSpec, $s2, 1;
	}
	else {
	  if ($spec2 =~ /$spec1/) {
#	    print "$spec2 continas $spec1\n";
	    splice @foundSpec, $s1, 1;
	    $s1--;
	  }
	  $s2++;
	}
      }
      $s1++;
    }
	
#    for ($s1 = 0; $s1 <= $#foundSpec; $s1++ ) {
#      print "found spec: @foundSpec[$s1]\n";
#    }
    if ($#foundSpec == 0) {
      # handle (corner) condition of only one name
      $spec = @foundSpec[0];
      if ($name =~ /$spec\_$spec\_(\d+)_(\d+)/) {
	@res = ($spec, $spec, $1, $2);
      }
      else {
	$name =~ /$spec\_(\d+)_$spec\_(\d+)/;
	@res = ($spec, $spec, $1, $2);
      }
    }
    else {
      my $specA = @foundSpec[0];
      my $specB = @foundSpec[1];
      # determine the order of the species
      if ($name =~ /$specB.*$specA/) {
	my $tmpSpec = $specA;
	$specA = $specB;
	$specB = $tmpSpec;
      }
#      print "trying to match $name with $specA $specB\n";
      if ($name =~ /$specA\_$specB\_(\d+)_(\d+)/) {
	@res = ($specA, $specB, $1, $2);
#	print "got res $specA $specB $1 $2\n";
      }
      else {
	$name =~ /$specA\_(\d+)\_(\d+)_$specB/;
	@res = ($specA, $specB, $1, $2);
#	print "got qry res $specA $specB $1 $2\n";
      }
    }
  }
  return @res;
}

sub max {
  my ($a, $b) = @_;
  if ($a > $b) {
    return $a;
  }
  else {
    return $b;
  }
}

sub min {
  my ($a, $b) = @_;
  if ($a < $b) {
    return $a;
  }
  else {
    return $b;
  }
}


sub GetHitCoordinates {
  my ($qrySeq, $database, $eThresh, $hits) = @_;
  @$hits = ();
  my $blastReport;
  $blastQuery =  Bio::Tools::Run::StandAloneBlast->new('program'=>'blastn',
						       'database'=>$database,
						       'q'=>'-2',
						       'W'=>'10',
						       'e'=>"$eThresh",
							);
  $blastReport = $blastQuery->blastall($qrySeq);

  my ($rb, $re, $qb, $qe, $ev);
  $rb = -1; $re = -1; $qb = -1; $qe = -1;
  while ( $result = $blastReport->next_result ) {
    while ($hit = $result->next_hit ) {
      while ($hsp = $hit->next_hsp ) {
	($rb, $re) = $hsp->range('query');
	($qb, $qe) = $hsp->range('sbjct');
	$ev = $hsp->evalue();
	push @$hits, [($rb, $re, $qb, $qe, $ev)];
      }
    }
  }
  return scalar @$hits;
}

sub LocateSequence {
  my ($qrySeq, $sbjct, $eThresh) = @_;
  my $hsp = GetFirstHSP(\$blastRes);
  my ($qryStart, $qryEnd, $sbjctStart, $sbjctEnd) = (0,0,0,0);
  if ($hsp) {
    my ($start, $end);
    ($qryStart, $qryEnd) = hsp->range('query');
    ($sbjctStart, $sbjctEnd) = hsp->range('sbjct');
  }
  return ($qryStart, $qryEnd, $sbjctStart, $sbjctEnd);
}

sub CheckForDeletion {
  my $eThresh = 0.001;

  my ($refSeq, $qrySeqDB, $refStart, $refEnd, $window, $skip);
  ($refSeq, $qrySeqDB, $refStart, $refEnd, $window, $skip, $eThresh) = @_;

  # the base of the dbname should be the base of the query sequence
  my $blastQuery = Bio::Tools::Run::StandAloneBlast->new('program'=>'blastn',
							 'database'=>$qrySeqDB,
							 'q'=>'-2',
							 'W'=>'10',
							 'e'=>"$eThresh",
							);

#  print "blasting from $qrySeqDB\n";
  my ($upstrStart, $upstrEnd);
  my ($downstrStart, $downstrEnd);

  # extract the flanking regions 
  $upstrStart = max(0, $refStart - $window - $skip);
  $upstrEnd   = max(0, $refStart - $skip);
  $len = $refSeq->length();
#  print "creating substr from $upstrStart to $upstrEnd should be less than $len\n";
  my $refUpstrSeq = Bio::Seq->new(-seq=>$refSeq->subseq($upstrStart, $upstrEnd),
				  -display_id=>"upstream");

  my $refSeqLen = $refSeq->length();
  $downstrStart = min($refSeqLen-1, $refEnd + $skip);
  $downstrEnd   = min($refSeqLen-1, $refEnd + $skip + $window);

  my $refDownstrSeq = Bio::Seq->new(-seq=>$refSeq->subseq($downstrStart, $downstrEnd),
				   -display_id=>"downstream");

  # search for the two regions in the query genome.

  my ($qryUpstrStart, $qryUpstrEnd,
      $qryDownstrStart, $qryDownstrEnd);

  my ($upstrHitStart, $upstrHitEnd,
      $downstrHitStart, $downstrHitEnd);

  my ($upstrHSP, $downstrHSP);
  my $ul;
  $ul = $refUpstrSeq->length();
  my $tmpio = Bio::SeqIO->new(-file => ">tmpout.fasta",
			      -format => "fasta" );
  $tmpio->write_seq($refUpstrSeq);

  $blastRes = $blastQuery->blastall($refUpstrSeq);

  # get the coordinates of the best match
  $upstrHSP = GetFirstHSP(\$blastRes);
  if ($upstrHSP == 0) {
#    print "noupstr region found \n";
    return (0,0);
  }
  ($qryUpstrStart, $qryUpstrEnd) = $upstrHSP->range('sbjct');
  ($upstrHitStart, $upstrHitEnd) = $upstrHSP->range('query');
#  print "upstream result: $qryUpstrStart, $qryUpstrEnd $upstrHitStart, $upstrHitEnd \n";
  # search against the downstr region
  $dl = $refDownstrSeq->length();
#  print "blasting len: $dl\n";
  $downstrBlastRes = $blastQuery->blastall($refDownstrSeq);
  $downstrHSP = GetFirstHSP(\$downstrBlastRes);
  if ($downstrHSP == 0) {
#    print "no downstr region found \n";
    return (0,0);
  }
  ($qryDownstrStart, $qryDownstrEnd) = $downstrHSP->range('sbjct');
  ($downstrHitStart, $downstrHitEnd) = $downstrHSP->range('query');
#  print "downstream result $qryDownstrStart $qryDownstrEnd $downstrHitStart, $downstrHitEnd \n";
  # check the span of the blast hit coordinates
  # it's possible that the blast coordinates
  my ($upstrUnaligned, $downstrUnaligned);

  $upstrUnaligned = $refUpstrSeq->length() - $upstrHitEnd;
  $downstrUnaligned = $downstrHitStart;

  my ($upstrEndAlign, $downstrStartAlign);
  $upstrEndAlign = $qryUpstrEnd + $upstrUnaligned;
  $downstrStartAlign = $qryDownstrStart - $downstrUnaligned;
#  print "ueh from $qryUpstrEnd $upstrUnaligned, dsa: $qryDownstrStart $downstrUnaligned \n";
  my ($refDiff, $qryDiff);
  $refDiff = $refEnd - $refStart;
  $qryDiff = $downstrStartAlign - $upstrEndAlign;
#  print "resulting deletioncheck: $refEnd $refStart: $refDiff $downstrStartAlign $upstrEndAlign $qryDiff\n";
  return ($upstrEndAlign, $downstrStartAlign);
}

sub GetFirstHSP {
  my ($blastReport) = @_;
  my ($result, $hit, $hsp);
  if ($result = $$blastReport->next_result) {
    my $nh;
    $nh = $result->num_hits();
    if ($hit = $result->next_hit) {
      if ($hsp = $hit->next_hsp) {
	return $hsp;
      }
    }
  }
  return 0;
}

sub GetTotalScore {
  my ($blastReport, $eThresh, $orientation) = @_;
  $totalScore = 0;
  $numHsps = 0;
  while ( $result = $$blastReport->next_result ) {
    while ($hit = $result->next_hit ) {
      while ($hsp = $hit->next_hsp ) {
	($rBegin, $rEnd) = $hsp->range('query');
	($qBegin, $qEnd) = $hsp->range('sbjct');
	$score = $hsp->score();
	$evalue = $hsp->evalue();
	if ($evalue < $eThresh) {
#	  print "    got range: $rBegin, $rEnd, $qBegin, $qEnd $score $evalue\n";
	  $totalScore += $score;
	}
	$numHsps++;
      }
    }
  }
#  print "got $numHsps hsps\n";
  return $totalScore;
}

sub SameOrder {
  my ($val1, $val2, $fact) = @_;
  if ($val1 > $val2) { 
    $p = $fact * $val2;
#  print "so $val1, $val2, $p, $fact\n";
    if ( $p > $val1) {
      return 1;
    }
    else {
      return 0;
    }
  }
  else {
    $p = $fact * $val1;
#  print "so $val1, $val2, $p, $fact\n";
    if ( $p > $val2) {
      return 1;
    }
    else {
      return 0;
    }
  }
}

sub ParseBin {
  my ($bin) = @_;
  %species = ();
  @nothing = ();
  foreach $name (keys %{$bin}) {
    ($ref, $qry, $start, $end) = ParseName($name, \@nothing);
    $species{$ref} = 1;
    $species{$qry} = 1;
  }
  @k = keys %species;
  return @k;
}

sub PrintBin {
  my ($bin, $number) = @_;
  print "bin: $number ";
  foreach $name (keys %{$bin}) {
    print " $name";
  }
  print "\n";
}

sub GetOrthoPosOnly {
  my ($from, $to, $region, $pos, $strand, $database, $skipRevLookup);
  $skipRevLookup = 0;
  my ($from, $to, $region, $pos, $strand, $database, $skipRevLookup) = @_;
  my $defaultsFile = "";
  if (exists $ENV{"GCFG"} ) {
    $defaultsFile =" -f " . $ENV{"GCFG"};
  }

  if ($database eq "") {
	$database = "NONE";
  }
  if ($skipRevLookup) {
    $nopt =" -n ";
  }
  else { 
    $nopt = "";
  }
  $cmd = $ENV{"EUSRC"} . "/comparative/". $ENV{"MACHTYPE"} . "/dbloop $defaultsFile $nopt -d $database -l $from $to $region $pos $strand";
##  print "running $cmd\n";
  $res = `$cmd`;
  if ($res =~ /error/) {
    print "error looking up coordinates: $res\n";
    return -1;
  }
  @vals = split(/\s+/, $res);
#  print "got vals: @vals\n";
  return @vals[2];
}

sub GetOrthoPos {
  my ($from, $to, $region, $pos, $strand, $database, $skipRevLookup);
  $skipRevLookup = 0;
  my ($from, $to, $region, $pos, $strand, $database, $skipRevLookup) = @_;
  my $defaultsFile = "";
  if (exists $ENV{"GCFG"} ) {
    $defaultsFile = " -f " . $ENV{"GCFG"};
  }
  if ($database eq "") {
	$database = "NONE";
  }
  if ($skipRevLookup) {
    $nopt =" -n ";
  }
  else { 
    $nopt = "";
  }
##  print "runnign gop on $from,$to,$region,$pos,$strand, $database\n";
  my $defaultsFile = "";
  if (exists $ENV{"GCFG"}) {
    $defaultsFile = "-f " . $ENV{"GCFG"};
  }
  if ($database eq "") {
    $database = "";
  }
  $cmd = $ENV{"EUSRC"} . "/comparative/". $ENV{"MACHTYPE"} . "/dbloop $nopt $defaultsFile -d $database -l $from $to $region $pos $strand";
#  print "gop running $cmd\n";
  $res = `$cmd`;
  @vals = split(/\s+/, $res);
  foreach $val (@vals) {
    print "val: $val\n";
  }
  if (@vals[0] =~ "Could") {
    print "bad connection\n";
    exit(0);
  }
  return @vals;
}

sub GetLength {
  my ($dbh, $name, $region) = @_;
#  print "running getlength\n";
  my $queryStr = "select length from sequences where name='$name' and sequence='$region'";
  my $sth = $dbh->prepare($queryStr);
  $sth->execute();
#  print "done executing\n";
  @vals = $sth->fetchrow_array;
  return @vals[0];
}


sub GetOtherStrandPosition {
  my ($dbh, $position, $seq, $region) = @_;
  my $length = GetLength($dbh, $seq, $region);
  return $length - $position;
}

sub IsNameQuery {
  my ($name) = @_;
  if ($name =~ /([\w\_]+)_(\d+)_(\d+)_([\w\_]+)/) {
    return 0;
  }
  elsif ($name =~ /([\w\_]+)_([\w\_]+)_(\d+)_(\d+)/) {
    return 1;
  }
}

sub IsNamePair {
  my ($name) = @_;
  if ($name =~ /(\w+)_(\d+)_(\d+)_(\w+)/) {
    return 1;
  }
  elsif ($name =~ /(\w+)_(\w+)_(\d+)_(\d+)/) {
    return 1;
  }
  else {
    return 0;
  }
}

sub PairsToInvList {
  my ($refSpecies, $pairs, $invList) = @_;
  print "checking $#{$refSpecies} \n";
  my $n = 0;
  for ($n=0; $n <= $#{$refSpecies}; $n++) {
    print "$n: $$refSpecies[$n]";
  }
  print "\n";
  foreach my $refKey (@{$refSpecies}) {
    my @k = keys %{$$pairs{$refKey}};
    print  " $refKey this contains $#k elements \n";
    foreach my $qryKey (keys %{$$pairs{$refKey}}) {
      if ($refKey ne $qryKey) {
	print   "     and $refKey, $qryKey contains $#{$$pairs{$refKey}{$qryKey}} inversions \n";
	for ($pairIndex = 0; 
	     $pairIndex <=  $#{$$pairs{$refKey}{$qryKey}}; 
	     $pairIndex++) {
	  $refName = $$pairs{$refKey}{$qryKey}[$pairIndex][0] . "_" .
	    $$pairs{$refKey}{$qryKey}[$pairIndex][1] . "_" .
	    $$pairs{$refKey}{$qryKey}[$pairIndex][2] . "_" .
	    $$pairs{$refKey}{$qryKey}[$pairIndex][3];
	  $qryName = $$pairs{$refKey}{$qryKey}[$pairIndex][0] . "_" .
	    $$pairs{$refKey}{$qryKey}[$pairIndex][3] . "_" .
	    $$pairs{$refKey}{$qryKey}[$pairIndex][4] . "_" .
	    $$pairs{$refKey}{$qryKey}[$pairIndex][5];
	  push @{$invList}, [ @{$$pairs{$refKey}{$qryKey}}[$pairIndex],
			      $refName,$qryName,0 ];
#	  push @{$invList}, [ @{$$pairs{$refKey}{$qryKey}}[$pairIndex],
#			      $refName,$qryName,0 ];
	}
      }
    }
  }
}

sub ReadBinFile {
  my ($binFileName, $bins) = @_;
  @{$bins} = ();
  my $index = -1;
  my $line;
  open (BINFILE, $binFileName) or die "cannot open $binFileName\n";
  while (<BINFILE>) {
    $line = $_;
    if ($line =~ /bin/) {
      push @{$bins}, ();
      ++$index;
    }
    elsif ($line =~ /\S/) {
      $line =~ /(\S+)\s+(\d+)\s+(\d+).*/;
      push @{$$bins[$index]}, [$1, $2, $3];
    }
  }
}

sub ReadFiles {
  my ($refSpecies, $pairs, $invList, $nameToSequence, $nameToIndex, $dir) = @_;
#  print "rf got @{$refSpecies} \n";
#  print "rf checking in $dir";
  foreach my $refKey (@{$refSpecies}) {
    foreach my $qryKey (keys %{$$pairs{$refKey}}) {
      if ($refKey ne $qryKey) {
	my $refSeqFile = "$dir/$refKey.$qryKey.ref";
	my $qrySeqFile = "$dir/$refKey.$qryKey.qry";
#	print "reading $refSeqFile $qrySeqFile\n";
	my $refSeqio = Bio::SeqIO->new(-file => $refSeqFile,
				       -format => "fasta" );
	my $qrySeqio = Bio::SeqIO->new(-file => $qrySeqFile,
				       -format => "fasta" );
	my $pairIndex = 0;
	while ($refSeq = $refSeqio->next_seq and
	       $qrySeq = $qrySeqio->next_seq) {
	  my $refTitle = $refSeq->display_id();
	  my $qryTitle = $qrySeq->display_id();
	  # store complete information about each inversion in invlist
	  push @{$invList}, [ @{$$pairs{$refKey}{$qryKey}}[$pairIndex],
			      $refTitle,
			      $qryTitle,
			      $refSeq ];
	  $$nameToSequence{$refTitle} = $refSeq;
	  $$nameToIndex{$refTitle} = $#{$invList};
	  push @{$invList}, [ @{$$pairs{$refKey}{$qryKey}}[$pairIndex],
			      $refTitle,
			      $qryTitle,
			      $qrySeq ];
	  # have a convenient map from name to sequence
	  $$nameToSequence{$qryTitle} = $qrySeq;
	  $$nameToIndex{$qryTitle} = $#{$invList};
	  ++$pairIndex;
	}
	if ($pairIndex != ($#{$$pairs{$refKey}{$qryKey}} + 1)) {
	  print "Number of sequences in file $dir/$refKey.$qryKey.ref does not match " .
	    "the number of inversions \n";
	  print "@{$$pairs{$refKey}{$qryKey}} $pairIndex\n";
	  exit(0);
	}
      }
    }
  }
}


sub ReadInvFile {
  my ($invFileName, $refSpecies, $pairs) = @_;
  print "reading $invFileName \n";
  open(INVFILE, $invFileName) or die "cannot open $invFileName\n";

  # parse the output of dbinvcheck
  @$refSpecies = ();
  %{$pairs} = ();
  my ($refSpec, $qrySpec, $refSeq, $qrySeq, $refFile, $qryFile);
  my %specHash = ();
  while (<INVFILE>) {
    my $line = $_;
    # look to see if this line is a header
    if ($line =~ /\// or $line =~ /\.lav/ or $line =~ /\A[^\d].*/) {
      $refSpec = ""; $qrySpec = "";
      $refSeq  = ""; $qrySeq =  "";
      if ($line =~ /\//) {
	if ($line =~ /\.*\/((\w+)\.(\w+)\.(\w+))\.((\w+)\.(\w+)\.(\w+))\..*/) {
          $refSpec = $2; $qrySpec = $6;
	  $refSeq = $3; $qrySeq = $7;
	  $refFile = $1; $qryFile = $5;
        }
        elsif ($line =~/\.*((\w+)\.(\w+))\.((\w+)\.(\w+))\..*/) {
	  $refSpec = $2; $qrySpec = $5;
          $refFile = $1; $qryFile = $4;
        }

      }
      elsif ($line =~/([^\d]\w+)\s+(\w+)\s+(\w+)/) {
	$refSpec = $1; $qrySpec = $2; $refSeq = $3; $qrySeq = $refSeq;
      }
      else {
	$line =~ /((\w+)\.(\w+)\.(\w+))\.((\w+)\.(\w+)\.(\w+))\..*/;
	$refSpec = $2; $qrySpec = $6;
	$refSeq = $3; $qrySeq = $7;
	$refFile = $1; $qryFile = $5;
      }
      if ($#{$refSpecies} < 0 or @{$refSpecies}[-1] ne $refSpec) {
	if (! exists $refHash{$refSpec}) {
	  $refHash{$refSpec} = 1;
	  push @{$refSpecies}, $refSpec;
	}
      }
    }
    else {
      $line =~/(\d+) (\d+) (\d+) (\d+).*/;
      $refStart = $1; $refEnd = $2;
      $qryStart = $3; $qryEnd = $4;
      if (!exists $$pairs{$refSpec}) {
	$$pairs{$refSpec} = ();
      }
      if (!exists ${$pairs}{$refSpec}{$qrySpec}) {
	@{${$pairs}{$refSpec}{$qrySpec}} = ();
      }
#      print "read $refSpec $refStart $refEnd $qrySpec $qryStart $qryEnd\n";
      push @{${$pairs}{$refSpec}{$qrySpec}}, [ $refSpec,
					       $refStart,
					       $refEnd,
					       $qrySpec,
					       $qryStart,
					       $qryEnd ];
      my $sz;
      $sz = $#{@{${$pairs}{$refSpec}{$qrySpec}}} + 1;
    }
  }
  close INVFILE;
}

sub ReadDataFile {
  my ($dataFileName, $data) = @_;
  open (DATA, "$dataFileName") or die "cannot open $dataFileName\n";
  while (<DATA>) {
    $line = $_;
    @values = split(/\s+/, $line);
    $key = shift @values;
    @{$$data{$key}}= (@values);
  }
}
    

sub CreateOrthoCoordinates {
  my ($dbh, $orthoPos, $invList, $regionName, $speciesList, $dbName) = @_;
  %{$orthoPos} = ();
  my $invIndex = 0;
  while ( $invIndex < $#{$invList}) {
    my ($refName, $qryName);
    $refName = $$invList[$invIndex][-3];
    $qryName = $$invList[$invIndex][-2];
#    my ($rSpec, $qSpec, $startPos, $endPos) = common::ParseName($refName, $speciesList);
    my ($rSpec, $qSpec, $startPos, $endPos) = ($$invList[$invIndex][0][0],
					       $$invList[$invIndex][0][3],
					       $$invList[$invIndex][0][1],
					       $$invList[$invIndex][0][2]);

    my ($fromSeq, $toSeq, $humanStartPos, $humanEndPos, $toStrand);
    ($fromSeq, $toSeq, $humanStartPos, $toStartStrand) =
      common::GetOrthoPos($rSpec, "human", $regionName, $startPos, 0, $dbName);
    
    ($fromSeq, $toSeq, $humanEndPos, $toEndStrand) =
      common::GetOrthoPos($rSpec, "human", $regionName, $endPos, 0, $dbName);

    if ($toStartStrand == 1 || toEndStrand == 1) {
      $humanLength = common::GetLength($dbh, "human", $regionName);
      if ($toStartStrand == 1 && toEndStrand == 1) {
	$temp = $humanEndPos;
	$humanEndPos = $humanLength - $humanStartPos;
	$humanStartPos = $humanLength - $humanEndPos;
      }
      elsif ($toStartStrand == 1) {
	my $length = $endPos - $startPos;
	$humanStartPos = $humanLength - $humanStartPos - $length;
      }
      elsif ($toEndStrand == 1) {
	my $length = $endPos - $startPos;
	$humanEndPos = $humanLength - $humanEndPos + $length;
      }
    }
    # it's possible that one of the coordinates mapped to a crazy
    # location, likely in the reverse strand code (since that's not 
    # bullet proof right now).  The fix is is the distance between
    # the orthologous points is way off, repeat the search using
    # only the primary net.

    $qryDiff = abs($endPos - $startPos);
    $humanDiff = abs($humanEndPos - $humanStartPos);
    if ($humanDiff > 3 * $qryDiff) {
      ($fromSeq, $toSeq, $humanStartPos, $toStartStrand) =
	common::GetOrthoPos($rSpec, "human", $regionName, $startPos, 0, $dbName, 1);
      
      ($fromSeq, $toSeq, $humanEndPos, $toEndStrand) =
	common::GetOrthoPos($rSpec, "human", $regionName, $endPos, 0, $dbName, 1);
    }

    if ($humanEndPos eq "common") {
      print "bad coordinate \n";
      exit(0);
    }
    @{${$orthoPos}{$refName}} = ($humanStartPos, $humanEndPos);
    # parse the query sequence
#    ($rSpec, $qSpec, $startPos, $endPos) = common::ParseName($qryName, $speciesList);

    my ($rSpec, $qSpec, $startPos, $endPos) = ($$invList[$invIndex][0][0],
					       $$invList[$invIndex][0][3],
					       $$invList[$invIndex][0][4],
					       $$invList[$invIndex][0][5]);

    my $qryLength = common::GetLength($dbh, $qSpec, $regionName);
    my ($forwardStart, $forwardEnd);
    $forwardStart = $qryLength - $endPos;
    $forwardEnd   = $qryLength - $startPos;
    ($fromSeq, $toSeq, $humanStartPos, $toStrand) =
      common::GetOrthoPos($qSpec, "human", $regionName, $forwardStart, 0, $dbName);
    ($fromSeq, $toSeq, $humanEndPos, $toStrand) =
      common::GetOrthoPos($qSpec, "human", $regionName, $forwardEnd, 0, $dbName);

    $qryDiff = abs($forwardStart - $forwardEnd);
    $humanDiff = abs($humanEndPos - $humanStartPos);
    if ($humanDiff > 3 * $qryDiff) {
      ($fromSeq, $toSeq, $humanStartPos, $toStartStrand) =
	common::GetOrthoPos($qSpec, "human", $regionName, $forwardStart, 0, $dbName, 1);
      
      ($fromSeq, $toSeq, $humanEndPos, $toEndStrand) =
	common::GetOrthoPos($qSpec, "human", $regionName, $forwardEnd, 0, $dbName, 1);
    }
    
    @{${$orthoPos}{$qryName}} = ($humanStartPos, $humanEndPos);
    $invIndex++;
  }
}

sub ReadOrthoCoordinates {
  my ($orthoFileName, $orthoPos) = @_;

  open(ORTHO, "$orthoFileName") or die "cannot open $orthoFileName\n";
  while (<ORTHO>) {
    my $line = $_;
    @vals = split(/\s+/, $line);
    @{$$orthoPos{@vals[0]}} = (@vals[1], @vals[2]);
  }
  close ORTHO;
}

sub ReadCharFile {
  my ($charFileName, $binList, $posList, $species) = @_;
  open (CHARFILE, $charFileName) or die "cannot open $charFileName \n";
  my $firstLine = <CHARFILE>;
  my $line = "";
  @$species = ();
  if ($firstLine !~ /bin/) {
    @$species = split(/\s+/, $firstLine);
  }
  else {
    $line = $firstLine;
  }

  $curBin = -1;
  while (($line ne "") or (($line = <CHARFILE>) ne "")) {
    $humanStart = 0;
    $humanEnd   = 0;
    if ($line =~ /query:/) {
      $line =~ s/query\:/0\t/g;
    }
    if ($line =~ /\Abin:\s+(\-?\d+)\s+(\-?\d+)\s+to\s+(\d+)/) {
      # reading a title line
      $binNumber = $1;
      $humanStart = $2;
      $humanEnd   = $3;
      push @$posList, [ ($humanStart, $humanEnd, $binNumber) ];
      ++$curBin;
      ${$binList}[$curBin], [()];
    }
    elsif ($line =~ /\Abin:\s+(-?\d+)\s+to\s+/) {
      $binNumber = $1;
      $humanStart = 0;
      $humanEnd = 0;
      push @$posList, [ ($humanStart, $humanEnd, $binNumber) ];
      ++$curBin;
      ${$binList}[$curBin], [()];
    }
    else {
      $line =~ /\s*(\w+)\s+(\-?\d+)\s+(\-?\d+)\s+(\d+)\s+(.*)/;
      my $spec= $1;
      my $start = $2;
      my $end = $3;
      my $nchar = $4;
      my $chars = $5;
      my @charArray = split(/\s+/, $chars);
      push @{${$binList}[$curBin]}, [($spec, $start, $end, $nchar, [@charArray])];
    }
    $line = "";
  }
}


sub PrintCharFile {
  my ($outCharFile, $species, $bins, $boundaries) = @_;

  open (CHAROUT, ">$outCharFile") or die "cannot open $outCharFile\n";

  if ($#{$species} >= 0) {
    print CHAROUT "@{$species}\n";
  }

  for my $binNum (0 .. $#$bins) {
    print CHAROUT "bin: $$boundaries[$binNum][2] $$boundaries[$binNum][0] to $$boundaries[$binNum][1]\n";
    for my $binIdx (0 .. $#{$$bins[$binNum]}) {
      print CHAROUT "\t ";
      for my $in (0 .. $#{$$bins[$binNum][$binIdx]} - 2) {
	print CHAROUT " $$bins[$binNum][$binIdx][$in]\t";
      }
      $len = $#{$$bins[$binNum][$binIdx][4]} + 1;
      print CHAROUT " $len ";
      print CHAROUT " @{$$bins[$binNum][$binIdx][4]}\n";
    }
  }
  close CHAROUT;
}



sub ReadSpeciesFile {
  my ($speciesFileName, $species) = @_;
  open(SPECIES, $speciesFileName) or die "cannot open $speciesFileName\n";

  $line = <SPECIES>;
  @$species = split(/\s+/, $line);
}
  return 1;

sub GetBinFromFileName {
  my ($name) = @_;
  $name =~ /\.(\d+)\d/;
  $value = $1;
}

sub BinFilesToHash {
  my ($pattern, $hash) = @_;
  @files = glob($pattern);
  foreach $file (@files) {
    $number = GetBinFromFileName($file);
    $$hash{$number} = $file;
  }
}
    

sub GetSequenceLengths {
  my ($pattern, $sequenceLengths) = @_;
  @lengthFiles = glob($pattern);
  foreach $lf (@lengthFiles) {
    open (LF, $lf) or die "cannot open $lf\n";
    $line = <LF>;
    $line =~ /(\S+)\s+(\d+)/;
    $$sequenceLengths{$1} = $2;
  }
}

sub ReadOrthData {
  my ($orthDataName, $orthData) = @_;
  if (! -e "$orthDataName") {
    print "$orthDataName does not exist\n";
    exit(0);
  }
  open (ORTHIN, "$orthDataName") or die "cannot open $orthDataName\n";
  while (<ORTHIN>) {
    my $line = $_;
    $line =~ /(\w+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/;
    my $spec = $1; my $s = $2; my $e = $3;
    my $hs = $4; my $he = $5;
    $$orthData{$spec} = [ $s, $e, $hs, $he ];
  }
}

sub OrderedKeys {
  my ($hash) = @_;
  my @unsKeys = keys %$hash;
  if ($#unsKeys == 0) {
    return @unsKeys;
  }
  else {
    if ($unsKeys[0] == 0 and $unsKeys[0] ne "0") {
      @sKeys = sort @unsKeys;
    }
    else {
      @sKeys = sort {$a <=> $b} @unsKeys;
    }
  return @sKeys;
  }
}

sub ComputeOverlap {	
  my ($as, $ae, $bs, $be) = @_;
  $overlap = 0;
  if ($ae < $be) {
    if ($ae < $bs) {
      $overlap = 0;
    }
    else {
      $overlap = $ae - max($as, $bs);
    }
  }
  else {
    if ($as > $be ) { 
      $overlap = 0;
    }
    else {
      $overlap = $be - max($as, $bs);
    }
  }
  return $overlap;
}
	

sub ReadCoordsFile {
  my ($coordsFileName, $coords, $bins) = @_;
  open (CF, $coordsFileName) or die "cannot open $coordsFileName\n";

  my $bin = -1;
  while (<CF>) {
    my $line = $_;
    if ($line =~ /bin: (\d+)/) {
      $binNum =$1;
      $bin++;
      @{$coords}[$bin] = [];
      push @{$bins}, $binNum ;
    }
    else {
      if ($bin == -1) {
	print "coords file must start with a bin\n";
	exit(0);
      }
      # trim leading whitespace
      $line =~ /\s*(\S.*)/;
      $line = $1;
      my @vals  = split(/\s+/, $line);
      push @{@{$coords}[$bin]}, [@vals];
    }
  }
}
