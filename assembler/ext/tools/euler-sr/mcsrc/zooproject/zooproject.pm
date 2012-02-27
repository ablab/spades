use Bio::DB::GenBank;

package zooproject;

sub ParseNameTable {
  my ($tableFileName, $niscToTrunc, $niscToId, $nameToId, $idToName) = @_;
  open(NAMETAB, $tableFileName) or die "cannot open $tableFileName\n";
  while (<NAMETAB>) {
    $line = $_;
    chomp $line;
    if ($line !~ /^#/) {
      @values = split(/\t+/, $line);
      # map from NISC
      ${$niscToTrunc}{$values[0]} = $values[1];
      ${$niscToId}{$values[0]} = $values[2];
      @{$idToName}[$values[2]] = $values[1];
      @{$nameToId}{$values[1]} = $values[2];
    }
  }
}


sub ParseCloneLine {
  my ($cloneLine, $cloneValues) = @_;
  my @values;
  @values = split(/\t+/, $cloneLine);
  my $jv = join(",", @values);
  ${$cloneValues}{"NISC Code"} = $values[0];
  ${$cloneValues}{"Organism"} = $values[1];
  ${$cloneValues}{"CloneName"} = $values[2];
  ${$cloneValues}{"GenbankId"} = $values[3];
  ${$cloneValues}{"Status"} = $values[4];
  if ($#values != 4 or $values[0] eq ""
      or $values[1] eq ""
      or $values[2] eq ""
      or $values[3] eq ""
      or $values[4] eq "") {
    return 0;
  }
  else {
    return 1;
  }
}

sub ParseTargetTable {
  my ($targetTableName, $cloneInfo) = @_;
  open(TARGETTAB, $targetTableName) or die "cannot open $targetTableName\n";
  # skip the header
  <TARGETTAB>;
  while(<TARGETTAB>) {
    $line = $_;
    chomp $line;
    my %clone;
    if (ParseCloneLine($line, \%clone) == 1) {
      push @$cloneInfo, {%clone};
#      print "added clone " . $clone{"Organism"} . " " . $clone{"Genbank Id"} . "\n";
    }
#    else {
#      print "didn't parse a clone line $line\n";
#    }
  }
}

sub FetchSequences {
  my ($target, $species, $accids) = @_;
  $db_obj = new Bio::DB::GenBank;
  foreach $a (0 .. $#$accids) {
    print "fetching $target $species $$accids[$a]: ";
    $seq_obj = $db_obj->get_Seq_by_acc($$accids[$a]);
    $seqFileOutName = ">" . $target . "." 
      . $species .  "." . $$accids[$a] . ".gbk";
    if ($seq_obj) {
      print "ok\n";
      $numOk++;
      $seqIO = Bio::SeqIO->new(-file => $seqFileOutName,
			       -format => 'genbank');
      $seqIO->write_seq($seq_obj);
    }
    else {
      print "FAILED\n";
    }
  }
  $numAccIds = scalar @$accids;
  return ($numOk, $numAccIds);
}

return 1;
