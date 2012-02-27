package grimm;


sub ReadSyntenyBlocks {
  my ($syntBlockFile, $blocks) = @_;
# parse the format
# genome1: genome1
# genome2: genome2
# block_id genome1_chr genome1_start genome1_len genome1_sign genome2_chr genome2_start genome2_len genome2_sign
#1 1 1 204542 + 1 2070 256962 +
  open(SB, $syntBlockFile) or die "cannot open $syntBlockFile\n";
  <SB>; <SB>;
  $b = 0;
  my @block;
  while(<SB>) {
    @block = split(/\s+/, $_);
    push @$blocks,  {};
    $$blocks[$b]{"blockid"} = $block[0];
    $$blocks[$b]{"chr1"} = $block[1];
    $$blocks[$b]{"start1"} = $block[2];
    $$blocks[$b]{"len1"} = $block[3];
    $$blocks[$b]{"sign1"} = $block[4];
    $$blocks[$b]{"chr2"} = $block[5];
    $$blocks[$b]{"start2"} = $block[6];
    $$blocks[$b]{"len2"} = $block[7];
    $$blocks[$b]{"sign2"} = $block[8];
    $b++;
  }
}

sub ParseMicroGrimmBlocks {
  my ($blockFile, $blocks, $order) = @_;
  $curBlock = -1;
  while ($#$blockFile >= 0) {
    $curBlock++;
    push @$blocks, ();
    unless ($$blockFile[0] =~
	    /begin_block (\d+) (\d) (\d+) (\d+) ([+\-]) (\d+) (\d+) (\d+) ([+\-])/) {
      print "expecting mgr_micro header \n";
      print "example:\n";
      print " # begin_block 1 1 3 1888127 + 1 3 1877420 +\n";
      print "got: \n";
      print $$blockFile[0];
      return 1;
    }
    shift @$blockFile;
    unless ($$blockFile[0] =~ /begin_anchors/) {
      print "expecting mgr_micro delineator\n";
      print "example: \n";
      print "# begin_anchors\n";
      print "got:\n";
      print $$blockFile[0];
      return 1;
    }
    shift @$blockFile;
    $curStrip = 0;
    # read in the anchors
    while ($#$blockFile >= 0 and
	   $$blockFile[0] !~ /end_anchors/) {
      unless ($$blockFile[0] =~ /(\d+) (\d) (\d+) (\d+) ([+\-]) (\d+) (\d+) (\d+) ([+\-])/) {
	print "malformatted anchor line: \n";
	print $$blockFile[0];
	return 1;
      }
      $strip = $1;
      if ($strip > $curStrip) {
	$curStrip = $strip;
	push @{$$blocks[$curBlock]}, ();
	$anchor = 0;
      }
      push @{$$blocks[$curBlock][$curStrip-1]}, {};
      @block = split(/\s+/, $$blockFile[0]);
      $$blocks[$curBlock][$curStrip-1][$anchor]{"blockid"} = $block[0];
      $$blocks[$curBlock][$curStrip-1][$anchor]{"chr1"} = $block[1];
      $$blocks[$curBlock][$curStrip-1][$anchor]{"start1"} = $block[2];
      $$blocks[$curBlock][$curStrip-1][$anchor]{"len1"} = $block[3];
      $$blocks[$curBlock][$curStrip-1][$anchor]{"sign1"} = $block[4];
      $$blocks[$curBlock][$curStrip-1][$anchor]{"chr2"} = $block[5];
      $$blocks[$curBlock][$curStrip-1][$anchor]{"start2"} = $block[6];
      $$blocks[$curBlock][$curStrip-1][$anchor]{"len2"} = $block[7];
      $$blocks[$curBlock][$curStrip-1][$anchor]{"sign2"} = $block[8];
      ++$anchor;
      shift @$blockFile;
    }
    shift @$blockFile;
    # parse the grimm
    unless ($$blockFile[0] =~ /begin_mgr/) {
      print "expected begin_mgr\n";
      print "got: \n";
      print $$blockFile[0];
      return 1;
    }
    shift @$blockFile;
    push @$order, ();

    for ($i = 1; $i <= 2; $i++ ) {
      unless ($$blockFile[0] =~ />genome$i/) {
	print "expected >genome$i\n";
	print "got: \n";
	print $$blockFile[0];
	return 1;
      }
      shift @$blockFile;
      # read the order
      @orderl = split(/\s+/, $$blockFile[0]);
      $orderli = join(",", @orderl);
      push @{$$order[$curBlock]}, [@orderl];
      shift @$blockFile;
    }
    unless ($$blockFile[0] =~ /end_mgr/) {
      print "expected end_mgr\n";
      print "got: \n";
      print $$blockFile[0];
      return 1;
    }
    shift @$blockFile;
    unless ($$blockFile[0] =~ /end_block/) {
      print "expected end_block\n";
      print "got: \n";
      print $$blockFile[0];
      return 1;
    }
    shift @$blockFile;
  }
  return 0;
}

sub FindInversions {
  my($ord1, $ord2, $inv) = @_;

  # make sure ord1 is the identity 
  my $i;
  for ($i = 0; $i < $#$ord1; $i++ ) {
    if ($$ord1[$i] != $i+1) {
      print "order1 should be the identity\n";
      return 1;
    }
  }
  $os = scalar @{$ord2};
  for ($i = 0; $i <= $#{$ord2}; $i++ ) {
    if (-1*$$ord2[$i] == ($i + 1)) {
      push @$inv, $i;
    }
  }
}

sub GetInvertedAnchors {
  my ($blocks, $inv, $anchors) = @_;
  $ninv = scalar @$inv;
  for ( $i = 0; $i <= $#{$inv}; $i++ ) {
    push @$anchors, {};
    $block = $$inv[$i];
    $minStart1 = $$blocks[$block][0]{"start1"};
    $maxEnd1   = $$blocks[$block][0]{"len1"} + $minStart1;
    $minStart2 = $$blocks[$block][0]{"start2"};
    $maxEnd2   = $$blocks[$block][0]{"len2"} + $minStart2;
    for ($b = 0; $b <= $#{$$blocks[$block]}; $b++) {
      $s1 = $$blocks[$block][$b]{"start1"};
      if ($$blocks[$block][$b]{"start1"} < $minStart1) {
	$minStart1 = $$blocks[$block][$b]{"start1"};
      }
      if ($$blocks[$block][$b]{"start1"} +
	  $$blocks[$block][$b]{"len1"} > $maxEnd1) {
	$maxEnd1 = $$blocks[$block][$b]{"start1"} +
	  $$blocks[$block][$b]{"len1"};
      }
      if ($$blocks[$block][$b]{"start2"} < $minStart2) {
	$minStart2 = $$blocks[$block][$b]{"start2"};
      }
      if ($$blocks[$block][$b]{"start2"} +
	  $$blocks[$block][$b]{"len2"} > $maxEnd2) {
	$maxEnd2 = $$blocks[$block][$b]{"start2"} +
	  $$blocks[$block][$b]{"len2"};
      }
    }
    $$anchors[$i]{"start1"} = $minStart1;
    $$anchors[$i]{"len1"} = $maxEnd1 - $minStart1;
    $$anchors[$i]{"start2"} = $minStart2;
    $$anchors[$i]{"len2"} = $maxEnd2 - $minStart2;
  }
}

return 1;
