#!/usr/bin/env perl

if ($#ARGV < 1) {
  print "usage: $0 mapfile output [-g gappenalty] [-t threshold]\n";
  print "   Joins together blast hits that should be contiguous, but were split\n";
  print "   because blast didn't like the gap, or because there were N's, etc.\n";
  exit(0);
}

$infile = shift @ARGV;
$outfile = shift @ARGV;
$gapPenalty = 2;
$gapThreshold = 0;
if ($#ARGV >= 0 ) {
  $opt = shift @ARGV;
  if ($opt eq "-g") {
    $gapPenalty = shift @ARGV;
  }
  if ($opt eq "-t") {
    $gapThreshold = shift @ARGV;
  }
}
open(MAP, "$infile") or die "cannot open $infile\n";
open(OUT, ">$outfile") or die "cannot open $mapfile\n";
$first = 1;
while (<MAP>) {
  $line = $_;
  if ($line =~ /^>(.*)/) {
    # reset the hits list
    #    print "creating new hits\n";
    $nextTitle = $1;
    $prevTitle = $title;
    $title = $nextTitle;
    if ($#hits >= 1) {
#      print "joining hits\n";
      JoinHits(\@hits);
#      print "done joining hits\n";
    }
    if ($first == 0) {
      print OUT ">$prevTitle\n";
      for ($h = 0; $h <= $#hits; $h++) {
	$line = join("\t", @{$hits[$h]});
	print OUT "$line\n";
      }
    }
    $first = 0;
    @hits = ();
  } else {
    @vals = split (/\s+/, $line);
    push @hits, [@vals];
  }
}
if ($first == 0) {
  if ($#hits >= 1) {
    JoinHits(\@hits);
  }
}


sub JoinHits {
  my ($hits) = @_;
  $mergeFound = 1;

  while ($mergeFound and
	 scalar @{$hits} > 1) {
    $nhits = scalar @{$hits};
#    print "checking $nhits hits\n";
    $minDist = 999999999999999;
    $mini = -1;
    $minj = -1;
    # find the two closest hsps
    $mergeFound = 0;
    for ($i = 0; $i <= $#{$hits}-1; $i++) {
      $j = $i + 1;
      while ($j <= $#{hits}) {
	# determine the physical order
	$merged = 0;
	$first = $i;
	$second = $j;
	if ($$hits[$j][1] < $$hits[$i][1]) {
	  $first = $j;
	  $second = $i;
	}
	# Make sure that the hsps do not overlap
	if ($$hits[$first][1] < $$hits[$second][1] and
	    $$hits[$first][2] < $$hits[$second][1]) {
	  # Make sure the query hits are in the proper orientation
	  # and do notoverlap as well
#	  print "comparing : \n";
#	  print "@{$hits[$first]}\n";
#	  print "@{$hits[$second]}\n";
	  if ($$hits[$first][6] == $$hits[$second][6] and 
	      ($$hits[$first][6] == 1 and 
	       $$hits[$first][3] < $$hits[$second][3] and
	       $$hits[$first][4] < $$hits[$second][3]) or 
	      ($$hits[$first][6] == -1 and
	       $$hits[$second][3] < $$hits[$first][3] and
	       $$hits[$second][4] < $$hits[$first][3])) {
	    # Find the score of the alignment based on matches on the subject
	    $scoreFirst  = (($$hits[$first][2] - $$hits[$first][1])*$$hits[$first][9]);
	    $scoreSecond = (($$hits[$second][2] - $$hits[$second][1])*$$hits[$second][9]);

	    # find the distance between subject hits (likely
	    # to be further apart than query, the read

	    if ($$hits[$first][6] == 1) {
	      $dist = $$hits[$second][3] - $$hits[$first][4];
	    }
	    else {
	      $dist = $$hits[$first][3] - $$hits[$second][4];
	    }
	    $gap  = $dist * $gapPenalty;
	    $score = $scoreFirst + $scoreSecond;

	    if (($gap - $score) < $gapThreshold) {
	      # create a copy of everything
	      $firstQStart = $$hits[$first][1];
	      $secondQEnd  = $$hits[$second][2];

	      if ($$hists[$first][6] == 1) {
		$firstSStart = $$hits[$first][3];
		$secondSEnd  = $$hits[$second][4];
	      }
	      else {
		$firstSStart = $$hits[$second][3];
		$secondSEnd  = $$hits[$first][4];
	      }

	      $$hits[$i][1] = $firstQStart;
	      $$hits[$i][2] = $secondQEnd;
	      $$hits[$i][3] = $firstSStart;
	      $$hits[$i][4] = $secondSEnd;
#	      print "joining: $prevTitle \n";
#	      print "@{$hits[$first]}\n";
#	      print "@{$hits[$second]}\n";
	      splice(@{$hits}, $j, 1);
	      $mergeFound = 1;
	      $nhits = scalar @$hits;
	      $merged = 1;
	    }
	  }
	}
	if ($merged == 0) {
	  $j++;
	}
      }
    }
  }
#  exit(0);
}

