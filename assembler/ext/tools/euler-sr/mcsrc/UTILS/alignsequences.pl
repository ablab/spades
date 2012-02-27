#!/usr/bin/perl 
use strict;
use POSIX;
use FileHandle;

STDOUT->autoflush;
die "usage: alignsequences correct modified [original] \n" if (@ARGV != 3);

open (CORRECT, @ARGV[0]) or die "cannot open @ARGV[0]\n";
open (MODIFIED, @ARGV[1]) or die "cannot open @ARGV[1]\n";
open (ORIGINAL, @ARGV[2]) or die "cannot open @ARGV[2]\n";


my $l1;
my $l2;

my $debug = 0;

my $continue = 1;
my ($seq1, $title1, $index1, $continue1);
my ($seq2, $title2, $index2, $continue2);
my ($seq3, $title3, $index3, $continue3);

my $k = 11;
my $maxlen = 500;

my ($i, $j);

my @scoreMat = ();
my @pathMat = ();
my @newMat = ();
my ($imax, $jmax);
foreach $i (0..($k*2+2)) {
  push @newMat, 0;
}
  
my($seqa, $seqb, $orig);
$seqb = "cttttttttttttttttttggtttttcgagacagggtttctctgtatatccctggctgtcctggaa";
$seqa = "atgcagtcttttttttttttttttttttttttttggcttttcgagacagggtttctctgtatatccctggctgtcctggaa";
$orig = "atgcagtcttttttttttttttggcttttcgagacagggtttctctgtatatccctggctgtcctggaa";


#$seqa = "tggaccacaacttgtgagggcagcctatacaaggacatggaagaaggaagattgctcttcacctgtttgctctatctctt";
#$seqb = "tggaccacaacttgtgagggcagcactatacaaggacatggaagaaggaagattgctcttcacctgtttgctctatctctt";


#$seqb = "tggtctgtatctacctagcgcgtgcacagactcgtgaatgaatgagaaggcaaactcagcagatgcaggatgctctgctctaagatgccaccataactacagaa";
#$seqa =      "tgtatctacctagcgcgtgcacagactcgtgaatgaatgagaaggcaaactcagcagatgcaggatgctctgctctaagatgccaccataactacagaa";

#$seqb = "tggtctgtatctacctagcgcgtgcaca";
#$seqa =      "tgtatctacctagcgcgtgcaca";


#$seqb = "tggtctgtatctacctagcgcgtgcaca";
#$seqa = "tggtctgtatctacctagcgcgtgcaca";

#$seqb = "aacatcgctttacccttacccggttcaagtttgacgacctttgaaatgccaaaggtccgctcgacttcttagtaacggtaacaccaga";
#$seqa = "acatcgctttacccttacccgaagtttgacgacctttgaaatgccaaaggtccgctcgacttcttagtaatcggtaacaccaga";

#$seqa =
#"tagtctacgacgctcgatggtaagaatatcgcaccccccgtagtagaagtagaataccatgacgaacggatagggtggttatgccagtacgttattgccgcagac"; 
#$seqa = "gtagtctacgacgtcgatggtatgaatatcgcatccccccgttagta";

foreach $i (0..$maxlen-1) {
 $scoreMat[$i] = [@newMat];
 $pathMat[$i] = [@newMat];
}

my ($modAlign, $origAlign);
my ($modImax, $origImax);
my ($err, $corr);
my ($bcopya, $bcopyo);
if ($debug) {
  ($i, $modImax, $jmax)  = kbandAlign($seqa, $seqb, $k, \@scoreMat, \@pathMat);
  for $i (0.. length($seqb)) {
    for $j (0.. $#{$pathMat[$i]}) {
      print "$pathMat[$i][$j] ";
    }
    print "\n";
  }
  ($modAlign, $seqa, $bcopya) = tracealign(\@pathMat, $modImax, $jmax, $seqa, $seqb);

  print "got score: $i\n";
  print "blah!\n";
  ($i, $origImax, $jmax)  = kbandAlign($orig, $seqb, $k, \@scoreMat, \@pathMat);
  ($origAlign, $orig, $bcopyo) = tracealign(\@pathMat, $origImax, $jmax, $orig, $seqb);

  for $i (0.. length($seqb)) {
    for $j (0.. $#{$pathMat[$i]}) {
      print "$pathMat[$i][$j] ";
    }
    print "\n";
  }

  print "got score: $i\n";
  print "got mod  align align: $modAlign\n";
  print "got orig align align: $origAlign\n";
  print "staring mod imax: $modImax,  orig imax: $origImax \n";
  print "seqa: '$seqa', orig: '$orig' \n";
#  ($err, $corr) = getAlignStats($seqb, $orig, $seqa, $origAlign,
#  $modAlign, $origImax, $modImax); 
  ($err, $corr) = compareAlignments($origAlign, $modAlign, $origImax, $modImax);
  print "got $err errors, $corr corruptions\n";
#  exit(0);
}

$continue1 = getRead(\*CORRECT, \$title1, \$seq1);
$continue2 = getRead(\*MODIFIED, \$title2, \$seq2);
$continue3 = getRead(\*ORIGINAL, \$title3, \$seq3);
if (not $continue1) {
  print "no reading done\n";
  exit(0);
}

my ($prev1, $prev2, $prev3);
$index1 = getIndex($title1);
$index2 = getIndex($title2);
$index3 = getIndex($title3);
$prev1 = $index1;
$prev2 = $index2;
$prev3 = $index3;

my ($score, $errors, $corruptions);
my ($modErrors, $origErrors) = (0,0);
my $netCorruptions = 0;
my $netsize = 0;
my ($modImax, $origImax);
my ($modScore, $origScore);
my ($seqCopy1, $seqCopy2);
my $i = 0;
my $skipped = 0;
while ($continue1) {
  $i = $i + 1;
#  if ($i % 100 == 0) {
#  }
  $index1 = getIndex($title1);
  $index2 = getIndex($title2);
  $index3 = getIndex($title3);
  if ($index1 % 100 == 0) {
#    printf("$index1\n");
  }
  # Advance correct, modified, and original to point at the same reads.
  $continue1 = getEqualReads(\*CORRECT, \*MODIFIED, \*ORIGINAL, \$title1, \$title2, \$title3, \$seq1, \$seq2, \$seq3);
  $netsize = $netsize + length($seq1);

  # Align the modified sequence to the corrrect one.  Count the number of errors here.
  resetMat(\@scoreMat, length($seq1), $k*2+2);
  resetMat(\@pathMat, length($seq1), $k*2+2);
  ($modScore, $modImax, $jmax) = kbandAlign($seq2, $seq1, $k, \@scoreMat, \@pathMat);
  if ($modScore < 0)  {
    $score = abs(length($seq1) - length($seq2));
    $modErrors = $modErrors + $score;
    $skipped = $skipped + 1;
    next;
  }
  else {
    $score = $modScore;
    ($modAlign, $seqCopy1, $seqCopy2) = tracealign(\@pathMat, $modImax, $jmax, $seq1, $seq2);
    if($debug) {
      for $i (0.. length($seqb)) {
	 for $j (0.. $#{$pathMat[$i]}) {
	   print "$pathMat[$i][$j] ";
	 }
	 print "\n";
      }
    }
  }
  $modErrors = $modErrors + $score;

  # Align the original sequence to the correct one.  The number of
  # errors reported here will be the total.
  resetMat(\@scoreMat, length($seq1), $k*2+2);
  resetMat(\@pathMat, length($seq1), $k*2+2);

  ($origScore, $origImax, $jmax) = kbandAlign($seq1, $seq3, $k, \@scoreMat, \@pathMat);
  $origErrors = $origErrors + $origScore;
  ($origAlign, $seqCopy1, $seqCopy2) = tracealign(\@pathMat, $origImax, $jmax, $seq1, $seq3);
    print "seq: $i $modScore $origScore\n";

  ($errors, $corruptions) = compareAlignments($origAlign, $modAlign, $origImax, $modImax);
  $netCorruptions = $netCorruptions + $corruptions;
  if ($score > 0) {
    if (0) {
    print "$index1 got score: $score\n";
    print "seq: $title1";
    print "aligning: $seq1\n";
    print "          $seq2\n";
    print "orig align: $origAlign\n";
    print "mod  align: $modAlign\n";
    print "got $errors errors and $corruptions corruptions, while score was: $score, orig: $origScore\n";
    print "\n";
  }
  }
  if ($prev1 + 1 < $index1) {
    foreach $i ($prev1 + 1 ... $index1 -1) {
#      print "skipped lhs $i\n";
    }
  }
  $prev1 = $index1;
  $prev2 = $index2;
}
print "skipped : $skipped\n";
my $pwd = `pwd`;
chomp($pwd);
print "$pwd: mod_errors: $modErrors  net_ntds: $netsize net_orig_errors: $origErrors net_corruptions: $netCorruptions\n";

sub getRead {
  my($file, $title, $seq) = @_;
  my $line = "";
  while (<$file>) {
    if ($_ !~ />/) {
      chomp;
      $line  = $line . $_;
    }
    else {
      last;
    }
  }
  $$seq  = $line;
  if ($_ =~ />/) {
    $$title = $_;
    return 1;
  }
  else { return 0; }
}

  
sub getIndex {
  my($title) = @_;
  $title =~ />.*\D(\d+)\b/;
  return $1;
}


sub getEqualReads {
  my($FILE1, $FILE2, $FILE3, $title1, $title2, $title3, $seq1, $seq2,$seq3) = @_;

  # FILE1 and FILE 3 should be incremented at the same time.
  my($index1, $index2,$index3);
  my($continue1, $continue2,$continue3);

  $index1 = getIndex($$title1);
  $index2 = getIndex($$title2);
  $index3 = getIndex($$title3);

  $continue1 = 1;
  $continue2 = 1;
  $continue3 = 1;
  while ($continue1 and $continue2 and   $continue3) {
    while ($index1 < $index2 && $continue1) {
      $continue1 = getRead($FILE1, $title1, $seq1);
      $continue3 = getRead($FILE3, $title3, $seq3);
      $index1 = getIndex($$title1);
      $index3 = getIndex($$title3);
    }
    while ($index2 < $index1 && $continue2) {
      $continue2 = getRead($FILE2, $title2, $seq2);
      $index2 = getIndex($$title2);
    }
    if ($index1 eq $index2) {
      $continue1 = getRead($FILE1, $title1, $seq1);
      $continue2 = getRead($FILE2, $title2, $seq2);
      $continue3 = getRead($FILE3, $title3, $seq3);
      return $continue1 && $continue2;
    }
  }
  return 0;
}

sub max {
  my(@args) = @_;
  my $i;
  my $max_ = 0;
  for $i (0 ... $#args) {
    if ($max_ < @args[$i]) {
      $max_ = @args[$i];
    }
  }
  return $max_;
}

sub min {
  my(@args) = @_;
  my $i;
  my $min_ = 9999999999;
  for $i (0 ... $#args) {
    if ($min_ > @args[$i]) {
      $min_ = @args[$i];
    }
  }
  return $min_;
}

sub resetMat  {
  my($mat, $rows, $columns) = @_;
  my($i, $j);
  for $i (0 .. $rows - 1) {
    for $j ( 0 .. $columns - 1) {
      ${$mat}[$i][$j] = 0;
    }
  }
}


sub kbandAlign {
  my ($seqA, $seqB, $k, $scoremat, $pathmat) = @_;
  my $temp;
#  if (length($seqA) > length($seqB)) {
#    $temp = $seqA;
#    $seqA = $seqB;
#    $seqB = $temp;
#  }
  my ($gap, $mismatch);
  $gap = 1;
  $mismatch = 1;
  my ($m);
  my ($i, $j);
  my ($lb, $ub);
  # score match, score gap up (insert) score gap left (delete), placeholder for match
  my ($sm, $sgu, $sgl, $tempsm);
  my $lenA = length($seqA) - 1;
  my $lenB = length($seqB) - 1;
  my $pathmatstr = "";
  my ($c, $offset);
  $offset = 0;
  # make sure the alignment will fit in the band.
  if ( abs ($lenA  - $lenB) > $k) {
    return (-1, 0,0);
  }
  my $p;
  for $i (0.. $lenB) {
    if ($i - $k < 0) { $lb = -1*$i; } else {$lb = -1*$k;};
    if ($i + $k > $lenA) { $ub = $lenA - $i;} else {$ub = $k}
    $p = 0;
    for $j ($lb .. $ub) {
      if (substr($seqB, $i,1) eq substr ($seqA, $i + $j,1)) {$sm = 0;} else {$sm = 1}
      $tempsm = $sm;
      # simulate diagonal of large matrix on smaller matrix.  
      $sm = $sm + ${$scoremat}[$i+1-1][$k+$j-1+1];
      $sgu = 1 + ${$scoremat}[$i+1-1][$k+$j+1];
      $sgl = 1 + ${$scoremat}[$i+1][$k+$j-1];
      if ($j == -1*$k) {
	  $m = min($sm, $sgu);
	  $sgl = -1;# prevent any spurrious paths
      }
      elsif ($j == $ub) {
#	print "getting m from $sm, $sgl\n";
	$m = min($sm, $sgl);
	$sgu = -1; # prevent any spurrious paths
      }
      else {
	# allow for nonlocal alignments, reset cost to 0 if within a certain distance
	# from the beginning.
	$m = min($sm, $sgl, $sgu);
	if ($i < 5 && $m > 0) {
#	  printf("m: %d  sm: %d  tempsm:  %d\n", $m, $sm, $tempsm);
#	  printf("setting pathmat %d %d to *\n", $i +1, $k + $j);
	  ${$pathmat}[$i+1][$k + $j] = "*";
	  next;
	}
      }
      if ($m == $sm) {
	  if ($tempsm == 0) { ${$pathmat}[$i+1][$k + $j] = "\\";}
          else { ${$pathmat}[$i+1][$k + $j] = "/";}
        }
	elsif ($m == $sgl) {
	  ${$pathmat}[$i+1][$k + $j] = "-";
        }
        elsif ($m == $sgu) {
          ${$pathmat}[$i+1][$k + $j] = "|";
        }
      else {
	print "no minimum found here $i $k + $j\n";
	${$pathmat}[$i+1][$k + $j] = "*";
      }
    $pathmatstr = $pathmatstr . " " . ${$pathmat}[$i + 1][$k + $j];

    if ($debug) {
#      printf("(%s, %s)%d ", substr($seqB, $i,1) , substr ($seqA, $i + $j,1), $m);
    }

#     printf(" : %d ", $m);
#printf("(%d %d) ", $i, $i + $j);
      $p = $p + 1;
      ${$scoremat}[$i + 1][$k + $j] = $m;
      }
     $pathmatstr = $pathmatstr . "\n";
      for $c (0.. $i*2) {
	$pathmatstr = $pathmatstr . " " ;
      }
     $offset = $offset + 1;
#if ($debug) {  printf("\n");}
  }
#if ($debug) { print "$pathmatstr \n"; }

# anchor the solution to the first spot.
${$pathmat}[0][$k] = "*";
#printf("final alignment position: %d, %d\n", $lenB, $lenA - $lenB+ $k);
#  printf("final alignment score:%d\n", ${$scoremat}[$lenB+1][$lenA - $lenB + $k]);
my $minscore = 1000;
my $index;
my ($mini, $minj);
foreach $i (0 ... $k) {
  # Search along the diagonal here, or 
  if ($lenA - $lenB + $i + $k > (2*$k)) {
    $index = $2*$k;
  }
  # bottom here.
  else {
    $index = $lenA - $lenB + $i + $k
  }
#  if ($debug){  printf("checking out score: ${$scoremat}[$lenB+1-$i][$index] %d %d\n",$lenB+1-$i, $index) ;}
  if (${$scoremat}[$lenB+1-$i][$index] < $minscore) {
    $mini = $lenB + 1 - $i;
    $minj = $index;
    $minscore = ${$scoremat}[$lenB+1-$i][$index];
  } 
}
#  print "kband align returning: $minscore, $mini, $minj\n";
  return ($minscore, $mini, $minj);
#  return ${$scoremat}[$lenB+1][$lenA - $lenB + $k];
}


sub tracealign {
  my($alignmat, $iopt, $jopt, $stra, $strb) = @_;
#  print "in tracealign!";
  my ($i,$j) = ($iopt, $jopt);
  my ($i2, $j2);
  my $optalign = "";
  my $optdisplay = "";
  my ($optstra, $optstrb) = ("", "");
  my ($copya, $copyb) = ($stra, $strb);
#print "TA: iopt: $iopt, jopt: $jopt, char: ${$alignmat}[$iopt][$jopt] b: $b \n";
  my $optpath = "";

# build the resulting alignment string by tracing back a path.
# Strings are built by prepending appropriate alignment characters 
# since tracing backwards in path.
while (${$alignmat}[$i][$j] ne "*" && $i > 1 && $j >= 0 && ord(${$alignmat}[$i][$j]) != ord(0)) { 
  $optpath = ${$alignmat}[$i][$j] . $optpath;
  if ($debug) { $optalign = ${$alignmat}[$i][$j] . $optalign;}
  if (${$alignmat}[$i][$j] eq "\\") {
    # match 
    $i = $i - 1;
    if ($debug) {
      $optdisplay = "|" . $optdisplay ; 
      $optstra = chop($stra) . $optstra ; 
      $optstrb = chop($strb) . $optstrb;
    }
  }
  elsif (${$alignmat}[$i][$j] eq "|") {
    #delete string b
    $i = $i - 1; 
    $j = $j + 1;
    if ($debug) {
#      print "putting a space in a\n";
      $optdisplay = " " . $optdisplay; 
      $optstra = " " .  $optstra;  
      $optstrb = chop($strb) . $optstrb;
    }
  }
  elsif (${$alignmat}[$i][$j] eq "-") {
    #insert into string b
    $j = $j - 1;
    if ($debug) {
      $optdisplay = " " . $optdisplay;
      $optstra = chop($stra) . $optstra;
      $optstrb = " " .  $optstrb;
    }
  }
  elsif (${$alignmat}[$i][$j] eq "/") {
    # mismatch
    $i = $i - 1; 
    if ($debug) {
      $optdisplay = " " . $optdisplay ; 
      $optstra = chop($stra) . $optstra ; 
      $optstrb = chop($strb) . $optstrb;
    }
  }
}
  # i ranges from length of the string down to 1,
  my $endj;
  if ($i < $k) {
    # Since the diagonal is shifted in the matrix representation, the
    # real position in the end alignment is shifted to the left by k,
    # then back to the right the further down.  
    # since $i starts at 1, subtract 1 for correct index.
    $endj = $k - $j + $i - 1;
  }
  else {
    $endj = $j;
  }

#  print "endj : $endj\n";

#  print "done tracing path, endi: $i, endj: $j \n";
  if ($endj > $i && $i > 0) {
    for $i2 ($i .. $endj - 1) {
      $optpath = "-" . $optpath;
      $copya = " " . $copya;
    }
  }
  elsif ($i > $endj && $endj > 0) {
    for $j2 ($j .. $i) {
      $optpath = "|" . $optpath;
      $copyb = " " . $copyb;
    }
  }
  return ($optpath, $copya, $copyb);

#if ($debug) { print "o:  $optalign \n"; print "a: $optstra\n"; print "d: $optdisplay\n"; print "b: $optstrb\n";}
## prepend the first character (at the optimal position).
#$optstra = chop($stra) . $optstra;
#$optstrb = chop($strb) . $optstrb;
#$optdisplay = " " . $optdisplay;
## Terminating $i early skips a portion of b.  Shift output to right.
#print "fixing a: 1 .. $i\n";
#for $i2 (1.. $i-1) {
#  $optstra = " " . $optstra;
#  $optstrb = chop($strb) . $optstrb;
#  $optdisplay = " " . $optdisplay;
#}
## Terminating $j late skips a portion of a.
#print "fixing b: $k .. $j\n";
#for $j2 ($k.. $j-1) {
#  print "adding to optstring b\n";
#  $optstrb = " " . $optstrb;
#  $optstra = chop($stra) . $optstra;
#  $optdisplay = " " . $optdisplay;
#}
#if ($debug) { print "o:  $optalign \n"; print "a: $optstra\n"; print "d: $optdisplay\n"; print "b: $optstrb\n";}
#
}


sub topcorner {
 my($i, $j, $width) = @_;
 
 if ($i == 1) { return 1; }
 
 if ($i <= $width && $i == ($j+1)) { return 1; }
 
 return 0;
}


sub getAlignStats {
  # refStr - the reference (correct) string from which orig and mod
  #    are aligned against.  origStr - the original (unmodified)
  #   string.  This should be the erroneous string before error
  #   correction.  
  # modStr - the modified string (post error
  #   correction).  
  # *path  - paths in the {orig,mod} versus ref
  #   alignment.  
  # *max   - positions in the {orig,mod} string where the
  #   local alignment ended.

  my ($refStr, $origStr, $modStr, $origPath, $modPath, $origI, $modI) = @_;
  my $refPos;
  my $errors = 0;
  my $corruptions = 0;
  my ($s, $o, $m);
  my $modRefI;
  my $origRefI;
  my $step;
  # Start the comparison of original and modified strings at the same position.
  $origI = $origI - 1;
  $modI  = $modI  - 1;
  if ($origI > $modI) {
    $refPos = $modI;
    ($origPath, $origRefI, $origI) = advanceToPos($origPath, $refPos, $origI, $origI);
    $modRefI = $refPos;
  }
  else {
    $refPos = $origI;
    ($modPath, $modRefI, $modI) = advanceToPos($modPath, $refPos, $modI, $modI);
    $origRefI = $refPos;
  }
  print "starting out at : $origRefI, $modRefI \n";
  printf("len op: %d  len mp: %d\n", length($origPath), length($modPath));
  # look for errors by comparing the aligned positions.
  while (length($origPath) > 0 && length($modPath) > 0) { 
    # look for errors left in error correction.
    $s = substr($refStr, $refPos, 1); 
    $m = substr($origStr, $refPos, 1);
    print "original string at $refPos: $s\n";
    printf("mod string at: $modI : %s\n", substr($modStr, $modI, 1));
    printf("ref string at: $origI: %s\n", substr($origStr, $origI, 1));
    if ($s ne substr($modStr, $modI, 1)) {
      $errors = $errors + 1;
      print "error at $refPos\n";
      # look for data corruption (correct to incorrect)
      if ($s eq substr($origStr, $origI, 1)) {
	$corruptions = $corruptions + 1;
	print "corruption at $refPos\n";
      }
    }
    else {
      print "all ok\n"; 
    }
    $refPos = $refPos - 1;                       # path      ref-dest ref-cur      #other-cur 
    $step = "";
    while ($step ne "\\" && length($origPath) > 0) {
      ($origPath, $origRefI, $origI) = advanceToPos($origPath, $refPos, $refPos + 1, $origI);
    }
    $step = "";
    while ($step ne "\\" && length($modPath) > 0) {
      ($modPath,  $modRefI,  $modI)  = advanceToPos($modPath,  $refPos, $refPos + 1, $modI);
    }
    print "new positions: $origI, $modI\n";
    if ($origRefI != $modRefI) {
      print "somehow $origRefI != $modRefI\n";
    }
  }
  return ($errors, $corruptions);
}



sub advanceToPos {
  # advance pointer in reference string from curPos to tarPos.  str is
  # dependent on the alignment of ref, so there is no target position
  # for str.  Still, it needs a current position that will be traced
  # along the alignment path until refCurPos = refTarPos
  my ($pathStr, $refTarPos, $refCurPos, $strCurPos) = @_;

  print "ATP: advancing to $refTarPos from $refCurPos\n";
  my $step;
  while($refCurPos > $refTarPos && $refCurPos >= 0 && $strCurPos >= 0 && length($pathStr) > 0) { 
    $step = chop($pathStr);
    print "ATP:  got step: $step\n";
    if ($step eq "\\" || $step eq "/") {
      $refCurPos = $refCurPos - 1;
      $strCurPos = $strCurPos - 1;
    }
    elsif ($step eq "|") {
      $refCurPos = $refCurPos - 1;
    }
    elsif ($step eq "-") {
      $strCurPos = $strCurPos - 1;
    }
  }
  if ($refCurPos  < 0 || $strCurPos < 0) {
    print "$refCurPos < 0, or $strCurPos < 0, somehow there were problems!\n";
  }
  return ($pathStr, $refCurPos, $strCurPos, $step);
}

sub compareAlignments {
  my($origAlign, $modAlign, $endo, $endm) = @_;
  my $len;
  my ($leno, $lenm) = (length($origAlign),length($modAlign));
  my ($poso, $posm) = ($leno - 1, $lenm - 1);
  my ($skipo, $skipm, $skipped);
  my ($errors, $corruptions) = (0,0);
  # Make sure the two sequences begin by comparing the same
  # nucleotide.
  while (substr($origAlign, $poso, 1) eq "-") {
    $poso = $poso - 1;
  }
  while (substr($modAlign, $posm, 1) eq "-") {
    $posm = $posm - 1;
  }

  if ($poso > $posm) {
#    printf("backing up original by %d\n", $poso - $posm);
    ($poso, $skipo) = backup($origAlign, $poso, $poso - $posm);
  }
  elsif ($posm > $poso) {
#    printf("backing up mod by %d\n", $posm - $poso);
    ($posm, $skipm) = backup($modAlign, $posm, $posm - $poso);
  }
#  print "starting at $poso, $posm\n";
  # determine when the alignments start, since it is possible the
  # local alignment skipped part of the beginning.
  my($starto, $startm);
  $starto = index($origAlign, "\\");
  $startm = index($modAlign,  "\\");

  while ($poso >= $starto && $posm >= $startm) {
    if (substr($modAlign, $posm, 1) ne "\\") {
      $errors = $errors + 1;
      if (substr($origAlign, $poso, 1) ne substr($modAlign, $posm, 1)) {
	$corruptions = $corruptions + 1;
      }
    }
    ($poso, $skipped) = backup($origAlign, $poso, 1);
    $skipo = $skipo + $skipped;
    ($posm, $skipped) = backup($modAlign, $posm, 1);
    $skipm = $skipm + $skipped;
  }
  return ($errors, $corruptions);
}

sub backup {
  my($alignStr, $pos, $len) = @_;
  my $i;
  if (length($alignStr) < $len) {
    return (-1, 0);
  }
  my $skipped = 0;
  for $i (reverse(1 .. $len)) { 
    # skip past any insertions.
    while (substr($alignStr, $pos, 1) eq "-" && $pos >= 0) {
      $pos = $pos - 1;
      $skipped = $skipped + 1;
    }
    $pos = $pos - 1;
  }
  return ($pos, $skipped);
}

