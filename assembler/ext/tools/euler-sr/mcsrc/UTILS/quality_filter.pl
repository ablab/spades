#!/usr/bin/env perl
use strict;  
use POSIX;

if ($#ARGV < 3) {
  print "Usage $0 RAW qual threshold window\n";
  exit(-1);
}

my($readsName);
my($qualName); 
my($threshold, $window);
my($curCoverage, $simLength, $simCoverage);
my($outname, $outMatesName);

my $addErrors = 0;
my $debug = 0;
my $basesChanged = 0;
$readsName = shift @ARGV;
$qualName  = shift @ARGV;
$threshold = shift @ARGV;
$window    = shift @ARGV;

if ($#ARGV >= 0) {
  if (@ARGV[0] == "adderror") {
    $addErrors = 1;
  }
  elsif (@ARGV[0] == "noadderror") {
    $addErrors = 0;
  }
  else {
    printf("last argument should be 'adderror' or 'noadderror'\n");
    exit(-1);
  }
}
    


open(READFILE, "$readsName") || die "cannot open $readsName\n";
open(QUALFILE, "$qualName") || die "cannot open $qualName\n";

my(@reads);
my(@titles, $title);
my(@qualTitles);
my(@qualities);
$title = <READFILE>;
<QUALFILE>;

chomp $title;
push @titles, $title;

#######################################################
#  Read the input and the quality files.
#######################################################
my($line);
my($index, $length);

$index = 0;
do {
  $line = getRead(\*READFILE, \@reads, \@titles, \$index);
  $length = $length + ($#{$reads[$index-1]} + 1);
  if ($index % 100 == 0) {
    print "read $index\n";
  }
}  while ($line ne "");

$index = 0;
do {
  $line = getRead(\*QUALFILE, \@qualities, \@qualTitles, \$index);
} while ($line ne "");

close(READFILE);
close(QUALFILE);

#######################################################
# Preprocess reads.
#   - Remove tail ends of low quality.
#   - Remove masked off regions.
#######################################################

#------------------------------------------------------
#  Locate the first regions from the start and the end 
#  that are of average quality > threshold
#------------------------------------------------------
my($i, $j, $start, $end);
my(@starts, @ends);
@starts = ();
@ends = ();
my($length);
for $i (0.. $#qualities) {
  if ($i % 100 == 0 && $debug) {
    printf("trimming read%d\n", $i);
  }
  push @starts,int(getStartIndex(\@{$qualities[$i]}, $threshold, $window));
  push @ends  ,int(getEndIndex(\@{$qualities[$i]}, $threshold, $window));
  if ($debug ){ 
    print "$i has start pos: @starts[$i] end: @ends[$i]\n";
  }
  if (@starts[$i] >= @ends[$i]) {
    $#{$reads[$i]} = -1;
    $#{$qualities[$i]} = -1;
  }
  $length = $#{$reads[$i]};
  trimFront(\@{$reads[$i]}, @starts[$i]);
  trimFront(\@{$qualities[$i]}, @starts[$i]);
  trimBack(\@{$reads[$i]}, $length - @ends[$i]);
  trimBack(\@{$qualities[$i]}, $length - @ends[$i]);
}

my($start, $end);
for $i (0.. $#qualities) {
  # find regions that are mainly N's
  ($start, $end)  = findMaskedRegions(\@{$reads[$i]}, 20, 0.9);
  $length = $#{$reads[$i]};
  if ($debug) {  print "$i has start pos: $start, end: $end\n";}
  if ($start >= $end) {
    $#{$reads[$i]} = -1;
    $#{$qualities[$i]} = -1;
  }
  trimFront(\@{$reads[$i]}, $start);
  trimFront(\@{$qualities[$i]}, $start);
  trimBack(\@{$reads[$i]}, $length - $end);
  trimBack(\@{$qualities[$i]}, $length - $end);
  if ($#{$reads[$i]} > 0) {
    replaceUncertanties(\@{$reads[$i]}, \@{$qualities[$i]});
  }

}

#######################################################
#
# Print the reads with the specified quality windows.
#
#######################################################

for $i (0 .. $#reads) {
  print "@titles[$i]\n";
  for $j (0 .. $#{$reads[$i]}) {
    if (($j % 50 == 0) and ($j > 0) and ($j != $#{$reads[$i]})) {
      print "\n";
    }
    print @{$reads[$i]}[$j];
  }
  print "\n";
}


sub getRead {
  my($readfile, $reads, $titles, $index) = @_;
  my($line);
  my(@values);
  my($val);
  $line = <$readfile>;
  while ($line ne "") {
    chomp($line);
    if ($line !~ />/) {
      $line =~ s/[nN]/n /g;
      $line =~ s/[xX]/x /g;
      $line =~ s/[gG]/g /g;
      $line =~ s/[aA]/a /g;
      $line =~ s/[cC]/c /g;
      $line =~ s/[tT]/t /g;
      push(@values, split(' ',$line));
    }
    elsif ($line =~ />/) {
      push @{$titles}, $line;
      if ($#values < 0) { return $line; }


      for $val (@values) {
	push(@ {${$reads}[$$index]}, $val);
      }
      $$index = $$index + 1;
      return $line;
    }
    $line = <$readfile>;
  }

for $val (@values) {
  push(@ {${$reads}[$$index]}, $val);
}
$$index = $$index + 1;
return $line;

}

sub getStartIndex {
  my($qualities, $thresh, $window) = @_;
  
  my($i, $j, $avg);

  for $i (0 .. ($#{$qualities} - $window + 1)) {
    $avg = 0;
    for $j ($i .. ($i + $window - 1)) {
      $avg = $avg + @{$qualities}[$j];
    }
    $avg = (1.0 * $avg) / $window;
    if ($avg > $threshold ) { 
      return $i; 
    }
  }

  return $#{$qualities};
}


sub getEndIndex {
  my($qualities, $thresh, $window) = @_;
  my($i, $j, $avg,$len);
  $len = $#{$qualities}; 
  for $i (($window - 1) .. $len) {
    $avg = 0;
    for $j (($len - $i) .. ($len - $i + $window - 1)) {
      $avg = $avg + @{$qualities}[$j];
    }
    $avg = (1.0 * $avg) / $window;
    if ($avg > $threshold ) { 
      return $len - $i;
    }
  }

  return  0;
}

sub makeSimulatedReads {
  my($rawReads, $rawMatePairs, $numSimReads,
     $simLength, $simReads, $simMatePairs, $qualities) = @_;

  my($i, $j);
  my($readIndex, $mateIndex);
  my(%readsSampled, $numReads);
  %readsSampled = ();
  $numReads = $#{$rawReads} + 1;
  my($length, $offset);
  $i = 0;
  my($range);
  my(@temparray);
  while ($i < $numSimReads) {
    $readIndex = POSIX::floor(rand($#{$rawReads}));
    $length = $simLength + POSIX::floor(rand(21)) - 10;
    if ($#{@{$rawReads}[$readIndex]} < $length) {

if ($debug) {      printf("not sampling %d\n", $readIndex);}

      next;
    }
    $readsSampled{int($readIndex)} = 1;
    $offset = POSIX::floor(rand($#{@{$rawReads}[$readIndex]} + 1 - $length));

if ($debug){    printf("sim read: %d, from raw read: %d, offset: %d length: %d\n",  $i, $readIndex, $offset, $length);}

    # copy slice
    @{$simReads}[$i] = ();

    @{$simReads}[$i] =  [ @{@{$rawReads}[$readIndex]}[$offset .. ($offset + $length - 1)]];
    if ($addErrors) {
#      print "looking at @{$simReads[$i]}\n";
      introduceErrors(\@{$$simReads[$i]}, \@{$$qualities[$readIndex]}, $offset, $offset+$length-1);
    }

#    print "added simulated read: @{${$simReads}[$i]}\n";
#    print "checking mate index: ${$rawMatePairs}{$readIndex}\n";
    $i = $i + 1;

if ($debug) {    printf("Sim read sample index: %d, offset: %d, length: %d\n", $readIndex,
	   $offset, $length);}

    if (exists($$rawMatePairs{$readIndex}) &&
	$$rawMatePairs{$readIndex} > 0) {
      $mateIndex = $$rawMatePairs{$readIndex};
      if ($#{@{$rawReads}[$mateIndex]} < $length) {
if ($debug){	printf("not sampling mate %d\n", $mateIndex);}
	next;
      }

      $readsSampled{$mateIndex} = 1;
      $length = $simLength + POSIX::floor(rand(21)) - 10;
#      $offset = POSIX::floor(rand(@{$endPos}[$mateIndex] - @{$startPos}[$mateIndex] - $length));
      
      $offset = POSIX::floor(rand($#{@{$rawReads}[$mateIndex]} + 1 - $length));
      # copy slice
      @{$simReads}[$i] =  [@{@{$rawReads}[$mateIndex]}[$offset .. ($offset + $length - 1)]];
      if ($addErrors) {
	introduceErrors(\@{$$simReads[$i]}, \@{$$qualities[$mateIndex]}, $offset, $offset+$length-1);
      }

#      print "added simulated read: @{@{$simReads}[$i]}\n";
if ($debug) {    printf("Sim read sample index: %d, offset: %d, length: %d\n", $mateIndex,
	   $offset, $length);}

      $i = $i + 1;
#      print "adding $readIndex, $mateIndex to sim mate pairs\n";
      push @$simMatePairs, $readIndex, $mateIndex;
    }
  }
  for $i (0 .. $numReads)  {
    if (not exists($readsSampled{$i})) {
if($debug){	print "Not present at index: $i\n";}
      }
  } 
}
    
    
sub printReadsFile {
  my($filename, $reads) = @_;
  open(READS, ">$filename") or die "can't open $filename\n";
  my ($i, $j);
  for $i (0.. $#{$reads}) {
    printf(READS ">sim_read%d\n", $i);

#    print "index: $i $#{${$reads}[$i] } \n";
    $j = 0;
    for $j (0 .. $#{${$reads}[$i] }) {
      printf(READS "%s", ${$reads}[$i][$j]);
    if ($j % 70 eq 69) {
      print READS "\n";
    }  
    } # end for
if (  (scalar @{${$reads}[$i]}) % 70 ne 0) {
    printf(READS  "\n");
  }
  }
close(READS);
}


sub printMatesFile {
  my ($filename, $mates) = @_;
  
  open(MATES, ">$filename") or die "can't open matefile $filename\n";
  
  my ($i, $j, $len);
  $len = $#{$mates} + 1;
  $i = 0;
  while ($#{$mates} >= 0) {
    printf(MATES "sim.con%d %d %d 0\n", $i, @$mates[0], @$mates[1]);
    $i = $i + 1;
    shift @$mates;
    shift @$mates;
  }
  close(MATES);
}

sub trimFront {
  my ($array, $length) = @_;
  my($i);
  for $i (0 .. ($length-1)) {
    shift @$array;
  }
}

sub trimBack {
  my ($array, $length) = @_;
  my($i);
  for $i (0.. ($length-1)) {
    pop @$array;
  }
}


sub findMaskedRegions {
  my ($read, $window, $ratio) = @_;

  my ($i,$j, $nucCount);
  my ($start, $end);
  my ($seq);
  $start = $#{$read};  # Initialize at the end of sequence to not assume there is an ok region.
  for $i (0 .. ($#{$read} - $window) ) {
    $nucCount = 0;
    $seq = "";
    for $j ($i .. $i + $window) {
      $seq = $seq . @{$read}[$j];
    if ( @{$read}[$j] eq 'g' || 
	 @{$read}[$j] eq 'a' || 
	 @{$read}[$j] eq 'c' || 
	 @{$read}[$j] eq 't') {
        $nucCount = $nucCount + 1;
      }
   }
   if ($nucCount > ($ratio * $window)) {
     $start = $i;
     last;
   }
  }

$end = 0; # initialize at beginning of sequence to show no good region has yet been found.
for $i (($window - 1) ..($#{$read})) {
  $nucCount = 0;
  for $j (($#{$read} - $i) .. ($#{$read} - $i + $window - 1)) {
    if ( @{$read}[$j] eq 'g' || 
	 @{$read}[$j] eq 'a' || 
	 @{$read}[$j] eq 'c' || 
	 @{$read}[$j] eq 't') {
      $nucCount = $nucCount + 1;
   }
  }
  if ($nucCount > $ratio * $window) {
    $end = $#{$read} - $i + $window;
    last;
  }
}

return ($start, $end);
}

sub replaceUncertanties {
  my ($read, $quality) = @_;
  my ($i);
  my ($replace);
  my @nucs = ('g', 'a', 'c', 't');
  $i  = 0;
  while ($i <= $#{$read} ) {
    if ((${$read}[$i] ne 'a') && (${$read}[$i] ne 'c') && 
        (${$read}[$i] ne 't') && (${$read}[$i] ne 'g')) {
      if (rand(1) > 0.5) {
       splice @$read, $i, 1;
       splice @$quality, $i, 1;
      }
      else {
        ${$read}[$i] = @nucs[POSIX::floor(rand(4))];
        $i = $i + 1;
      }
    }
    else { $i = $i + 1}
  }
}
  
sub introduceErrors {
  my ($read, $qualities, $start, $end) = @_;
  
  my ($i, $j);
  my @nucs = ('g', 'a', 'c', 't');
  my $randnum;
  my $prob;
  my $nuc;
  for $i (0.. ($end - $start)) {
    $randnum = rand(1);  
    $prob = POSIX::pow(10, -0.1*@{$qualities}[$i + $start]);
    if ($randnum < $prob  ) {
      $basesChanged = $basesChanged + 1;
      ${$read}[$i] = @nucs[POSIX::floor(rand(4))];
    }
  }
}
  
