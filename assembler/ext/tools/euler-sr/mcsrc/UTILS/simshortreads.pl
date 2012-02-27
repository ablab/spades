#!/usr/bin/env perl
use strict;  
use POSIX;

if ($#ARGV < 8) {
  print "Usage $0 RAW qual threshold window cur_cov sim_len sim_cov sim_seqfile sim_matefile\n";
  exit(-1);
}

my($readsName);
my($qualName); 
my($threshold, $window);
my($curCoverage, $simLength, $simCoverage);
my($outname, $outMatesName);

my $addErrors = 0;
my $debug = 1;
my $basesChanged = 0;
$readsName = shift @ARGV;
$qualName  = shift @ARGV;
$threshold = shift @ARGV;
$window    = shift @ARGV;
$curCoverage= shift @ARGV;
$simLength = shift @ARGV;
$simCoverage  = shift @ARGV;
$outname   = shift @ARGV;
$outMatesName   = shift @ARGV;

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
my(@qualities);
<READFILE>;
<QUALFILE>;

#######################################################
#  Read the input and the quality files.
#######################################################
my($line);
my($index, $length);

$index = 0;
do {
  $line = getRead(\*READFILE, \@reads, \$index);
  $length = $length + ($#{$reads[$index-1]} + 1);
#  print "got read @{$reads[$index-1]} \n";
}  while ($line ne "");

$index = 0;
do {
  $line = getRead(\*QUALFILE, \@qualities, \$index);
} while ($line ne "");

close(READFILE);
close(QUALFILE);

printf("estimated genome length: %d, total read length: %d\n", $length / $curCoverage, $length);

my($newReadCount); 
$newReadCount =   $length* $simCoverage*1.0/($curCoverage*$simLength);

#######################################################
# Preprocess reads.
#   - Remove tail ends of low quality.
#   - Remove masked off regions.
#######################################################

#------------------------------------------------------
#  Locate the first regions from the start and the end 
#  that are of average quality > threshold
#------------------------------------------------------
print "new read count: $newReadCount\n";
printf("average read length: %d\n", $length / $#reads);
my($i, $j, $start, $end);
my(@starts, @ends);
@starts = ();
@ends = ();
my($length);
for $i (0.. $#qualities) {
  if ($i % 100 == 0) {
    printf("trimming read%d\n", $i);
  }
  push @starts,int(getStartIndex(\@{$qualities[$i]}, $threshold, $window));
  push @ends  ,int(getEndIndex(\@{$qualities[$i]}, $threshold, $window));
#  print "$i has start pos: @starts[$i] end: @ends[$i]\n";
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
  if ($i % 100  ==0) {
  printf("removing masked regions: %d\n", $i);
}
  ($start, $end)  = findMaskedRegions(\@{$reads[$i]}, 20, .9);
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


printf("done preprocessing reads\n");
#######################################################
# Extract mate pair information from the reads file. Mate
#  pairs are defined by two reads that have the same initial
#  label, and 'b'/'g' second labels, or 'x'/'y' second labels.
#  >read10.g
#  >read10.b  - read10.g, read10.b are mate-pairs.
#  >read30.x
#  >read30.y  - read30.x, read30.y are mate-pairs.
#  >read50.s  - no mate pair.
#######################################################
open(READFILE, "$readsName") || die "cannot open $readsName\n";
my %mateIndices = ();

getMatePairIndices(\*READFILE, \%mateIndices);
close(READFILE);
printf("done getting mate pair indices\n");
#my $key;
#for $key (keys(%mateIndices)) {  
#  printf("%d has mate at %d\n", $key, $mateIndices{$key});
#}


#######################################################
# simulate the short reads from the longer ones.
#######################################################
my(@simReads, @simMatePairs);
@simReads = ();
@simMatePairs= ();
printf("making simulated reads\n");
makeSimulatedReads(\@reads, \%mateIndices, 
		   $newReadCount, $simLength, \@simReads, \@simMatePairs, \@qualities);

printf("total bases changed: $basesChanged\n");
printReadsFile($outname, \@simReads);
#printf("got %d sim mate pairs\n", $#simMatePairs);
#print " they are: @simMatePairs\n";
printMatesFile($outMatesName, \@simMatePairs);


sub getMatePairIndices {
  my($readfile, $indices) = @_;

  my($line, $index);
  my %labelIndices = ();
  $index = -1;
  $line = <$readfile>;
  while ($line ne "") {
    if ($line =~ />(\w+\.\w+)/) {
      $index = $index + 1;
#      print "got label: $1\n";
      $labelIndices{$1} = $index;
    }
    $line = <$readfile>;
  }
  my($key, $prevKey, $prevLabel, $curLabel, $prevType, $curType);
  $prevKey   = "";
  $prevType  = "";
  $prevLabel = "";
  foreach $key (sort(keys(%labelIndices))) {
    if ($key eq "") { next;}
if ($debug) {    print "looking at key: $key , prev: $prevKey \n";}
    $key =~ /(\w+)\.(\w)/;
    $curLabel = $1;
    $curType  = $2;
#    if ($curLabel eq $prevLabel) { print "labels are equal\n"; }
    if ($curLabel eq $prevLabel) {
      if((($curType eq 'g') &&  ($prevType eq 'b')) ||
	 (($curType eq 'y') &&  ($prevType eq 'x')))
	{
if ($debug) {	  printf("found mate pair for $key and $prevKey\n"); }
	  $$indices{$labelIndices{$key}} = $labelIndices{$prevKey};
	  $$indices{$labelIndices{$prevKey}} = $labelIndices{$key};
	}
      else {
if($debug){	print "no match for $key\n";}
	$$indices{$labelIndices{$key}} = -1;
      }
   }   
    
    $prevType  = $curType;
    $prevLabel = $curLabel;
    $prevKey   = $key;
  }
}



sub getRead {
  my($readfile, $reads, $index) = @_;
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
  
