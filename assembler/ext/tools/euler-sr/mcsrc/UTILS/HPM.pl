#!/usr/bin/env perl
use strict;
use POSIX;

if ($#ARGV != 4) {
  printf("usage: $0 infile starterrors poly_bias random_error_rate outfile\n");
  exit(0);
}
my ($infile, $startErrors, $miscountProb, $randomErrorProb, $outfile) = @ARGV; 

open(INFILE, "$infile") or die "cannot open $infile\n";
open(OUTFILE, ">$outfile") or die "cannot open $outfile\n";

my %polyACount = ();
my %polyCCount = ();
my %polyTCount = ();
my %polyGCount = ();

my @reads = ();
my $read = "";
my ($pa, $pc, $pg, $pt);
my $line;
my $title;
$title = <INFILE>; # get rid of the first line
chomp $title;

my ($ea, $ec, $et, $eg, $rand);
my $netErrors = 0;
while (<INFILE>) {
  chomp;
  $line = lc($_);
  if ($line =~ />/) { 
    ($read, $ea) = ProcessRead($read, "a", $startErrors, \%polyACount, 1);
    ($read, $ec) = ProcessRead($read, "c", $startErrors, \%polyCCount, 1);
    ($read, $et) = ProcessRead($read, "t", $startErrors, \%polyTCount, 1);
    ($read, $eg) = ProcessRead($read, "g", $startErrors, \%polyGCount, 1);
    ($read, $rand) = IntroduceRandomErrors($read, $randomErrorProb);
    $netErrors = $netErrors + $ea + $ec + $et + $eg + $rand;
    print OUTFILE "$title\n";
    while (length($read) > 0) {
      print OUTFILE substr($read, 0, 60);
      print OUTFILE "\n";
      $read = substr($read, 60);
    }
    $title= $line;
    $read = "";
    next;
  }
  $read = $read . $line;
}
($read, $ea) = ProcessRead($read, "a", $startErrors, \%polyACount, 1);
($read, $ec) = ProcessRead($read, "c", $startErrors, \%polyCCount, 1);
($read, $et) = ProcessRead($read, "t", $startErrors, \%polyTCount, 1);
($read, $eg) = ProcessRead($read, "g", $startErrors, \%polyGCount, 1);
($read,$rand) = IntroduceRandomErrors($read, $randomErrorProb);

$netErrors = $netErrors + $ea + $ec + $et + $eg;

 print OUTFILE "$title\n";
    while (length($read) > 0) {
      print OUTFILE substr($read, 0, 60);
      print OUTFILE "\n";
      $read = substr($read, 60);
    }

print "introduced $netErrors errors\n";
my $key;

for $key (sort {$a <=> $b} keys %polyACount) {
  print "$key\t$polyACount{$key}\n";
}
for $key (sort {$a <=> $b} keys %polyCCount) {
  print "$key\t$polyCCount{$key}\n";
}
for $key (sort {$a <=> $b} keys %polyGCount) {
  print "$key\t$polyGCount{$key}\n";
}
for $key (sort {$a <=> $b} keys %polyTCount) {
  print "$key\t$polyTCount{$key}\n";
}


sub ProcessRead {
  my ($read, $nuc, $thresh, $counts, $modify) = @_;
  my $pos = 0;
  my $lenMatch = 0;
  my $newlen;
  my $newstr;
  my $matchstr;
  my $polystr;
  my $oldread = $read;
  my $modified = 0;
  my $netErrors = 0;
  while ($pos >= 0 && $pos < length($read)) {
    $lenMatch = 0;
    $matchstr = "$nuc\{$thresh,\}";
    if (substr($read, $pos) =~ /($matchstr)/) {  
      $polystr = $1;
      $lenMatch = length($polystr);
      $\{$counts}{$lenMatch} =  $\{$counts}{$lenMatch} + 1;
      $pos = index($read, $polystr, $pos);
      $newlen = $lenMatch; 
      
      if ($modify) {
	$modified = 1;
	$newlen  = NewLength($lenMatch, $thresh,1 - $miscountProb );
	$netErrors = $netErrors + abs($lenMatch - $newlen);
	$newstr = BuildStr($nuc, $newlen);
	# only change the first occurrence of the string
	$matchstr = "$nuc\{$lenMatch\}";
	substr($read,$pos) =~ s/$matchstr/$newstr/;
#	print "changed $matchstr \n";
#	print "        $newstr\n";
      }
      $pos = $pos + $newlen;
    }
    else { 
      $pos = -1;
    }
  }
  return ($read, $netErrors);
}

sub BuildStr {
  my ($template, $mult) = @_;
  my $i;
  my $result = "";
  for $i (0..$mult-1) {
    $result = $result . $template;
  }
  return $result;
}
sub NewLength {
  my($oldwidth, $threshold, $bias) = @_;
  my $len = $oldwidth - $threshold;
  # for now use even distribution
  my $i;
  my $count = 0;
  my $stepleft = 0;
  my $stepright = 0;
  for $i (1.. $len) {
    if (rand(1) < $bias) {
      $stepleft = $stepleft + 1;
    }
    else {
      $stepright = $stepright + 1;
    }
  }
  return $oldwidth + $stepright - $stepleft;
}

sub IntroduceRandomErrors {
  my ($string, $errorrate) = @_;
  my ($i, $len, $copy);
  $copy = $string;
  $string = reverse($string);
  my $nuc;
  
  my @actg = ();
  push @actg, 'a','c','t','g';
  my %actgh = ();
  $actgh{'a'} = 0;
  $actgh{'c'} = 1;
  $actgh{'t'} = 2;
  $actgh{'g'} = 3;
  my $err = "";
  my ($r, $num);
  my $newString =  "";
  my $numerrors = 0;
  $len = length($string);
  for $i (0.. $len-1) {
    # should introduce error
    $nuc = lc(chop($string));
    $r = rand(1);
    if ($r < $errorrate) {
      $numerrors = $numerrors + 1;
      if ($r < 0.3333) {
	# mutation
	$num = POSIX::floor(rand(3));
	$nuc = @actg[($actgh{$nuc} + 1) % 4];
	$newString = $newString . $nuc;
      }
      elsif($r < 0.6666) {
	# insertion
	$newString = $newString . $nuc;
	$nuc =  @actg[POSIX::floor(rand(4))];
	$newString = $newString . $nuc;
      }
      else {
	# deletion, don't bother appending the real char.
      }
    }
    else {
      $newString = $newString . $nuc;
    }
  }
  if ($numerrors > 0 && $len > 0) {
#    printf("error rate was : %f\n", $numerrors / $len);
#    print "$newString\n";
#    print "$copy\n";
  }
  else {
#    print "length of '$string' is <= 0\n";
  }
  return ($newString, $numerrors);
}

	




  
  
