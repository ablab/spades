#!/usr/bin/env perl

use POSIX;

$in = shift @ARGV;

$seq = "";
open (IN, "$in") or die "cannot open $in\n";
# discard fasta title
$title = <IN>;

chomp($title);
while (<IN>) {
  $line = $_;
  chomp $line;
  $seq = $seq . $line;
}

if ($#ARGV >= 0) {
  $sdf = shift @ARGV; # start divergence after fraction
}
else {
  $sdf = 1.0;
}
if ($#ARGV >= 0) {
  $tdf = shift @ARGV; # totally divergent after fraction
}
else {
  $tdf = 1.0;
}

if ($#ARGV >= 0) {
  $od = shift @ARGV;  # original divergence
}
else {
  $od = 0;
}


$len = length($seq);


$sdp = POSIX::floor($len * $sdf);
$tdp = POSIX::floor($len * $tdf);
$pos = 0;
$newseq = "";
@nucs = ('a', 'c', 't', 'g');
$errp = $od;
$gapc = 0.50;
$pctind = 0.10;
$gapbias = 0.5;
$mat = 0;
$mut = 1;
$gap = 2;
$del = 3;

$pe = $mat; # no previous error
$errp = ((4/3)*(1-$pctind) + ($pctind))*$errp; # adjust error frequency for mut to original nucleotide
#print "using: $errp $gapBias \n";
while ($pos < $len) {
  if ($pe == $mat || $pe == $mut || (rand() > $gapc)) {
    if (rand() < $errp) {
      if (rand() < $pctind) {
	if (rand() < $gapbias) {
	  $pe = $gap;
	  $rp = POSIX::floor(rand(4));
#	  print "gapping:  $rp  @nucs[$rp]\n";
	  $newseq = $newseq . @nucs[$rp];
	}
	else {
	  $pe = $del;
	  $ch = substr($seq, $pos, 1);
#	  print "deleting: $ch\n";
	  $pos++;
	}
      }
      else {
	# mutate sequence at pos:
	$pe = $mut;
	$ch = substr($seq, $pos, 1);
	$rp = POSIX::floor(rand(4));
#	print "mutating $ch :$rp  @nucs[$rp]\n";
	$newseq = $newseq . @nucs[$rp];
	$pos++;
      }
    }
    else {
      $pe = $mat;
      $newseq = $newseq . substr($seq, $pos, 1);
      $pos++;
    }
  }
  else {
    if ($pe == $gap) {
      $newseq = $newseq . @nucs[POSIX::floor(rand(4))];
    }
    else {
      $pos++;
    }
  }
#  print "$tdp $sdp $pos\n";
  if ($pos >=  $sdp && $pos < $tdp && $sdp < $tdp) {
    $frac = ($pos - $sdp) / ($tdp - $sdp);
#    print "$pos $errp $frac\n";
    $errp =  $od + (1-$d) * $frac ;
#    print "errp: $errp\n";
  }
}

print "$title (corrupt $sdp $tdp)\n";
while (length($newseq) > 0) {
  $sub = substr($newseq, 0, 60);
  $newseq = substr($newseq, 60);
  print "$sub\n";
}
