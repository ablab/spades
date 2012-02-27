#!/usr/bin/env perl

if ($#ARGV < 1) {
  print "usage: $0  traceSequenceFile blatFile [minscore] [maxId] [matesId]  \n";
  exit(0);
}



# parse all sorts of input.
$traceFasFileName = shift @ARGV;
$blatFileName     = shift @ARGV;

$minScore = 0;
if ($#ARGV >= 0) {
  $minScore         = shift @ARGV;  # if blat can't do it, I'll filter myself!
}

$maxId = 0;
$matesId = 0;
if ($#ARGV >= 0) {
  $maxId   = shift @ARGV;
}
if ($#ARGV >= 0) {
  $matesId = shift @ARGV;
}
$lbstar = 9999999999999999;
$ubstar = 0;
if ($#ARGV >= 0) {
  $lbstar = shift @ARGV;
}

if ($#ARGV >= 0) {
  $ubstar = shift @ARGV;
}

open(TFF, "$traceFasFileName") or die "cannot open $traceFasFileName\n";

@traceFasFile = <TFF>;

@titles = grep(/^>.*/ , @traceFasFile);

%names = {};
%mates = {};
%mate_names = {};
%directions = {};
@ids = ();
foreach $title (@titles) {
  if ($title =~ /^>gnl\|ti\|(\S*) name:(\S+) mate:(\S+) mate_name:(\S+) .* end:([FR])/) {
    push @ids, $1;
    $names{$1} = $2;
    $mates{$1} = $3;
    $mate_names{$1} = $4;
    $directions{$1} = $5;
  }
}


open(BF, "$blatFileName" ) or die "cannot open $blatFileName\n";

# store all hits in a hash of (of arrays).
%scores = {};
while (<BF>) {
  $line = $_;      #name   contig   pct    len     mismch   gap    q.s     q.e     sub.s    sub.e  eval    score
  $line =~ /gnl\|ti\|(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)/;
  @entry = ($2, $3, $4, $5, $6, $7, $8,  $9, $10, $11, $12);
  $key = $1;
#  print "$key  $9 $10 \n";
#  print "@entry\n";
  if (! exists($scores{$key})) {
    $scores{$key} = ();
  }
  push @{$scores{$key}}, [ @entry ]    ;
}



# look at mate-pair infomration to see if there is anything fantastic about the mate pair distances.

foreach $id (@ids) {
  $mateid  = $mates{$id};
  $dir     = $directions{$id};
  $matedir = $directions{$mateid};
  # There are probably several hits for this id
  $numid = $#{$scores{$id}};
  if ($numid > $maxId) {
    $numid = $maxId;
  }
  $nummate = $#{$scores{$mateid}};
  if ($nummate > $matesId) {
    $nummate = $matesId;
  }
  print "$id ($dir) - $mateid: ($dir) ";
  for $hitIndex (0.. $numid) {
    if (${$scores{$id}}[$hitIndex][10] > $minScore) {
    # Look at all the positions of the mates
    print "(";
    for $mateHitIndex (0 .. $nummate) {
      $if = ${$scores{$id}}[$hitIndex][7];
      $ib = ${$scores{$id}}[$hitIndex][8];
      $mf = ${$scores{$mateid}}[$mateHitIndex][7];
      $mb = ${$scores{$mateid}}[$mateHitIndex][8];
 
      if ($ib < $if) {
    	  $tmp = $ib;
    	  $ib  = $if;
    	  $if  = $tmp;
    	}
      if ($mb < $mf) {
    	   $tmp = $mb;
    	   $mb  = $mf;
    	   $mf  = $tmp;
      }


 $dist = 0;
 if ($ib < $mf) {
   $dist = $mf - $ib;
   if ($mf > $lbstar and $ib < $ubstar) {# and $ib > (31942740 - 2000)) {
     print "*";
   }
 }
  
 if ($mb < $if) {
   $dist = $if - $mb;
   if ($if > $lbstar and $if < $ubstar) {# and $ib > (31942740 - 2000)) {(31942823+2000) and $mb > (31942740 - 2000)) {
     print "*";
   }
 }
     

#      if ($dir eq "F") {
#	 $dist = ${$scores{$mateid}}[$mateHitIndex][7] - ${$scores{$id}}[$hitIndex][8];
#  }
#      else {
#	$dist = ${$scores{$mateid}}[$hitIndex][7] - ${$scores{$id}}[$mateHitIndex][8];
#}
      print "$dist $if-$ib $mf-$mb  ";
}
    print ")";
}
}
  print "\n";
}
