#!/usr/bin/env perl

if ($#ARGV < 1) {
  print "Usage: MapReads.pl readTitlesFile blastTabFile\n";
  print "Both the read titles file and the blast tab file should be sorted\n";
  print "lexicographically.  The blasttabfile should only contain the \n";
  print "highest hit for each read.\n";
  exit 0;
}
$readsFile = shift @ARGV;
$blastIn   = shift @ARGV;


open(RF, "$readsFile") or die "cannot open $readsFile\n";
open(BI, "$blastIn") or die "cannot open $blastIn\n";


$done = 0;
$title = "init";
$searchBlastMatch = 0;
$searchReadTitle = 0;
$fetchReadLine = 1;
$fetchBlastLine = 1;
while ($readsDone == 0 && $blastDone == 0) {
  if ($fetchReadLine == 0 && $fetchBlastLine == 0) {
    print "ERROR! I'm supposed to read at least one line!\n";
    exit(1);
  }
#  if ($fetchReadLine == 0) {
#    print "not trying a read this time around\n";
#  }
  $searchReadTitle = 1;

  # if we are suposed to search for a read, look for it
  # break when no more input is available ($readsDone == 1)
  if ($fetchReadLine == 1) {
    if (($read = <RF>)) {
      #print "checking read $read\n";
      if ($read =~ />(\S+)/) {
	$title = $1;
	$title =~ /gnl\|ti\|(\d+)/;
	$titleNumber = $1;
#	print "found title: $titleNumber\n";
	$searchReadTitle = 0;
      }
    }
    else {
#      print "finished processing reads\n";
      $readsDone = 1; 
    }
  }
  # Stop processing if the read input is empty
  if ($readsDone == 1) {
#    print "done processing read\n";
    next;
  }
  # advance to the next blast hit
  if ($fetchBlastLine == 1) {
    if (($blastLine = <BI>)) {
      #print "checked blast line : $blastLine\n";
      if ($blastLine =~/(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)/) {
	$hitTitle = $1;
	$qstart = $7;
	$qend   = $8;
	$sstart = $9;
	$send   = $10;
	$strand = 0;
	if ($send < $sstart) {
	  $strand = 0;
	  $temp = $sstart;
	  $sstart = $send;
        $send   = $temp;
	}
	$hitTitle =~ /gnl\|ti\|(\d+)/;
	$hitTitleNumber = $1;
#	print "got new blast title: $hitTitle\n";
      }
      else {
	print "ERROR Parsing blast line!\n";
	print "The offending line is:\n";
	print $blastLine;
	exit(0);
      }

    }
    else {
#      print "done reading blast\n";
      $blastDone = 1;
    }
  }
  # stop processing if the blast input is finished
  if ($blastDone == 1) {
    next;
  }
#  print "will compare '$title' '$hitTitle'\n";
  # Now compare the two titles.  If they are the same, print the 
  # mapping.  Otherwise, set a file to advance in the file that 
  # has the lower read title.
  # if the blast hit is the first one for this read, print it
  if ($hitTitle eq $title) {
    print "$title $sstart $send\n";
    $searchBlastMatch = 0;
    $fetchBlastLine = 1;
    $fetchReadLine  = 1;
  }
  elsif ($hitTitle lt $title) {
    # the blast hit is greater than the read, so we need to 
    # advance in reads before blasts
#    print "blast went ahead, $title $hitTitle\n";
    $fetchReadLine = 0;
    $fetchBlastLine  = 1;
  }
  elsif ($title lt $hitTitle) {
    # The read title is greater than the blast title.
    # Advance blast hits instead
#    print "read went ahead, $title, $hitTitle\n";
    $fetchReadLine = 1;
    $fetchBlastLine = 0;
  }
}  
