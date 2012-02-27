#!/usr/bin/env perl
if ($#ARGV < 0) {
  print " usage: PrintGridResults.pl prefix suffix [files] [-T fileOfFileNames]\n";
}
@files = ();
$prefix = shift @ARGV;
$suffix = shift @ARGV;
$fileOfFileNames = "";
while ($#ARGV >= 0) {
    $arg = shift @ARGV;
    if ($arg =~ /^\-/) {
       $opt = $arg;
       if ($opt eq "-T") {
          $fileOfFileNames = shift @ARGV; 
       }
    }
    else {
				push @files, $arg;
   }
}
if ($fileOfFileNames ne "") {
  open(FFN, "$fileOfFileNames") or die "cannot open filename file $fileOfFileNames\n";
  @fileNames = <FFN>;
  chomp @fileNames;
  while ($#fileNames >= 0) {
    push @files, shift @fileNames;
  }
}
%kmers = ();
%cloneSizes = ();
$prefix =~ s/\./\\\./g;
foreach $file (@files) {
		if ($file =~ /$prefix.*/) {
				if ($file =~ /.*$suffix/) {
       		if ($file =~ /$prefix\.(\d+)\.(\d+)\.$suffix*/) {
  				  $k = $1;
            $d = $2;
        		$kmers{$k}      = 1;
        		$cloneSizes{$d} = 1;
        		$res = `cat $file`;
        		chomp $res;
      	  	$results{$1}{$2} = $res; 
          }
          else {
          }
       }
   }
}

@kmerList = sort {$a <=> $b} keys(%kmers);
@cloneSizeList = sort {$a <=> $b } keys(%cloneSizes);


foreach $k (@kmerList) {
		foreach $c (@cloneSizeList) {
				if (exists $results{$k}) {
						if (exists $results{$k}{$c}) {
								print "$results{$k}{$c} ";
						}
						else {
								print "0 ";
						}
				}
				else {
						print "0 ";
				}
		}
		print "\n";
}
