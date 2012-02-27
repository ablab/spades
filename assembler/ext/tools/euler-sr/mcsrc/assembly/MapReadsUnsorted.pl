#!/usr/bin/env perl
if ($#ARGV < 1) {
		print "usage: MapReadsUnsorted readsFile blastFile (tabular)\n";
		exit(1);
}

$readsFile = shift @ARGV;
$blastFile   = shift @ARGV;


%titles = {};
%pos    = {};


open(RF, "$readsFile") or die "cannot open $readsFile\n";
$title = "";
$seqLen = 0;
@titles = ();
@lengths = ();
while(<RF>) {
		$line = $_;
		if ($line =~ /^>(\S+)/) {
				$title = $1;
				if ($title ne "") {
						push @titles, $title;
						push @lengths ,$seqLen;
						$seqLen = 0;
						$title = $1;
				}
		}
		else {
				chomp $line;
				$ll  = length($line);
				$seqLen += length($line);
		}

}
if ($title != "") {
		push @titles, $title;
	  push @lengths, $seqLen;
}


open(BF, "$blastFile" ) or die "cannot open $blastFile\n";

while (<BF>) {
		$line = $_;
		chomp $line;
		if ($line =~ /(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)/) {
				$hitTitle = $1;
				$qStart = $7;
				$qEnd   = $8;
				$sStart = $9;
				$sEnd   = $10;
				$strand = 0;
				if ($sEnd < $sStart) {
#						print "strand 1!!!\n";
						$strand = 1;
						$temp = $sStart;
						$sStart = $sEnd;
						$sEnd   = $temp;
				}
#				print "parsed $hitTitle $qStart $qEnd $sStart $sEnd $strand\n";
#				@{$pos{$hitTitle}} = ($qStart, $qEnd, $sStart, $sEnd, $strand);
				print "$sStart $sEnd $hitTitle $strand\n";
		}
}
#print "done reading blast\n";



#for ($s = 0; $s <= $#titles; $s++) {
##		print "checkin $titles[$s] $lengths[$s]\n";
#		if (exists $pos{$titles[$s]}) {
#				$refStart = @{$pos{$titles[$s]}}[2];
#				$refEnd = @{$pos{$titles[$s]}}[3];
#				if (exists($pos{$title})) {
#						print "$refStart $refEnd $titles[$s] @{$pos{$titles[$s]}}[4]\n";
#				}
#		}
#}

