#!/usr/bin/env perl

$sourceDir = shift @ARGV;
$destDir   = shift @ARGV;

@fileTypes = (".cpp", ".h", "Makefile");

@sourceFiles = ();
@destFiles   = ();
foreach $ft (@fileTypes) {
		@sft = glob("$sourceDir/*$ft");
		@dft = glob("$destDir/*$ft");
		while ($#sft >= 0) {push (@sourceFiles, $sft[0]); shift @sft};
		while ($#dft >= 0) {push (@destFiles, $dft[0]); shift @dft;}
}

for ($s = 0; $s <= $#sourceFiles; $s++) {
		$source = $sourceFiles[$s];
		if ($source =~ /\/([^\/]+)$/) {
				$source = $1;
				$sourceFiles[$s] = $source;
		}
}

for ($d = 0; $d <= $#destFiles; $d++ ){
		$dest = $destFiles[$d];
		if ($dest =~ /\/([^\/]+)$/) {
				$dest = $1;
				$destFiles[$d] = $dest;
		}
}
		

@sortedSource = sort @sourceFiles;
@sortedDest   = sort @destFiles;

@missingDest = ();
@missingSource = ();
while ($#sortedSource >= 0 && $#sortedDest >= 0) {
		if ($sortedSource[0] eq $sortedDest[0]) {
				shift @sortedSource;
				shift @sortedDest;
		}
		elsif ($sortedSource[0] lt $sortedDest[0]) {
				push @missingDest, $sortedSource[0];
				shift @sortedSource;
		}
		else {
				push @missingSource, $sortedDest[0];
				shift @sortedDest;
		}
}

while ($#sortedSource >= 0) {
		push @missingDest, shift @sortedSource;
}
while ($#sortedDest >= 0) {
		push @missingSource, shift @sortedDest;
}

if ($#missingDest < 0) {
		print "No files need an update.\n";
}
else {
		print "missing file are: \n";
		for $md (@missingDest) {
				print "$md\n";
		}
}
