#!/usr/bin/env perl
use FwgLib;

@inFiles = ();
$FASTA = 1;
$FASTQ = 2;
$fileType = $FASTA;
$outFile = "";
while($#ARGV >= 0) {
		$opt = shift @ARGV;
		if ($opt !~ /^\-/) {
				push @inFiles, $opt;
		}
		elsif ($opt eq "-fastq") {
				$fileType = $FASTQ;
		}
		elsif ($opt eq "-outFile") {
				$outFile = shift @ARGV;
		}
		else {
				print "unknown option: $opt\n";
				exit(1);
		}
}

if ($outFile eq "") {
		print "You must specify an outfile as '-outFile f'\n";
		exit(1);
}

$EUSRC = FwgLib::CrucialGetEnv("EUSRC");
$MACHTYPE = FwgLib::CrucialGetEnv("MACHTYPE");
$cmd = "$EUSRC/assembly/$MACHTYPE/qualityTrimmer";

# run a check on all files first.
if ($fileType == $FASTA) {
		foreach $file (@inFiles) {
				if (! -e "$file.qual") {
						print "ERROR. Each fasta file 'f.fasta' must have a corresponding quality file 'f.fasta.qual -minQual 15 -span 4'\n";
						exit(1);
				}
		}
}
@trimmedFiles = ();
foreach $file (@inFiles) {
		if ($fileType == $FASTA) {
				$execCmd = "$cmd -fastA $file -qual $file.qual -outFasta $file.trimmed -minQual 10 -span 3";

		}
		elsif ($fileType == $FASTQ) {
				# is always the case if the first if isn't true for now.
				$execCmd = "$cmd -fastQ $file -outFasta $file.trimmed";
		}
		FwgLib::RunCommand($execCmd);
		push @trimmedFiles, "$file.trimmed";
}

system("cat @trimmedFiles > $outFile");
