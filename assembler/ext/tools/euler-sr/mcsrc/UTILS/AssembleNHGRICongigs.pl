#!/usr/bin/env perl

if (@ARGV < 1) {
  print "usage: assemble_contigs.pl species threshold\n";
  exit(0);
}


$species = shift @ARGV;
$threshold = 1000000;
if ($#ARGV >= 0) {
  $threshold = shift @ARGV;
}

@speciesFiles = glob "$species.*";

`mkdir $species`;

chdir($species);

for $i (0 .. $#speciesFiles-1) {
  for $j ($i+1 .. $#speciesFiles) {
    `blastz ../@speciesFiles[$i] ../@speciesFiles[$j] > $species.$i.$j.lav `;
  }
}

@lavFiles = glob "*.lav";

foreach $lavFile (@lavFiles) {

  open (LAV, "$lavFile");
  @lav = <LAV>;
  @scores = ();
  grep( (/  s (\d+)/ && push(@scores, $1)) , @lav);
  @scores = sort {$a <=> $b}  @scores;
  $maxScore= pop @scores;
  if ($maxScore > $threshold) {
    for $i (0 .. $#lav) {
      if (@lav[$i] =~ /s $maxScore/) {
	@lav[$i+1] =~ /b (\d+) (\d+)/;
	$br = $1; $bq = $2;
	@lav[$i+2] =~ /e (\d+) (\d+)/;
	$er = $1; $eq = $2;
	print "$lavFile $br $er $bq $eq\n";
      }
    }
  }
}

