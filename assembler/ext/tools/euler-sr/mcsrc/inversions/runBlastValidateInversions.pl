#!/usr/bin/env perl
print "startint rbfi\n";
@dirs = ();
while ($#ARGV >= 0) {
  push @dirs, shift @ARGV;
}

$srcDir  = "~/projects/mcsrc";
$compDir = "$srcDir/comparative/powerpc";
$invDir  = "$srcDir/inversions";
$cwd = $ENV{"PWD"};
$dbName = "july_alignments";
$seqName = "";
foreach $dir (@dirs) {
#  $runEs = "~/projects/encode/utils/runExtractSeq.pl $dir";
#  system($runEs);

  #`pushd $dir`;
  $sequenceDir = "$cwd/$dir";

  $invFile = "~/projects/inversions/ENCODE/$dir/$dir.new.inversions";
  $seqName = $dir;
  chdir("$dir");
  $socket = "";
  if (exists  $ENV{"GSCK"}) {
    $socket = $ENV{"GSCK"};
  }

#  `cat *.ref *.qry > $dir.sequences`;
#  `formatdb -p F -i $dir.sequences`;

#  $command = "$invDir/printCoordinates.pl $dbName $invFile $seqName $dir.ortho . $socket \n";
#     print "running $command\n";    
#  system($command);
  `mkdir bins`;
#  $command = "$invDir/blastValidateInversions.pl $invFile . $dir.sequences $dir $dir.ortho bins/bin $socket";
#  print "runing '$command'\n";
#  system("$command");
  chdir("bins");
  @binFiles = glob("*.fasta");
  @binBase = ();
  foreach $bf (@binFiles) {
    $bf =~ /(.*)\.[^\.]+\Z/;
    $base = $1;
    #      push @binBase, $base;
    print "bf: $bf base: $base\n";
#    system("$invDir/checkBin.pl $dbName $base.fasta ../sequences.txt $dir $sequenceDir $sequenceDir/blastdb ../$dir.ortho $base.confirmed.fa 0.001 -d ../ENm001.data");
  }

  $command = "$invDir/makeAllBins.pl \"*.confirmed.fa\" ../sequences.txt ../$dir.ortho ../ $dir.2.chars";
  system($command);
  chdir("../..");
}

