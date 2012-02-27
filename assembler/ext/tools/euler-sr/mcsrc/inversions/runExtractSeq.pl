#!/usr/bin/env perl

@dirs = ();
while ($#ARGV >= 0) {
  push @dirs, shift @ARGV;
}

$cwd = $ENV{"PWD"};
$m = $ENV{"MACHTYPE"};
$compDir = "~/projects/mcsrc/comparative/$m";
foreach $dir (@dirs) {
    print "going to $dir \n";
    $sequenceDir = "$cwd/$dir";
    print "sequencedir: $sequenceDir\n";
    $invFile = "$sequenceDir/$dir.inversions";
    @files = glob("$sequenceDir/*.fasta");
#    print "got files @files\n";
    @species = ();
    @sequences = ();
    foreach $file (@files) {
      $file =~ /.*\/([^\/]+\Z)/;
      $base = $1;
      push @sequences, $1;
      $base =~ /(\A[^.]+)\..*/;
      $spec = $1;
      push @species, $1;
    }
    print "go $#files files\n";
    foreach $idx (0 .. $#files) {
      foreach $jdx (0 .. $#files) {
	if ($idx != $jdx) {
	  #
	  # Code to extract the inversion sequences.
	  #
	  $outName = @species[$idx] .".". @species[$jdx] . ".ref";
	  print "ref outName: $outName\n";
	  $command = "$compDir/extinv $invFile @species[$idx] @species[$jdx] @files[$idx] $dir/$outName -n";
	  system($command);
	  print "running $command\n";
	  $outName = @species[$idx] .".". @species[$jdx] . ".qry";
	  print "qry outName: $outName\n";
	  $command = "$compDir/extinv $invFile @species[$idx] @species[$jdx] @files[$jdx] $dir/$outName -q -n";
	  system($command);
	}
      }
#      print "got base $base\n";
    }
}
