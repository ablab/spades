#!/usr/bin/env perl
use POSIX;

if ($#ARGV != 9) {
  print "usage: SimulateAssembly.pl genome readlength cloneLength coverage muationRate indelRate cloneUncertainty workDir resultsDir jobid\n";
  exit(0);
}
$argstr = join(" ", @ARGV);

$genome    = shift @ARGV;
$readLength = shift @ARGV;
$cloneLength = shift @ARGV;
$coverage   = shift @ARGV;
$mutationRate = shift @ARGV;
$indelRate    = shift @ARGV;
$cloneUncertainty = shift @ARGV;
$workDir   = shift @ARGV;
$resultsDir = shift @ARGV;
$jobID = shift @ARGV;


$outputDir = "$workDir.$$";

open(RESULT, ">$resultsDir/$jobID.result");
open(LOG, ">$resultsDir/$jobID.log");
$time = `date`;
print LOG "Running Simulate assembly at $time\n";
print LOG "genome: $genome\n";
print LOG "readLength: $readLength\n";
print LOG "cloneLength: $cloneLength\n";
print LOG "coverage: $coverage\n";
print LOG "mutationRate: $muationRate\n";
print LOG "indelRate: $indelRate\n";
print LOG "cloneUncertainty: $cloneUncertainty\n";
print LOG "workDir: $workDir\n";
print LOG "resultsDir: $resultsDir\n";
print LOG "jobID: $jobID\n";

if (! -e $genome) {
  print "Genome file not found! $genome\n";
  exit(0);
}

if (! -e $outputDir) {
  print LOG "making dir: $outputDir\n";
  `mkdir -p $outputDir`;
}


chdir ($outputDir);
print LOG "working dir: $outputDir\n";
# Link the genome here for reference
if ($genome =~ /\/([^\/]+)$/) {
  $localGenome = $1;
}
else {
  $localGenome = $genome;
}

`ln -s $genome ./$localGenome`;

$reads = "$localGenome.reads";
$mates = "$reads.mates";

$mach   = $ENV{"MACHTYPE"};
$asmexe = $ENV{"EUSRC"} . "/assembly/$mach";
$asmscr = $ENV{"EUSRC"} . "/assembly";
$eulexe = $ENV{"EUSRC"} . "/euler/";



# First create the reads
$cmd = "$asmexe/readsimulator $genome $reads -rl $readLength -cloneLib $cloneLength 0 10 -mutate $mutationRate -indel $indelRate -pairFile $mates -coverage $coverage";
RunCmd($cmd);


# Now fix errors
$cmd = "$asmscr/FixErrorsSerial.pl $reads $coverage $outputDir";
RunCmd($cmd);



# Now assemble
$fixed = "$reads.fixed";
$cmd = "$asmscr/assemble.pl $fixed -dir $outputDir";
RunCmd($cmd);

# Simplify the graph
$cmd = "$asmexe/simplifyGraph $fixed $fixed.simple.1 -minComponentSize 500 -minEdgeLength 80 -removeSimpleBulges 80";
RunCmd($cmd);
$cmd = "$asmexe/simplifyGraph $fixed.simple.1 $fixed.simple.2 -removeBulges 120";
RunCmd($cmd);
$cmd = "$asmexe/simplifyGraph $fixed.simple.2 $fixed.simple.3 -removeLowCoverage 5";
RunCmd($cmd);

# prep the reads for equivalent transformation
$cmd = "$asmscr/ReorderAndLink.pl $fixed.simple.3 $fixed  -dir $outputDir";
RunCmd($cmd);

$readsr = "$fixed.simple.3.r";

# run equivalent transformation on the resulting graph 
$cmd = "$eulexe/euler_et -s $readsr -x 20 -o $readsr.et.dot";
RunCmd($cmd);

# now, if we ran euler_db on all mate-pairs, that would be way too many.
# figure out how many mate pairs there are, and just use the first 10%

$matesReduced = "$mates.reduced";
$nmates = `wc -l $mates`;
$fewer = POSIX::floor($nmates * 0.10);
$cmd = "head -$fewer $mates > $matesReduced";
RunCmd($cmd);

# create the rule file
$ruleFile = "name.rul.$cloneUncertainty";

$cmd = "$asmscr/PrintNameRul $cloneLength $cloneUncertainty $readLength > $ruleFile";
RunCmd($cmd);

# Run euler_db
$cmd = "$eulexe/euler_db -i $readsr -x 20 -o $readsr.db.dot -w 2 -m $matesReduced -r $ruleFile ";
RunCmd($cmd);


# Now, if everything worked, gather some summary statistics.
$nETgt500=  `pcl $readsr.et.contig -minLength 500 | wc -l`;
chomp $nETgt500;
$nDBgt500 = `pcl $readsr.db.contig -minLength 500 | wc -l`;
chomp $nDBgt500;
$ETN50 = `pcl $readsr.et.contig | N50`;
chomp $ETN50;
$DBN50 = `pcl $readsr.db.contig | N50`;
chomp $DBN50;
$ndbGraphEdge = `grep -c ">" $fixed.edge`;
chomp $ndbGraphEdge;
$nSimple1GraphEdge = `grep -c ">" $fixed.simple.1.edge`;
chomp $nSimple1GraphEdge;
$nSimple2GraphEdge = `grep -c ">" $fixed.simple.2.edge`;
chomp $nSimple2GraphEdge;
$nSimple3GraphEdge = `grep -c ">" $fixed.simple.3.edge`;
chomp $nSimple3GraphEdge;
$nFixed = `grep -c ">" $fixed`;
chomp $nFixed;
$nOrig  = `grep -c ">" $reads`;
chomp $nOrig;

print RESULT "$jobID, SUCCESS, $nDBgt500, $DBN50, $nETgt500, $ETN50, $ndbGraphEdge, $nSimple1GraphEdge, $nSimple2GraphEdge, $nSimple3GraphEdge, $nFixed, $nOrig, $argstr\n";
print "made it to end, worked\n";
`rm -rf $workDir.*`;

sub RunCmd {
  my ($cmd) = @_;
  `$cmd`;
  $status = $?;
  print LOG "$status $cmd\n";
  if ($status != 0) {
    $host = $ENV{"HOST"};
    $dir  = $workDir;
    print RESULT "$jobID FAILED at $host, $dir, command $cmd : $argstr\n";
    exit(0);
  }
}





