#!/usr/bin/env perl


if ($#ARGV <0 ) {
  print "usage: AssembleAndCorrect.pl readsFile \n";
  print "            [-matesFile matesFile] [-vertexSize v] [-exeDir e]\n";
  exit(0);
}
$vertexSize = 20;
$machtype = $ENV{"MACHTYPE"};
$srcDir = $ENV{"EUSRC"};
$exeDir = "$srcDir/assembly/$machtype";

$readsFile = shift @ARGV;
$matesFile = "";

while ($#ARGV >= 0) {
  $option = shift @ARGV;
  if ($option eq "-matesFile") {
    $matesFile = shift @ARGV;
  }
  if ($option eq "-vertexSize") {
    $vertexSize = shift @ARGV;
  }
}
system("$srcDir/assembly/assemble.pl $readsFile -vertexSize $vertexSize");

system("$srcDir/assembly/$machtype/simplifyGraph $readsFile $readsFile.simple -minComponentSize 500   -minEdgeLength 50 -removeBulges 80 -removeLowCoverage -vertexSize $vertexSize");
system("$srcDir/assembly/$machtype/joinsas $readsFile.simple 10 $readsFile.simple.j -vertexSize $vertexSize");
#system("dot -Tps $readsFile.simple.dot -o $readsFile.simple.ps; ps2pdf $readsFile.simple.ps");
system("$srcDir/assembly/ReorderAndLink.pl $readsFile.simple.j $readsFile -vertexSize $vertexSize");
if ($matesFile ne "") {
  system("$srcDir/euler/euler_et -s $readsFile.simple.j.r -x $vertexSize -o $readsFile.simple.r.et.dot -S ");
  system("$srcDir/euler/euler_db -i $readsFile.simple.j.r -x $vertexSize -o $readsFile.simple.r.db.dot -X -M 0 -m $matesFile");
}
else {
  system("$srcDir/euler/euler_et -s $readsFile.simple.j.r -x $vertexSize -o $readsFile.simple.r.et.dot -S -X");
}
