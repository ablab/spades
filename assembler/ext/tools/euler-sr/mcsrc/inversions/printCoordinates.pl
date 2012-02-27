#!/usr/bin/env perl

use Bio::SeqIO;
use Bio::Seq;

use common;
use DBI;

if ($#ARGV < 2 ) {
  print "usage: $0 dbName invFileName seqname orthoPosOut [seqdir] [socket]\n";
  print "example: $0 inv/ENm001.inversions inv ENm001.all.fasta ENm001 ENm001.orth\n";
  exit(1);
}

$dbName = shift @ARGV;
$invFileName = shift @ARGV;
$regionName  = shift @ARGV;
$orthoPosOut = shift @ARGV;
$dbDriver = "dbi:mysql:$dbName";
#if ($#ARGV >= 0) {
#  $dataDir = shift @ARGV;
#  $dbDriver .= ";data_dir=$dataDir";
#}
print STDERR "pc: $dbName $invFileName $regionName $orthoPosOut $dbDriver \n";
$seqDir = ".";
if ($#ARGV >= 0) {
  $seqDir      = shift @ARGV;
}
print "$dbName $invFileName $regionName $orthoPosOut, $seqDir, \n";
print STDERR "args remaining -1 : $#ARGV\n";
if ($#ARGV >= 0) {
  $socket = shift @ARGV;
  print "using  socket: $socket\n";
  $dbDriver .= ";mysql_socket=$socket";
}
print "driveR: $dbDriver\n";
$dbh = DBI->connect($dbDriver, "mchaisso", "", );
if ($dbh == undef) {
  print "could not connect do $dbDriver\n";
  exit(0);
}

$qry = "show tables";
$stmt = $dbh->prepare($qry); $stmt->execute();

$binDir = $ENV{"HOME"}. "/projects/mcsrc/". $ENV{"MACHTYPE"} . "/bin";


%pairs = ();
@refSpecies = ();

common::ReadInvFile($invFileName, \@refSpecies, \%pairs);
print "got ref species @refSpecies \n";

# Make a directory where all of the results will go

if (! -e "tmp") {
  `mkdir tmp`;
}

# Read in all of the sequences corresponding to the inversions.
@invList = ();
$invIndex = 0;

# invSet is a map from title to sequence.
%nameToSequence = ();
%nameToIndex  = ();
#common::ReadFiles(\@refSpecies, \%pairs, \@invList, \%nameToSequence, \%nameToIndex, $seqDir);
common::PairsToInvList(\@refSpecies, \%pairs, \@invList);

print "read $#invList files \n";
# create a list of orthologous positions to speed things up
#for ($i = 0; $i <= $#invList; $i++ ) {
#  for ($j = 0; $j <= $#{$invList[$i][0]}; $j++) {
#    print "$invList[$i][0][$j] ";
#  }
#  print ".\n";
#}
print "done with list, $#invList\n";

%orthoPos = ();
common::CreateOrthoCoordinates($dbh, \%orthoPos, \@invList, $regionName, \@refSpecies, $dbName);


open (ORTHO, ">$orthoPosOut") or die "cannot open $orthoPosOut \n";

foreach $seq (keys %orthoPos) {
  print ORTHO "$seq @{$orthoPos{$seq}} \n";
}
close ORTHO;

# find out what sequences are present:
%species = ();
foreach $invIdx (0 .. $#invList) {
#  print "considering $invList[$invIdx][0][0]  $invList[$invIdx][0][3]  \n";
  $species{$invList[$invIdx][0][0]} = 1;
  $species{$invList[$invIdx][0][3]} = 1;
}
		 
open(SEQ, ">sequences.txt");
@spec = sort keys %species;
foreach $si (0 .. $#spec-1){
  print SEQ "@spec[$si] ";
}
print SEQ "@spec[-1]\n";
close SEQ;
