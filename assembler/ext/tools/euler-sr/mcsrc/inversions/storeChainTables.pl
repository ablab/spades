#!/usr/bin/env perl

if ($#ARGV != 4) {
  print "usage: $0 filename dbname refname qryname seqname\n";
  exit(0);
}

use DBI;

$fileName = shift @ARGV;
$dbName   = shift @ARGV;
$refSeq   = shift @ARGV;
$qrySeq   = shift @ARGV;
$seqName  = shift @ARGV;

print "option fileName $fileName\n";
print "option dbName $dbName\n";
print "option refSeq $refSeq\n";
print "option qrySeq $qrySeq\n";
print "option seqName $seqName\n";

# fix the queryid in the table

$refSeq =~ /([^.]+)\.?.*/;
$refBase = $1;
$qrySeq =~ /([^.]+)\.?.*/;
$qryBase = $1;

print "rb: $refBase qb: $qryBase \n";
$tableName = "$refBase\_$qryBase\_$seqName\_chain";

$tableFileName = "$refSeq\_$qrySeq\_$seqName\_chain" . ".txt";

$dbConnect = "dbi:mysql:$dbName";
if (exists $ENV{"GSCK"}) {
  $socket = $ENV{"GSCK"};
  $dbConnect .= ";mysql_socket=$socket";
  print "connect now: $dbConnect\n";
}

$dbh = DBI->connect($dbConnect);

$qryStr = "DROP TABLE IF EXISTS $tableName  ;";
print "query: $qryStr\n";
$stmt = $dbh->prepare($qryStr); $stmt->execute();

$qryStr = "CREATE TABLE $tableName (id INT AUTO_INCREMENT PRIMARY KEY, " .
  "lavid INT, " .
  "chainid INT) " ;
print "running $qryStr\n";
$stmt = $dbh->prepare($qryStr); $stmt->execute();

$qryStr = "ALTER TABLE $tableName add index (chainid)";
$stmt = $dbh->prepare($qryStr); $stmt->execute();

$qryStr = "ALTER TABLE $tableName add index (lavid)";
$stmt = $dbh->prepare($qryStr); $stmt->execute();

print "importing\n";
`cp $fileName $tableFileName `;
print "mysqlimport $dbName $tableFileName\n";
`~/fwg_scripts/mysqlimportl.csh -L $dbName  $tableFileName`;
if ($tableFileName ne $fileName) {
  `rm $tableFileName`;
}
print "done";
