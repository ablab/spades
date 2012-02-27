#!/usr/bin/env perl

if ($#ARGV < 4) {
  print "usage: $0 filename dbname refname qryname seqname [{global|local}]\n";
  exit(0);
}

use DBI;

$fileName = shift @ARGV;
$dbName   = shift @ARGV;
$refName   = shift @ARGV;
$qryName  = shift @ARGV;
$seqName  = shift @ARGV;
$scope = "local";
if ($#ARGV >= 0) {
  $scope = shift @ARGV;
}
print "option fileName $fileName\n";
print "option dbName $dbName\n";
print "option refName $refName\n";
print "option qryName $qryName\n";
print "option seqName $seqName\n";
print "option scope   $scope\n";
# fix the queryid in the table

$tableName = "$refName\_$qryName\_$seqName\_chain";
$tableFileName = $tableName . ".txt";

$dbConnect = "dbi:mysql:$dbName";
$importcommand = "";
if ($scope eq "global") {
  if (exists $ENV{"GSCK"}) {
    $socket = $ENV{"GSCK"};
    print "using global socket: $socket\n";
    $dbConnect .= ";mysql_socket=$socket";
 }
 $importcommand = "~/fwg_scripts/mysqlimportg.csh";
}
else {
  if (exists $ENV{"LSCK"}) {
    $socket = $ENV{"LSCK"};
    $dbConnect .= ";mysql_socket=$socket";
  }
  $importcommand = "~/fwg_scripts/mysqlimportl.csh";
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

print "importing $fileName\n";
`ln $fileName $tableFileName `;
$cmd= "$importcommand -L $dbName $tableFileName";
print "running $cmd\n";
system($cmd);
print "done";
