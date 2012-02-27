#!/usr/bin/env perl

use DBI;
if ($#ARGV < 5) { print "$#ARGV usage: $0 dbname defaultsFile socketFile refSeq qrySeq seqame \n"; exit(0);}
$scope = "local";
$dbName       = shift @ARGV;
$defaultsFile = shift @ARGV;
$socketFile   = shift @ARGV;
$refName       = shift @ARGV;
$qryName       = shift @ARGV;
$seqName      = shift @ARGV;
$fileName     = shift @ARGV;
$scope   = shift @ARGV;

print "$0 dbName: $dbName\n";
print "$0 defaultsFiel: $defaultsFile\n";
print "$0 socketFile $socketFile\n";
print "$0 refName $refName\n";
print "$0 qryName $qryName\n";
print "$0 seqName $seqName\n";
print "$0 fileName $fileName\n";
                                                                                                                                             
# fix the queryid in the table
#print "rs: $refSeq  qs: $qrySeq\n";

$tableName = "$refName\_$qryName\_$seqName";
$tableName =~ s/\./_/g;
$tableFileName = "$tableName.txt";

#print "using socket: $socketFile\n";
$dbh = DBI->connect("dbi:mysql:$dbName;socket=$socketFile;mysql_read_default_file=$defaultsFile", 
		    "mchaisso", "");

$qryStr = "DROP TABLE IF EXISTS $tableName  ;";
#print "query: $qryStr\n";
$stmt = $dbh->prepare($qryStr); $stmt->execute();

print "creating $tableName \n";
$qryStr = "CREATE TABLE $tableName (id INT AUTO_INCREMENT PRIMARY KEY, " .
  "tStart INT, " .
  "tEnd INT, " .
  "qStart INT, " .
  "qEnd INT, " .
  "strand INT, " .
  "chainId INT, " .
  "tid INT, " .
  "qid INT )";

$stmt = $dbh->prepare($qryStr); $stmt->execute();

$qryStr = "ALTER TABLE $tableName add index (tStart)";
$stmt = $dbh->prepare($qryStr); $stmt->execute();

$qryStr = "ALTER TABLE $tableName add index (tEnd)";
$stmt = $dbh->prepare($qryStr); $stmt->execute();

$qryStr = "ALTER TABLE $tableName add index (qStart)";

#print "using socket: $socketFile\n";
$dbh = DBI->connect("dbi:mysql:$dbName;socket=$socketFile;mysql_read_default_file=$defaultsFile", 
		    "mchaisso", "");

$qryStr = "DROP TABLE IF EXISTS $tableName  ;";
#print "query: $qryStr\n";
$stmt = $dbh->prepare($qryStr); $stmt->execute();

print "creating $tableName \n";
$qryStr = "CREATE TABLE $tableName (id INT AUTO_INCREMENT PRIMARY KEY, " .
  "tStart INT, " .
  "tEnd INT, " .
  "qStart INT, " .
  "qEnd INT, " .
  "strand INT, " .
  "chainId INT, " .
  "tid INT, " .
  "qid INT )";

$stmt = $dbh->prepare($qryStr); $stmt->execute();

$qryStr = "ALTER TABLE $tableName add index (tStart)";
$stmt = $dbh->prepare($qryStr); $stmt->execute();

$qryStr = "ALTER TABLE $tableName add index (tEnd)";
$stmt = $dbh->prepare($qryStr); $stmt->execute();

$qryStr = "ALTER TABLE $tableName add index (qStart)";
$stmt = $dbh->prepare($qryStr); $stmt->execute();

$qryStr = "ALTER TABLE $tableName add index (qEnd)";
$stmt = $dbh->prepare($qryStr); $stmt->execute();

$qryStr = "ALTER TABLE $tableName add index (chainId)";
$stmt = $dbh->prepare($qryStr); $stmt->execute();

#print "importing\n";
#`mv $tmpFileName $tableName.txt`;
print "importing $tableFileName\n";
if (-e $fileName and $fileName ne $tableFileName) {
  system("cp $fileName $tableFileName");
}
$importCommand = "~/fwg_scripts/mysqlimportl.csh";
if ($scope eq "global") {
  $importCommand = "~/fwg_scripts/mysqlimportg.csh";
}
elsif ($scope eq "local") {
  $importCommand = "~/fwg_scripts/mysqlimportl.csh";
}

$cmd  = "$importCommand -L --socket=$socketFile $dbName $tableFileName";
print "running $cmd\n";
system($cmd);
if ($fileName ne $tableFileName) {
#  system("rm $tableFileName");
}
