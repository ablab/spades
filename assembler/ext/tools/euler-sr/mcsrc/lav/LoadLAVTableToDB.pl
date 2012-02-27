#!/usr/bin/env perl

use DBI;

$dbName = shift @ARGV;
$defaultsFile = shift @ARGV;
$socketFile   = shift @ARGV;

%attr = {defaults-file=>$defaultsFile,
	 socket=>$socketFile};

$dbh = DBI->connect("dbi:mysql:$dbName", "mchaisso", "", \%attr);


  columns.push_back("id"); columns.push_back(" int AUTO_INCREMENT PRIMARY KEY ");
  columns.push_back("tStart"); columns.push_back("int");
  columns.push_back("tEnd");   columns.push_back("int");
  columns.push_back("qStart"); columns.push_back("int");
  columns.push_back("qEnd");   columns.push_back("int");
  columns.push_back("strand"), columns.push_back("int");
  columns.push_back("chainId"); columns.push_back("int");
  columns.push_back("tid"); columns.push_back("int");
  columns.push_back("qid"); columns.push_back("int");
  indices.push_back("tStart");
  indices.push_back("tEnd");
  indices.push_back("qStart");
  indices.push_back("qEnd");
  indices.push_back("chainId");

my $sth = $dbh->prepare($queryStr);
$sth->execute();


