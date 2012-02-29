#!/usr/bin/perl -w 

use strict;
use NaifDB;

#
# csv files should look like this
#
# key1,  key2,  ..., keyn
# val11, val12, ..., val1n
# val21, val22, ..., val2n
# ...
# valN1, valN2, ..., valNn
#


my $db = new_naif_db()->from_csv("test.csv");

printf "db size = %d\n", $db->size();

my $dbnull = new_naif_db()->from_csv("testnull.csv");


my $db123 = $db->record_grep("lib", "solexa123");

my $record = $db->record_get(-1);


my $db0 = new_naif_db();
$db0->to_csv("trete.csv");

$dbnull->record_push_back($record);

$dbnull->to_csv("testnull.csv");

$db123->record_set(0, $record);

$db123->to_csv("testout.csv");


