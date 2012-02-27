#!/usr/bin/env perl
use strict;

if ($#ARGV != 1) {
  print "usage: $0 source_dir new_dir \n";
  exit(-1);
}
my($sourceDir, $destDir) = @ARGV;
my @dbfiles = glob("$sourceDir/*.db.*");
my @basefiles = glob("$sourceDir/*.corr");
push @basefiles, glob("$sourceDir/*.mates.10");
push @basefiles, glob("$sourceDir/name.rul");

my @etfiles = ();
my @newbasefiles = ();
grep( (/.*\/(.*)db(.*)/ && push(@etfiles, $1 . "et" . $2)), @dbfiles);

grep ( (/.*\/(.*)/ && push(@newbasefiles, $1)), @basefiles);

print "dbfiles: @dbfiles \n";
print "etfiles: @etfiles \n";

print "new base files: @newbasefiles \n";

`pwd`;
`mkdir $destDir`;

my $i;
foreach $i (0 .. $#dbfiles) {
  `cd $destDir; ln -s ../$sourceDir/@dbfiles[$i] @etfiles[$i]  `;
}

foreach $i (0 ..$#basefiles) {
  `cd $destDir; ln -s ../$sourceDir/@basefiles[$i] @newbasefiles[$i]`;
}
# 
