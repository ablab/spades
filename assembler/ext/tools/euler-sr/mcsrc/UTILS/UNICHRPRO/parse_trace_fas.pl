#!/usr/bin/env perl

$traceFasFileName = shift @ARGV;

open(TFF, "$traceFasFileName") or die "cannot open $traceFasFileName\n";

@traceFasFile = <TFF>;

@titles = grep(/^>.*/ , @traceFasFile);

%names = {};
%mates = {};
%mate_names = {};
%directions = {};
foreach $title (@titles) {
  if ($title =~ /^>gnl\|ti\|(\S*) name:(\S+) mate:(\S+) mate_name:(\S+) .* end:([FR])/) {
    print "got $1 $2 $3 $4 $5\n";
    $names{$1} = $2;
    $mates{$1} = $3;
    $mate_names{$1} = $4;
    $directions{$1} = $5;
  }
}
    

