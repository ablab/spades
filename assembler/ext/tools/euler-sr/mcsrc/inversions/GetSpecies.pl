#!/usr/bin/env perl


$pattern = shift @ARGV;
@sf = glob($pattern);

@s = ();
foreach $file (@sf) {
  $file =~ /(\w+[^\.])\..*/;
  $sp = $1;
  push @s, $sp;
}


@s = sort @s;

print "@s\n";
