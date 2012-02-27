#!/usr/bin/env perl
use Bio::SeqIO;

$in = shift @ARGV;
$name = shift @ARGV;

$seqIn = new Bio::SeqIO('-file'=>$in);
if ($first = $seqIn->next_seq()) {
$s = $first->seq;
print "$s\n";
}
