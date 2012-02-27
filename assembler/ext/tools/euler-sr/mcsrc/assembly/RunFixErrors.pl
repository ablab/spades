#!/usr/bin/env perl

$readsFile    = shift @ARGV;
$spectrumFile = shift @ARGV;
$outputFile = shift @ARGV;

`fixErrorsSAP $readsFile $spectrumFile $outputFile -minMult 5 -maxGap 4 -discardFile $readsFile.discards -maxScore 10 -maxTrim 5`;

