#!/usr/bin/env perl

use ReadLibrary;

if ($#ARGV != 1) {
  print "\n usar: $0 reads_file.fasta xml\n\n";
  exit(0);
}


$in = shift @ARGV;
$xmlOut = shift @ARGV;

open (FASTAIN, $in) or die "no puedo abrir $in\n";
open (XML, ">$xmlOut") or die "no puedo escribir $xmlOut\n";
# Print xml header information
print XML "<?xml version=\"1.0\"?>\n";
print XML "<TRACE_VOLUME>\n";

$numParsed = 0;
$numFailed = 0;
%libraries = ();
%templates = ();
while (<FASTAIN>) {
  $line = $_;
  chomp $line;
  if ($line =~ />/) {
    ($template, $dir, $type, $name) = ReadLibrary::ParseABITitle($line);
    $largoClon = $ReadLibrary::MATE_MEAN[$type];
    $stddevClon = $ReadLibrary::MATE_STDEV[$type];

    $numParsed++;
    print XML "\t<TRACE>\n";
    print XML "\t\t<TRACE_NAME>$name</TRACE_NAME>\n";
    print XML "\t\t<PLATE_ID>$template</PLATE_ID>\n";
    print XML "\t\t<WELL_ID>A1</WELL_ID>\n";
    print XML "\t\t<TEMPLATE_ID>$template</TEMPLATE_ID>\n";
    print XML "\t\t<TRACE_END>$dir</TRACE_END>\n";
    print XML "\t\t<INSERT_SIZE>$largoClon</INSERT_SIZE>\n";
    print XML "\t\t<INSERT_STDEV>$stddevClon</INSERT_STDEV>\n";
    print XML "\t\t<LIBRARY_ID>library_$largoClon</LIBRARY_ID>\n";
    print XML "\t</TRACE>\n";
  }
}
print XML "</TRACE_VOLUME>\n";

close XML;
