#!/usr/bin/env perl

use ReadLibrary;
use findReadInfo;
use ABI;

if ($#ARGV != 1) {
  print "usar: $0 reads_file.fasta xml\n";
  exit(0);
}


$in = shift @ARGV;
$out = shift @ARGV;


open (FASTAIN, $in) or die "no puedo abrir $in\n";
open (XMLOUT, ">$out") or die "cannot open $out\n";
# Print xml header information


$numParsed = 0;
$numFailed = 0;
%templates = ();
while (<FASTAIN>) {
  $line = $_;
  chomp $line;
  if ($line =~ /^>(.*)/) {
    $title = $1;
    ($baseName, $dir, $type, $fullName) = ReadLibrary::ParseABITitle($title);
#    ($organismo, $largoClon, $read, $clon, $ext) = &findReadInfo($title);
    $numParsed++;
    $template = $baseName;
    if (!exists $templates{$template}) {
      $templates{$template} = 1;
    } else {
      $templates{$template}++;
    }
  }
}

print XMLOUT "<rule>\n";
print XMLOUT "  <name> all reads are paired production reads </name>\n";
print XMLOUT "  <match>\n";
print XMLOUT "    <match_field>trace_name</match_field>\n";
print XMLOUT "    <regex>.</regex>\n";
print XMLOUT "  </match>\n";
print XMLOUT "  <action>\n";
print XMLOUT "    <set>\n";
print XMLOUT "      <set_field>type</set_field>\n";
print XMLOUT "      <value>paired_production</value>\n";
print XMLOUT "    </set>\n";
print XMLOUT "  </action>\n";
print XMLOUT "</rule>\n";

@tk = keys %templates;
$ntk = $#tk;
foreach $template (keys %templates) {
  if ($templates{$template} != 1) {
print XMLOUT "<rule>\n";
print XMLOUT "   <name> all reads are paired production reads </name>\n";
print XMLOUT "   <match>\n";
print XMLOUT "      <match_field>template_id</match_field>\n";
print XMLOUT "      <regex>$template</regex>\n";
print XMLOUT "   </match>\n";
print XMLOUT "   <action>\n";
print XMLOUT "      <set>\n";
print XMLOUT "         <set_field>type</set_field>\n";
print XMLOUT "         <value>unpaired_production</value>\n";
print XMLOUT "      </set>\n";
print XMLOUT "   </action>\n";
print XMLOUT "</rule>\n";
  }
}

close XMLOUT;
