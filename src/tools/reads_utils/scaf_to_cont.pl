#!/usr/bin/perl
use strict;

my $file = $ARGV[0];
open(IN,$file) || die "Incorrect file $file. Exiting...\n";
#my $out_file = $ARGV[1];
#sysopen(OUT,$out_file,O_WRONLY|O_CREAT) || die "Incorrect file $out_file. Exiting...\n";

my ($seq, $name)=('','');
while(<IN>){
  chomp;
  my $line = $_;
  $seq.= uc($line) if(eof(IN));
  if (/\>(\S+)/ || eof(IN)){
    if($seq ne ''){
      my @seqgaps = split(/[N]{1,}/, $seq);
      if($#seqgaps > 0){
        my $ctgcount=0;
        foreach my $ctgseq (@seqgaps){
          $ctgcount++;
          print "$name contig$ctgcount (size=".length($ctgseq).")\n$ctgseq\n";
        }
      }else{
        print ">$name\n$seq\n";
      }
    }
    $seq='';
    $name = $_;
  }else{
    $seq.=uc($line);
  }
}
