#!/usr/bin/env perl
die "usage: $0 base_filename probe_len  mismatch_count outdir\n" unless ($#ARGV==3);

$base = @ARGV[0];
$k    = @ARGV[1];
$d    = @ARGV[2];
$outdir = @ARGV[3];

@files = glob("*.fa");
%numbers = ();
grep /\w([1-9XY][0-9]*)\./ && ($numbers{$1} = 1), @files;

@chrs = keys(%numbers);
$numchrs = $#chrs;
$i;
foreach $chr  (@chrs) {
  $chrstr = "$base\{";
  $i = 0;
  @subsets = glob("*chr$chr.*");
  foreach $chr2 (@chrs) {
    if ($chr2 ne $chr) {
      if ($i != ($numchrs-1)) {
	 $chrstr = $chrstr . "$chr2.*,";
	 $i = $i + 1;
      }
      else {
	 $chrstr = $chrstr . "$chr2.*}.fa";
      }
    }
  }
  foreach $file (@subsets) {
    print "cd /home/mchaisso/GENOMES/MM/SPLIT_FILES; nohup /home/mchaisso/PROBE_LOCATE-intel/unichr_select -d $d -k $k -u -o $outdir/$file.probes $file $chrstr > $outdir/$file.out &\n ";
  }
}

  

