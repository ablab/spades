#!/usr/bin/env perl

$id1 = shift @ARGV;
$id2 = shift @ARGV;

$grimmSyntenyDir = $ENV{"GSDIR"};

$logfile = "find_inversions.txt";
while ($#ARGV >= 0) {
  $option = shift @ARGV;
  if ($option eq "-gsdir") {
    $grimmSyntenyDir = shift @ARGV;
  }
  if ($option eq "-grimmdir") {
    $grimmDir = shift @ARGV;
  }
}

# fragment the alignment according to overlap
`mkdir -p fragmented`;
$mach = $ENV{"MACHTYPE"};
$srcdir = $ENV{"EUSRC"};
Exec("$srcdir/lav/$mach/fragmenter lav/$id1.$id2.lav fragmented/$id1.$id2.lav", $logfile);

# rescore the alignments, since sometimes the blastz ones don't
# make any sense
Exec("$srcdir/lav/$mach/rescore fragmented/$id1.$id2.lav $id1.fasta $id2.fasta fragmented/$id1.$id2.rescored.lav", $logfile);

# select the highest scoring framgents
Exec("mkdir -p selected", $logfile);
Exec("$srcdir/lav/$mach/blockselect fragmented/$id1.$id2.rescored.lav $id1.fasta $id2.fasta selected/$id1.$id2.lav", $logfile);

# turn the output into a grimm-ok format
Exec("mkdir -p grimm/$id1.$id2", $logfile);
Exec("$srcdir/lav/$mach/lav2grimm selected/$id1.$id2.lav grimm/$id1.$id2/coords.txt", $logfile);

# run grimm-synteny Anchors on the output
Exec("$grimmSyntenyDir/grimm_synt -f grimm/$id1.$id2/coords.txt -d grimm/$id1.$id2 -A", $logfile);

# run grimm-synteny on the output
Exec("$grimmSyntenyDir/grimm_synt -f grimm/$id1.$id2/unique_coords.txt -d grimm/$id1.$id2 -g 10000 -m 2000 -c", $logfile);

Exec("$srcdir/InvFinder/ParseMgrMicro.pl grimm/$id1.$id2/mgr_micro.txt", $logfile);

sub Exec {
  my ($cmd, $logfile) = @_;
  $date = `date`;
  chomp $date;
  system("$cmd");
  $status = $?;
  `echo "$date\t$cmd\tSTATUS: $status" >> $logfile`;
  if ($status != 0) {
    exit(0);
  }
}
