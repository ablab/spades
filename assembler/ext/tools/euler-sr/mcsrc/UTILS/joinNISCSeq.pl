#!/usr/bin/env perl

use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlast;

if (@ARGV < 2) {
  print "usage: joinNISCSeq giFile species_to_join \n";
  exit(1);
}


$giFile = shift @ARGV;
$species = shift @ARGV;

open (GI, "$giFile") or die "cannot open $giFile\n";

@fn = ();
@spec = ();
@acc = ();
@clone = ();
while (<GI>) {
  $line = $_;
  if ($line !~ /Redundant/ and $line !~ /Mapped/ and $line =~ /$species/) {
    $line =~ /([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)/;
    push @clones, $1;
    push @species, $2;
    push @accessions, $4;
    push @fileNames, "$species.$4.gbk";
    print "got values: $1 $2 $4 \n";
  }
}
close GI;

# Configure blast2seq
$bl2seqOut = "bl2seq.out";
$home = $ENV{'HOME'};
$BLASTDIR = $home . "/software/blast-2.2.13/bin/";
$blastSA = Bio::Tools::Run::StandAloneBlast->new('program'=>'blastn',
						 'outfile'=>$bl2seqOut,
						 'm'=>'T');

#$blastSA->PROGRAMDIR($BLASTDIR);
$blastSA->program_dir($BLASTDIR);
$path = $blastSA->program_dir;

print "bl2seqdir: $path\n"; 
#read in sequences
%sequences = ();
%ovpCoords = ();
%ovpAccns = ();
foreach $file (@fileNames) {
  if (-e $file) {
    $seqio_obj = Bio::SeqIO->new(-file => $file,
			       -format => "genbank" );

    $seq_obj = $seqio_obj->next_seq;
    
    $acc = $seq_obj->accession;
    $l = $seq_obj->length;
    $sequences{$acc} = $seq_obj;
    @{$refOvpCoords{$acc}} = ();
    @{$ovpCoords{$acc}} = ();
    @{$ovpAccn{$acc}} = ();
    for my $feat_object ($seq_obj->get_SeqFeatures) {
      $primaryTag = $feat_object->primary_tag;
      $refOvpStart = $feat_object->start;
      $refOvpEnd   = $feat_object->end;
      print "got primary tag: $primaryTag $refOvpStart $refOvpEnd\n";
      if ($primaryTag eq "misc_feature" and
	  $feat_object->has_tag("note")) {
	
	@notes = ();
	push @notes , $feat_object->get_tag_values("note");
	$value = @notes[0]; # don't get more than one note
	
	if ($value =~ /overlaps/) {
	  $value =~ /.*(AC\d+).*/;
	  $accnum = $1;
	  $start = -1; $end = -1;
	  if ($value =~/.*represents overlap with nucleotides\s*(\d+)\-(\d+).*/) {
	    $start = $1; $end = $2;
	  }
	  else {
	    if ($value =~/.*\(nucleotides\s*(\d+)\-(\d+)\).*/) {
	      $start = $1; $end = $2;
	    }
	  }
	  push @{$ovpAccn{$acc}}, $accnum;
	  push @{$refOvpCoords{$acc}}, $refOvpStart;
	  push @{$refOvpCoords{$acc}}, $refOvpEnd;
	  push @{$ovpCoords{$acc}}, $start;
	  push @{$ovpCoords{$acc}}, $end;
	  print "$acc stored $accnum, $refOvpStart, $refOvpEnd, $start, $end\n";
	}
      }
    }
  }
}


%cloneNextOverlap = ();
%cloneNextOverlapRefStart = ();
%cloneNextOverlapRefEnd = ();
%cloneNextOverlapQryStart = ();
%cloneNextOverlapQryEnd = ();

# initialize the ordering
for $accn (@accessions) {
  $cloneNextOverlap{$accn} = "";         
  $cloneNextOverlapRefStart{$accn} = -1;
  $cloneNextOverlapRefEnd{$accn} = -1;
  $cloneNextOverlapQryStart{$accn} = -1;
  $cloneNextOverlapQryEnd{$accn} = -1;
}

# determine the ordering of the clones
for $accn (@accessions) {

  print "starting contig with $accn @clones[$accIndex]\n";

  # don't process this clone if it's already been added to a contig
  
  for $i (0 .. $#{$ovpAccn{$accn}}) {
    # store coordinate information for overlaps
    $ovpAccn = @{$ovpAccn{$accn}}[$i];
    $refOvpStart = @{$refOvpCoords{$accn}}[2*$i];
    $refOvpEnd = @{$refOvpCoords{$accn}}[(2*$i)+1];
    $ovpStart = @{$ovpCoords{$accn}}[2*$i];
    $ovpEnd = @{$ovpCoords{$accn}}[2*$i+1];
    $refOvpLen = $refOvpEnd - $refOvpStart + 1;
    $refLen = length $sequences{$accn}->seq;
    $ovpLen = $ovpEnd - $ovpStart + 1;
#    if ($ovpStart == -1 ) {
      # no overlap text found, do it manually.
      print "finding overlap for $accn $ovpAccn\n";
    if (exists $sequences{$accn} and
	exists $sequences{$ovpAccn}) {
      ($refOvpStart, $refOvpEnd, $ovpStart, $ovpEnd) = 
	FindOverlap(\$sequences{$accn}, \$sequences{$ovpAccn});
      print "manually found overlap: $refOvpStart, $refOvpEnd, $ovpStart, $ovpEnd\n";
    }
    #    }

    print "$i $accn $refOvpStart .. $refOvpEnd ($refOvpLen) overlap: $ovpAcc $ovpStart .. $ovpEnd ($ovpLen) \n";

    if ($refOvpStart == 1) {
      # this overlap defines an overlap of a preceeding clone
      if ($cloneNextOverlap{$ovpAccn} eq "") {
	# no overlap has been defined for the preceeding clone
	# define it here
	$cloneNextOverlap{$ovpAccn} = $accn;
	$cloneNextOverlapRefStart{$ovpAccn} = $ovpStart;
	$cloneNextOverlapRefEnd{$ovpAccn} = $ovpEnd;
	$cloneNextOverlapQryStart{$ovpAccn} = $refOvpStart;
	$cloneNextOverlapQryEnd{$ovpAccn} = $refOvpEnd;
      }
    }
    else {
      # this overlap defines an overlap of this contig and the next
      # one
      $cloneNextOverlap{$accn} = $ovpAccn;
      $cloneNextOverlapRefStart{$accn} = $refOvpStart;
      $cloneNextOverlapRefEnd{$accn} = $refOvpEnd;
      $cloneNextOverlapQryStart{$accn} = $ovpStart;
      $cloneNextOverlapQryEnd{$accn} = $ovpEnd;
    }
  }
}

# print a summary of each overlap

for $accn (@accessions) {
  print "$accn overlaps with $cloneNextOverlap{$accn}, ($cloneNextOverlapRefStart{$accn}- ";
  print "$cloneNextOverlapRefEnd{$accn}) ($cloneNextOverlapQryStart{$accn}-";
  print " $cloneNextOverlapQryEnd{$accn}) \n";
}


$accIndex = 0;
@contigs = ();
$fileName = "$species.fasta";
$seqout = Bio::SeqIO->new(-file=>">$fileName",
			  -format=>"fasta");

$contig = "";
$contigNumber = 0;
for $accnIndex (0 .. $#accessions) {
  # get the sequene for this contig
  $accn = @accessions[$accnIndex];
  if (exists $sequences{$accn}) {
    if ($contig eq "") {
      $contig = $sequences{$accn}->seq;
    }
    $joinedContig = 0;
    if ($accnIndex < (@accessions-1)) {
      # this is not the last contig. try to append it
      if ($cloneNextOverlap{$accn} ne "" ) {
	$overlapAccn = $cloneNextOverlap{$accn};
	if (exists $sequences{$overlapAccn} ) {
	  $overlapLength = $cloneNextOverlapRefEnd{$accn} - $cloneNextOverlapRefStart{$accn} + 1;
	  $ls = length $sequences{$accn}->seq;
	  $lc = length $contig;
	  print "joining contigs: $accn $ls $overlapAccn\n";
	  $sa = substr($contig, 0, -$overlapLength);
	  $sb = $sequences{$overlapAccn}->seq;
	  $lsa = length $sa;
	  $lsb = length $sb;
	  #      print "lsa: $lsa, $lsb\n";
	  $contig = $sa . $sb;
	  $lc = length $contig;
	  #      print "lc: $lc\n";
	  $joinedContig = 1;
	}
      }
    }
    if ($accnIndex == $#acessions or 
	! $joinedContig) {
      # this contig isn't aligned with anything.  
      # output it's contig.
      if (-e $fileName) {
	`ls  -l $fileName`;
      }
      $bpseq = Bio::Seq->new();
      $bpseq->seq($contig);
      $bpseq->display_id($species."_".$contigNumber);
      print "outputting seq of len ", $bpseq->length, "\n";
      $seqout->write_seq($bpseq);
      $contig = "";
      ++$contigNumber;
    }
  }
}


sub FindOverlap {
  my ($s1, $s2) = @_;
  print "got values: $$s1, $$s2\n";
  $bl2seqReport = $blastSA->bl2seq($$s1, $$s2);
  $refOvpStart = -1;
  $refOvpEnd   = -1;
  $qryOvpStart = -1;
  $qryOvpStart = -1;
  # get the highest scoring hit
  if ( $result = $bl2seqReport->next_result ) {
    if ($hit = $result->next_hit ) {
      if ($hsp = $hit->next_hsp ) {
	($rBegin, $rEnd) = $hsp->range('query');
	($qBegin, $qEnd) = $hsp->range('sbjct');
	print "got blast result: $rBegin, $rEnd, $qBegin, $qEnd\n";
	# determine roughly the overlapping side
	$s1Length = $$s1->length;
	$s2Length = $$s2->length;
	if ($rBegin > $qBegin) {
	  # s2 overlaps the end of s1
	  $refOvpStart = $rBegin - $qBegin;
	  $refOvpEnd   = $s1Length;
	  $qryOvpStart = 1;
	  $qryOvpEnd   = $qEnd;
	}
	else {
	  # s2 overlaps the beginning of s1
	  $refOvpStart = 1;
	  $refOvpEnd   = $rEnd;
	  $qryOvpStart = $qBegin - $rBegin;
	  $qryOvpEnd   = $s2Length;
	}
      }
    }
  }
  return ($refOvpStart, $refOvpEnd, $qryOvpStart, $qryOvpEnd);
}




