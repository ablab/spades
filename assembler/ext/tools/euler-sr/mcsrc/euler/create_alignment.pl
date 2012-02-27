#!/usr/bin/perl -w

#***************************************************************************
# Title:          create_alignment.pl
# Author:         Vagisha Sharma
# Created:        Jun. 2002
# Last modified:  May. 2004
#
# Copyright (c) 2001-2004 The Regents of the University of California
# All Rights Reserved
# See file LICENSE for details.
#**************************************************************************/

use strict;
use Data::Dumper;

my $intvfile = shift;
my $readsfile = shift;
my $line_w = shift;
my $separation = shift;
my $mafile = shift;



$line_w = 80  if (!defined $line_w);
$separation = 3 if (!defined $separation);

if ((!-e $intvfile) || (!-e $readsfile)) {
  print "either $intvfile or $readsfile does not exist\n";
  exit;
}


# line separating blocks and contigs
my $spacer1 = "";
for (my $i = 0; $i < $line_w; $i++) {
$spacer1 .= "-";
}
my $spacer2 = "";
for (my $i = 0; $i <= $separation; $i++) {
$spacer2 .= "\n";
}

my $blankline = "";
for (my $i = 0; $i < $line_w; $i++) {
$blankline .= " ";
}
for (my $i = 0; $i < $separation; $i++) {
#$blankline .= "\n";
}


my $error_out = "";
#my $ma_file = "RAconverter.input";

my $allcontigs = {};
# read the ace file
read_acefile();
# read the sequence file
my $total_reads = read_readseqs();

#$total_reads--; # reads indices start at 0;
$error_out .= "total number of reads in the reads file: $total_reads\n";

order_reads();

print "total reads : $total_reads\n";
my $output = print_alignment();
open (OUT, ">$mafile") or die "could not open $mafile for writing\n";
print OUT $output;
#print $output;
close OUT;
########### error output
print $error_out;


exit;


##############################################################################
# read in the ace file and srote the read coordinates on each contig
##############################################################################
sub read_acefile {
	my $currcontig;
	my $maxreadnum = 0;
	open (INTV, "$intvfile") or die "could not open $intvfile for reading\n";
	## MODIFIED Oct 20
	my $maxconlen = 0;
	my $prevconnum = -1;

	## ADDED Oct 26
	my $readnum = 0;

	while (<INTV>) {
	  my $line = $_;
	  #print "$line\n";
	  chomp($line);
	  ## ADDED 'EDGE' on Oct 31 for running with chimp intv file
	  if ($line =~ m/Contig/ || $line =~ m/EDGE/) {
	    #print "$line\n";
	    
	    ## MODIFIED Oct 20
	    if ($prevconnum != -1) {
	      my $contig = $allcontigs->{"$prevconnum"};
	      $contig->{'conlen'} = $maxconlen;
	      $maxconlen = 0;
	    }


	    my @tokens = split (/\s+/, $line);
	    my $contignum = $tokens[1];
	    my $conlen = $tokens[3];
	    my $contig = {};
	    $contignum--; # numbers in ace file are incremented by 1
	    $prevconnum = $contignum; ## MODIFIED Oct 20
	    $allcontigs->{"$contignum"} = $contig;
	    # MODIFIED Oct 20 conlen should be determined from the read coordinates
	    #$contig->{'conlen'} = $conlen;
	    my $reads = {};
	    $contig->{'reads'} = $reads;
	    $currcontig = $contig;
	    #print Dumper($allcontigs->{"Contig$contignum"});
	
	    ## ADDED Oct 26
	    $readnum = 0 # reset the readnum
	  }
	
	  elsif ($line =~ m/^INTV/) {
	    #print "$line\n";
	    my @tokens = split(/\s+/, $line);
	    #my $readnum = $tokens[1];
	    my $intvrdnum = $tokens[1];
	    my $read_s = $tokens[2];
	    my $match_l = $tokens[3];
	    my $contig_s = $tokens[4];
	    my $read = {};
	    $read->{'read_s'} = $read_s;
	    $read->{'match_l'} = $match_l;
	    $read->{'contig_s'} = $contig_s;

	    ## added Oct 26 : Can't index reads by readnum in intv file because one read may appear > once
	    $read->{'intvrdnum'} = $intvrdnum;

	    if (defined $currcontig->{'reads'}->{"$readnum"}) {
	      print "$readnum found twice!!!\n";
	    }
	    $currcontig->{'reads'}->{"$readnum"} = $read;
	    #$maxreadnum = $readnum   if ($readnum > $maxreadnum);
	    $maxreadnum = $intvrdnum   if ($intvrdnum > $maxreadnum);

	    ## MODIFIED Oct 20:
	    my $conlen = $contig_s + $match_l - 1;
	    $maxconlen = $conlen    if ($maxconlen < $conlen);	
	
	    ## added Oct 26:
	    $readnum++;

	  }
	}
	
	close INTV;
	 ## MODIFIED Oct 20: set the maxlen for the last contig
	if ($prevconnum != -1) {
	  my $contig = $allcontigs->{"$prevconnum"};
	  $contig->{'conlen'} = $maxconlen;
	  $maxconlen = 0;
	}	
 	$error_out .= "maxreadnum in acefile: $maxreadnum\n";
	#print Dumper($allcontigs);
}
##############################################################################
# order reads such that each row has non overlapping reads
##############################################################################
sub order_reads {
  #my ($allcontigs) = @_;

  my @contigkeys = ();
  foreach my $key (keys(%$allcontigs)) {
    push (@contigkeys, $key)   if ($key =~ m/(\d+)/);
  }
  @contigkeys = sort {$a <=> $b} @contigkeys;


  foreach my $conkey (@contigkeys) {
    my $contig = $allcontigs->{"$conkey"};
    my $conlen = $contig->{'conlen'};
    my $reads = $contig->{'reads'};
    # sort the reads on contig start positions
    my @readkeys = sort {$reads->{$a}->{'contig_s'} <=> $reads->{$b}->{'contig_s'}} (keys(%$reads));
    my $rownum = 0;

    for (my $i = 0 ; $i < @readkeys; $i++) {
      my $readkey = $readkeys[$i];
      my $read = $reads->{"$readkey"};

      # go to next if this read is already in a row
      next   if (defined $read->{'rownum'});

      my $con_s = $read->{'contig_s'};
      my $len = $read->{'match_l'};
      my $con_e = $con_s+$len-1;

      $read->{'rownum'} = $rownum;

      my $curr_e = $con_e;
      #print "looking at $readkey $con_s $con_e $curr_e\n";
      for (my $j = $i+1; $j < @readkeys; $j++) {
	
	my $rk = $readkeys[$j];
	my $r = $reads->{"$rk"};

	# go to next if this row has already been assigned to a row
	next   if (defined $r->{'rownum'});

	my $cs = $r->{'contig_s'};
	
	# if this read overlaps the previous read
	#next if ($cs <= $curr_e+1);
	next if ($cs <= $curr_e+20);

	my $ce = $cs + $r->{'match_l'} -1;
	
	#print "\tread num: $rk $cs $ce $curr_e\n";
	$curr_e = $ce  if ($curr_e < $ce);
	$r->{'rownum'} = $rownum;
      }

      $rownum++;
    }# sort the reads according to row, then according to con_s
    my @sortedkeys = sort { $reads->{$a}->{'rownum'} <=> $reads->{$b}->{'rownum'} ||
			    $reads->{$a}->{'contig_s'} <=> $reads->{$b}->{'contig_s'}} @readkeys;
    foreach my $key (@sortedkeys) {
      my $read = $reads->{"$key"};
      my $con_s = $read->{'contig_s'};
      my $con_e = $con_s + $read->{'match_l'} -1;
      my $rownum = $read->{'rownum'};

      #print "$rownum: $key ($con_s - $con_e)\n"  if ($conkey eq '2');
    }
    #print "number of rows : $rownum\n";
    #last;
  }
}


##############################################################################
# read in the reads file and store the sequence of each read
##############################################################################
sub read_readseqs {
  #print "in readseqs...............\n";

	my $readseqs = {};
	$allcontigs->{'readseqs'} = $readseqs;
	  open (READS, "$readsfile") or die "could not open $readsfile for reading\n";
	  my $readnum = 0; ### what is the starting readnum????
	  my $readseq = "";
	  my $thisreadnum = -1;
	  my $found = 0;
	  while (<READS>) {
	    my $line = $_;
	    chomp ($line);
	    if ($line =~ m/^>/) {
	      #print "$readnum\n";
	      # print the sequence of the previously found read
	      if ($readseq ne "") {
		#print "$thisreadnum --> $readseq\n";
		$readseqs->{"$thisreadnum"} = $readseq;
	
	      }
	      $thisreadnum = $readnum;
	      $readnum++;
	      $readseq = "";
	    }
	    else {
	      $readseq .= $line;
	    }
	
	}
	close READS;
	# get the sequence for the last read in the file
	$readseqs->{"$thisreadnum"} = $readseq;
	
	return $readnum;
}

##############################################################################
# print the alignment of each contig in vertical (TODO) or horizontal format
##############################################################################
sub print_alignment {
  #print "in print_alignment\n";
  my $i = 0;
	my @contigkeys = ();
	foreach my $key (keys(%$allcontigs)) {
	  push (@contigkeys, $key)   if ($key =~ m/(\d+)/);
	}
	@contigkeys = sort {$a <=> $b} @contigkeys;
	my $output = "";
	foreach my $conkey (@contigkeys) {
	  my $contig = $allcontigs->{"$conkey"};
	  my $conlen = $contig->{'conlen'};
	  $output .= "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
	  $output .= "% Contig $conkey $conlen\n";
	  #print "contig number: $conkey...................\n";
	  $output .= print_horizalign($conkey);
	  $output .= $spacer1; ## midified Oct8
	  $output .= "\n";     ## modified Oct8
	  $i++;
 	  #last  if ($i == 1);
	}
	return $output;
}




##############################################################################
# print the alignment of a contig in horizontal format
##############################################################################
sub print_horizalign {
  my ($conkey) = @_;

  #print "in print_horizalign......\n";

  my $contig = $allcontigs->{"$conkey"};
  my $reads = $contig->{'reads'};

  #my @readkeys = sort {$reads->{$a}->{'contig_s'} <=> $reads->{$b}->{'contig_s'}} (keys(%$reads));
  my @readkeys = sort { $reads->{$a}->{'rownum'} <=> $reads->{$b}->{'rownum'} ||
			    $reads->{$a}->{'contig_s'} <=> $reads->{$b}->{'contig_s'}} (keys(%$reads));
#   my @readkeys = keys(%$reads);
  ## Modified Oct 8 : ltlimit is not always 0
  my $lt_limit = -1;
  #my $lt_limit = 0;
  #my $rt_limit = $line_w - $lt_limit -1;
  my $conlen = $contig->{'conlen'};
  my $output = "";

  
  my $num_conreads = @readkeys;
  $output.=  "% number of reads for contig: $num_conreads\n";
  $output .= "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";

  #my $buffer = "";
  my $maxrow = 0;
  my $rows = {};
  # split the reads into rows
  # Modified -- Oct 8 : also get the lt_limit: the leftmost starting point on the contig
  # this need not always be 0
  foreach my $readkey (@readkeys) {
    my $read = $reads->{"$readkey"};
    # get the lt_limit for this contig
    $lt_limit = $read->{'contig_s'}  if ($lt_limit == -1 || $read->{'contig_s'} < $lt_limit);
    my $readrow = $read->{'rownum'};
    if (!defined $rows->{"$readrow"}) {
      my $row = ();
      push (@$row, $readkey); 
      $rows->{"$readrow"} = $row;
    }
    else {
      my $row = $rows->{"$readrow"};
      push (@$row, $readkey);
      }
  }
  #print "start is $lt_limit\n";
  my $rt_limit = $line_w + $lt_limit -1;

  #print Dumper ($rows);
  # sort rows according to row number
  my @sortedrows = sort {$a <=> $b} (keys(%$rows));


  while ($lt_limit < $conlen) {
    my $row = 0;
    my $buffer = "";
    foreach my $rowkey (@sortedrows) {
      my $row = $rows->{"$rowkey"};
      #print "row num $rowkey\n";
      #print Dumper ($row);
      my $thisoutput = print_onerow($contig, $rows, $rowkey, $lt_limit, $rt_limit);
      if ($thisoutput eq "") {
	$buffer .= "$blankline\n";
	#$buffer .= "";
	
      }
      else {
	$output .= $buffer;
	$output .= $thisoutput;
	$output .= "\n";
	$buffer = "";
	}
      # end of the row; print newline
      #$output .= "\n";
      $row++;
    }

    $maxrow = $row  if ($row > $maxrow);

    $lt_limit = $rt_limit +1;
    $rt_limit += $line_w;
    #$output .= "$lt_limit\t";
    $output .= $spacer1;
    $output .= $spacer2;
    #$output .= "\t$rt_limit\n";
    #$output .= "----------------------------------------------------------------------\n";

     
  }
  #print "max row num: $maxrow\n";
  return $output;
}


##############################################################################
# print one row of the alignment  
##############################################################################
sub print_onerow {
  my ($contig, $rows, $rownum, $lt_limit, $rt_limit) = @_;

  #print "in print_onerow......\n";

  my $row = $rows->{"$rownum"};
  my $conreads = $contig->{'reads'};

  my $curr_e = $lt_limit;
  my $output = "";
  my $nums = "";
  foreach my $readnum (@$row) {
    my $read = $conreads->{"$readnum"};
    my $con_s = $read->{'contig_s'};
    my $con_e = $con_s + $read->{'match_l'} - 1;
    ## ADDED Nov3 
    my $intvrdnum = $read->{'intvrdnum'};

    next   if ($con_e < $lt_limit);
    last   if($con_s > $rt_limit);

    #print "printing one read\n";
    my $thisreadout = "";
    ($thisreadout, $curr_e) = print_read($contig, $readnum, $curr_e, $rt_limit);
    $output .= $thisreadout if ($thisreadout ne "");
    #$nums .= "$readnum  "   if ($thisreadout ne "");
    $nums .= "$intvrdnum  "   if ($thisreadout ne "");
    #last   if ($con_s > $rt_limit);
  }
  # print blanks at the end if required
  if ($curr_e < $rt_limit && $output ne "") {
    my $gaplen = $rt_limit - $curr_e +1;
    for (my $i = 0; $i < $gaplen; $i++) {
      $output .= " ";
    }
    #last   if ($con_s > $rt_limit);
  }
#  $output .= " $curr_e"   if ($output ne "");
#  $output .= " ($nums)"   if ($output ne "");
  #$output .= "\n"  if ($output ne "");;
  #print "one row output is $output\n";
  return $output;
}

##############################################################################
# print one row of the alignment  
##############################################################################
sub print_read {
  my ($contig, $readnum, $lt_limit, $rt_limit) = @_;

  #print "in print_read.........\n";
  my $allreadseqs = $allcontigs->{'readseqs'};
  
  my $line_w = $rt_limit - $lt_limit +1;
  my $reads = $contig->{'reads'};

  my $output = "";
  my $complement = 0;
  #my $complreadnum = $readnum; 

  #########################
  ## COMMENTED OUT Oct26
#  if ($readnum >= $total_reads) {
#    $complreadnum = $readnum - $total_reads;
#    $complement = 1;
#  }
  my $read = $reads->{"$readnum"};
 
  my $s = $read->{'read_s'};
  my $con_s = $read->{'contig_s'};
  my $len = $read->{'match_l'};
  my $con_e = $con_s+$len-1;
  my $seq = "";

  ## ADDED Oct 26
  my $intvrdnum = $read->{'intvrdnum'};
  my $complreadnum = $intvrdnum; 
  #print "INTV read num is $intvrdnum  for $readnum $total_reads\n";
  if ($intvrdnum >= $total_reads) {
    $complreadnum = $intvrdnum - $total_reads;
    $complement = 1;
  }

  $seq = $allreadseqs->{"$intvrdnum"}  if ($intvrdnum == $complreadnum);
  $seq = $allreadseqs->{"$complreadnum"}  if ($intvrdnum != $complreadnum);

  ##COMMENTED OUT Oct26
#  $seq = $allreadseqs->{"$readnum"}  if ($readnum == $complreadnum);
#  $seq = $allreadseqs->{"$complreadnum"}  if ($readnum != $complreadnum);
  
  if (!defined $seq ) {
	  $error_out .=  "readnum is $intvrdnum $readnum $complreadnum (sequence is not defined)\n";
    return $output;
  }
  
  my $readlen = length ($seq);

  # if complement read is required
  if ($complement == 1) {
	$seq = reverse ($seq);
    $seq =~ tr/AGCTagct/TCGAtcga/;
  }

  # reads are sorted according to start on contig
  # !!!! if this read has been printed before print blanks of width line_w
 # if ($con_e < $lt_limit) {
#   # my $gap = "";
##    for (my $i = 0; $i < $line_w; $i++) {
##      $gap .= " ";
##    }
##    $gap .= "\n";
#    $output = "blank";;
#    return $output;
#  }
  return ($output, $lt_limit)  if ($con_s > $rt_limit);
  return ($output, $lt_limit)  if ($con_e < $lt_limit);

  #$output .= "$lt_limit $rt_limit $con_s $con_e $s $len $readlen\n";
  

  my $lt_gap = 0;
  my $rt_gap = 0;
  my $l_space = "";
  my $r_space = "";
  my $start;
  my $offset;
  my $rowseq = "";

  # if con_s is within lt and rt limits
  if ($con_s >= $lt_limit && $con_s <= $rt_limit) {
    $lt_gap = $con_s - $lt_limit;
    $rt_gap = $rt_limit - $con_e;
    $start = $s;
    $offset = $con_e > $rt_limit ? $rt_limit-$con_s+1 : $con_e-$con_s+1;
    
  }
  # con_e is within the lt and rt limit
  elsif ($con_e >= $lt_limit && $con_e <= $rt_limit) {
    $lt_gap = $con_s - $lt_limit;
    $rt_gap = $rt_limit - $con_e;   
    #$start = $con_s < $lt_limit? $lt_limit-$con_s : $s;
    $start = $con_s < $lt_limit? $lt_limit-$con_s + $s: $s;
    #$offset = $len - $start;
    $offset = $con_e - $lt_limit + 1;
  }
  # lt and rt limit lie with the contig start and end of this read
  elsif ($con_s < $lt_limit && $con_e > $rt_limit) {
    my $loff = $lt_limit - $con_s;
    $start = $s + $loff;
    $offset = $line_w;
  }

  #$output .= "ltgap is $lt_gap\; rtgap is $rt_gapn";
  if ($lt_gap > 0) {
    for (my $i = 0; $i < $lt_gap; $i++) {
#       $l_space .= "x";
      $l_space .= " ";
    }
    $output .= "$l_space";
  }

  # TODO: sometimes the start coordinate is beyond the length of the read.. catch this error
  $rowseq = substr($seq, $start, $offset);
  if (!defined $rowseq) {
	  	#$error_out .=  "$seq\n";
	  	$error_out .= "Coordinates are not right for $readnum $complreadnum !!!!!!! $start $offset $len\n";
	  	return ("", $lt_limit);
  }

  # some coordinates are not right
  # check if substr returns error if the length argument to substr is more than the actual length of the string
  # ANS: no it does not; it just gets whateven it can
  my $rowseqlen = length ($rowseq);
  if ($rowseqlen < $offset) {
    my $padlen = $offset - $rowseqlen;
    for (my $i = 0; $i < $padlen; $i++) {
      $rowseq .= " ";
    }
  }

  $output .= $rowseq;

  ## rt_gap does not count if you do not add it to the curr_e
  ### CHECK THIS!! rt_gap is printed in print_onerow
  if ($rt_gap > 0) {
    for (my $i = 0; $i < $rt_gap; $i++) {
#       $r_space .= "+";
      $r_space .= " ";
    }
    #$output .= "$r_space";
  }
  #$output .= "\n";
  my $curr_e = $con_e > $rt_limit ? $rt_limit+1 : $con_e+1;
  #print "readseq is $output\n";
  return ($output, $curr_e);

}
   
