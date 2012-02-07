#!/usr/bin/perl
use warnings;
use strict;
use File::Basename;
use File::Spec;

#####################################################################
############################# PREP ##################################
#####################################################################
if(@ARGV != 6 and @ARGV != 4){ 
	my ($name, $path,$suff) = fileparse($0);
	my $usage  = "\n\tusage1: $name <reads1.fq> <reads2.fq> <contigs.fa> <out-dir> <min-len> <num-threads>\n";
	   $usage .= "\tOR\n";
	   $usage .= "\tusage2: $name <aln_pe.sam> <contigs.fa> <out-dir> <min-len>\n\n";
	   $usage .= "\tNote: for usage1 BWA (v0.6.*) executable must be placed in PATH\n\n";
	die $usage;
}
my $do_aln = 0;
if(@ARGV == 6){ $do_aln = 1; }

my $num_threads = 3;
my ($small_contigs, $no_reads_contigs, $logfile)  = ("contigsSmall", "contigsNoReads", "prep.log");
my ($tot_pe_reads, $tot_pe_aln, $flt_on_mult, $noX0) = (0, 0, 0, 0);
my %files; my %n_reads_c; my %contigs_h;
my ($tot_contigs, $no_reads, $too_small) = (0, 0, 0);

my ($reads1, $reads2, $contigs, $outdir, $sam, $min_len);
my ($bwa_aln1_cmd, $bwa_aln2_cmd);
my $bwa = "bwa";
if($do_aln){
	($reads1, $reads2, $contigs, $outdir,$min_len, $num_threads) = @ARGV;
	$reads1 = File::Spec->rel2abs($reads1);
    $reads2 = File::Spec->rel2abs($reads2);
	
	# BWA alignment commands (used to have: -e $max_ge, -o max_go, -L and -d 100)
	my ($max_ge, $max_go, $k, $n) = (100, 3, 3, 0.08);
	$bwa_aln1_cmd  = "$bwa aln $contigs $reads1 -t $num_threads -O 7 -E 2 -k $k -n $n -q 15";
	$bwa_aln2_cmd  = "$bwa aln $contigs $reads2 -t $num_threads -O 7 -E 2 -k $k -n $n -q 15";
}
else{ 
	($sam, $contigs, $outdir, $min_len) = @ARGV; 
	$sam = File::Spec->rel2abs($sam);
}

my $contigs_full = File::Spec->rel2abs($contigs);
my ($contigs_name, $path, $suf) = fileparse($contigs);

mkdir($outdir) or die "could not create $outdir, $!. Quitting...\n";
chdir($outdir) or die "could not change dir to $outdir, $!. Quitting...\n";
system("cp $contigs_full .");

#####################################################################
############################# MAIN ##################################
#####################################################################
write_log_startup();
if($do_aln){ 
	$sam = align_all_reads_to_all_contigs($reads1, $reads2, $contigs_name); 
}
sep_contigs($contigs_name);
filter_into_sam_files($sam);
close_files_write_stats();
clean_short_and_detect_empty();
write_summary();

#####################################################################
######################### SUBROUTINES ###############################
#####################################################################
sub write_log_startup{
	my $time = localtime;
	
	open(LOG, ">", $logfile) or die "could not open $logfile, $!. Quitting...\n";

	print LOG "\n";
    print LOG "$time --- Starting prep for NGS-Refine...\n\n";

	if($do_aln){
		print LOG "\tWill align reads:\n";
		print LOG 	"\t\t$reads1\n";
		print LOG 	"\t\t$reads2\n";
		print LOG "\tTo contigs:\n";
		print LOG 	"\t\t$contigs_full\n";
		print LOG "\n";
		print LOG "\tBWA parameters:\n";
		print LOG "\t\tbwa-aln-1: $bwa_aln1_cmd\n";
		print LOG "\t\tbwa-aln-2: $bwa_aln2_cmd\n";
		print LOG "\n";
	}

	print LOG "Minimum contig length: $min_len\n";
	close LOG;

	print "\n$time --- Starting NGS-Refine prep... \n\n";
}

#####################################################################
sub clean_short_and_detect_empty{
	# assumes contig and sam file names are of the form "ID.fa"
	mkdir($small_contigs) or die "could not mkdir $small_contigs, $!. Quitting...\n";
	mkdir($no_reads_contigs) or die "could not mkdir $no_reads_contigs, $!. Quitting...\n";
	my @dir_files = <*>;

	foreach my $file (@dir_files){
		if($file =~ /([0-9]+)\.fa/){
			my $id = $1;
			if(exists $contigs_h{$id} ){
				# long enough (check if sam exists, if not move to no-reads directory)
				if(not -e "$id.sam"){
					system("mv $file $no_reads_contigs/");
					$no_reads++;
				}
			}
			else{
				# too short (move contig and sam to small-contigs directory)
				system("mv $file $small_contigs/");
				if(-e "$id.sam"){ system("mv $id.sam $small_contigs/"); $too_small++; }
			}
		}
	}

	# clean up remaining SAM files
	@dir_files = <*>;
	foreach my $file (@dir_files){
		if($file =~ /([0-9]+)\.sam/){
			my $id = $1;
			if(not exists $contigs_h{$id}){
				# sam file of too short contig
				system("mv $file $small_contigs/");
				$too_small++;
			}
		}
	}

	# if alignment done, clean left overs to aln dir
	if($do_aln){
	    # clean up
	    mkdir("aln") or die "could not mkdir aln, $!. Quitting...\n";
	    system("mv $contigs_name.* *.sai aln/");
	}
}

#####################################################################
sub write_summary{

	open(LOG, ">>", $logfile) or die "could not open $logfile, $!. Quitting...\n";
	my $time = localtime;
	print LOG "\n$time --- Done! (total contigs: $tot_contigs, 0-reads-aligned contigs: $no_reads, too-short contigs: $too_small)\n\n";
	print LOG "\nSummary:\n";
	print LOG "\t$tot_pe_reads read-pairs\n";
	print LOG "\t$tot_pe_aln aligned\n";
	print LOG "\t", $tot_pe_aln - $flt_on_mult, " permissively aligned ($flt_on_mult filtered on mult-map, $noX0 kept although missing X0)\n";
	print LOG "\n";
	close LOG;
	
	print "\n$time --- Done! see $outdir/$logfile for details\n\n";
}

#####################################################################
sub close_files_write_stats{
	open(STATS, ">", "stats.txt") or die "could not open stats.txt, $!. Quitting...\n";
	print STATS "#contig\t#-of-reads-uniq-aligned\n";
	foreach my $c (sort keys %files){
		close $files{$c};
		print STATS "$c\t$n_reads_c{$c}\n";
	}
	close STATS;
}

#####################################################################
sub filter_into_sam_files{
	my ($sam) = @_;
	open(SAM,"<", $sam) or die "could not open $sam, $!. Quitting...\n";
	print "Processing reads into separate (contig) SAM files...\n\n";
	my $rewind = 0;
	# read SAM header lines
	while(my $l = <SAM>){
		$rewind = length $l;
		last if($l !~ /^@/);
	}
	seek(SAM, -1*$rewind , 1);

	# process alignments
	while(my $l1 = <SAM>){
		$tot_pe_reads++;
		print "\t$tot_pe_reads read-pairs processed...\n" if($tot_pe_reads % 1_000_000 == 0);
		
		my $l2 = <SAM>;
		my @a1 = split(/\t/, $l1);
		my @a2 = split(/\t/, $l2);
		my ($cname1, $cname2) = ($a1[2], $a2[2]);
		my ($bit1,$bit2) = ($a1[1], $a2[1]);
		
		# TODO changed the follwoing part so that more reads would be included
		#      only reads that explicitely aligned to multiple places are excluded
		#      otherwise the reads will be included
		#      this should include reads that "dangle" off the edges of contigs, and read where the mate aligned to another contig
		#      also - very important - changed so that if a reads goes to a contig - it's pair also goes, 
		#      regardless if it aligned to that contig!

		#next if( ($bit1 & 4) and ($bit2 & 4) ); # neither aligned
		next if( ($cname1 eq "*") and ($cname2 eq "*") ); # TODO comment this (and uncomment previous) if does not work well

		$tot_pe_aln++;
		$l1 =~ /X0:i:([0-9]+)/; # get X0 field
		my $X0 = $1;
		
		# determine whether to write read-pair
		my $w = 1;
		if(defined $X0){
			if($X0 == 1){
				$w = 1;
			}
			else{
				$w = 0;
				$flt_on_mult++;
			}
		}
		else{
			# write reads
			$w = 1;
			$noX0++;
		}

		if($w){
			my ($F1, $F2);
			if($cname1 ne "*"){
				$F1 = get_file_handle($cname1);
				print $F1 $l1;
			}
			if($cname2 ne "*"){
				$F2 = get_file_handle($cname2);
				print $F2 $l2;
			}

			# if they did not go to the same place, send the pair along..
			if($cname1 ne $cname2){
					print $F1 $l2 if($cname1 ne "*");
					print $F2 $l1 if($cname2 ne "*");
			}
		}
		
# TODO uncomment this if does not work well
#		if(defined $X0){
#			if($X0 == 1){
#				# write read1
#				my $F1 = get_file_handle($cname1);
#				print $F1 $l1;
#				
#				# write read2
#				if(not $bit2 & 4){
#					my $F2 = get_file_handle($cname2);
#					print $F2 $l2;
#				} 
#			}
#			else{
#				$flt_on_mult++;
#			}
#		}
#		else{	
#			$noX0++;
#			#die "\nexpecting X0 field in first alignment of pair..\n$l1\n$l2\n";
#		}
	}
	close SAM;
	print "\nDone.\n";
}

#####################################################################
sub get_file_handle{
	my ($name) = @_;
	my $F;
	if(exists $files{$name}){
        $F = $files{$name};
		$n_reads_c{$name}++;
    }
	else{
		my $fname = $name;
		if($name =~ /NODE_([0-9]+)_/){
			# Velvet contig name
			$fname = "$1.sam";
		}
		elsif($name =~ /Contig_([0-9]+)/){
			# SOAP contig name
			$fname = "$1.sam";
		}
		else{
			# Euler contig name
			$fname = $name.".sam";
		}
		open($F, ">", $fname) or die "could not open $fname, $!. Quitting...\n";
        	$files{$name} = $F;
		$n_reads_c{$name} = 1;
	}
	return $F;
}

#####################################################################
sub align_all_reads_to_all_contigs{
	my ($reads1, $reads2, $contigs) = @_;
	print "\nGenerating paired-end alignments...\n\n";
	
	# align with BWA
	my ($sai1, $sai2, $sam) = ("reads1_aln.sai", "reads2_aln.sai", "reads_aln.sam");
	my $cmd;
	$cmd = "$bwa index -a is $contigs 2>/dev/null";
	system($cmd);
	
	$cmd = $bwa_aln1_cmd . " > $sai1";
	system($cmd);
	$cmd = $bwa_aln2_cmd . " > $sai2";
	system($cmd);
	
	my $bwa_sampe_cmd = "$bwa sampe $contigs $sai1 $sai2 $reads1 $reads2 > $sam";
	system($bwa_sampe_cmd);
	
	print "\nDone!\n\n";
	return $sam;
}

#####################################################################
sub sep_contigs{
        print "Extracting contigs >= $min_len bp...";
        my ($f) = @_;
        open(CONTS, "<", $f) or die "Could not open $f, $!. Quitting...\n";
        while(my $l = <CONTS>){
                if($l =~ /^>/){
						$tot_contigs++;

						# get header
                        my $header = $l;
                        chomp $header;
                        $header =~ s/^>//;
						
						# read contig sequence
						my $hold_sep = $/;
						$/ = ">";
						my $seq = <CONTS>; # read sequence						
						chop $seq;
						my $cont_len = length($seq);
						$/ = $hold_sep;
						seek(CONTS, -1, 1);

						# get contig ID
						my $cont_id;
						my ($soap, $velv) = (0, 0);
			            if($header =~ /NODE_/){
							# Velvet
							my @line = split(/_/, $header);
							$cont_id = $line[1];
							$velv = 1;
						}
						elsif($header =~ /Contig_/){
							# SOAP
							my @line = split(/_/, $header);
							$cont_id = $line[1];
							$soap = 1;
						}
						else{
							# Euler
							my @line = split(/\s+/, $header);
							$cont_id = $line[0];
						}
						
                        if($cont_len >= $min_len){
							# long enough (do not put in short-contig directory)
							$contigs_h{$cont_id} = $header;
						}
                        my $file_name = $cont_id.".fa";
                        open(FA, ">", $file_name) or die "Could not open $file_name, $!. Quitting...\n";
#						if($velv){
#							# velvet
#	                        $l =~ s/_/ /g; # changes fasta header
#        	                $l =~ s/NODE //; # changes fasta header
#						}
#						elsif($soap){
#							$l =~ s/_/ /g; # changes fasta header
#							$l =~ s/Contig //; # changes fasta header
#						}
                        print FA $l, $seq;
                        close FA;
                }
        }
        close CONTS;
        print " done.\n\n";
}

