#!/usr/bin/perl -w
###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2011) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################
#
# ConvertToFastbQualb.pl
#
# Takes a bam, fastq, or fasta file and converts it to a fastb and qualb file
# (<OUT_HEAD>.{fastb,qualb}). 
#
# 2010-05   Josh Burton    
# 2010-06   Filipe Ribeiro    <crdhelp@broadinstitute.org>
#

use strict;
use FindBin;

# Local libraries, from elsewhere in the BroadCRD repository.
use lib "$FindBin::Bin/../";
use lib "$FindBin::Bin/";

use PerlRunTime;    # run_or_die
use PerlUtils;      # ISO_date, abort
use ArachneArgs;    # Command-line ARG=VALUE parser




sub Tag { return ISO_date() . " (CTFQ): "; }



# CONTROL BEGINS HERE
# Parse command-line options of the form KEY=value.
# This function comes from the ArachneArgs.pm module.

my %args = getCommandArguments
    (PICARD_TOOLS_DIR      => { value => undef,
                                help  => "The Picard tools directory." }, 
     TMP_DIR               => { value => "",
                                help  => "A temporary directory (should be large for large data sets)." },
     DATA_FILE             => { value => undef,
                                help  => "The data file(s) to convert ('*.{bam,fast{a,b,q}}')." },
     OUT_HEAD              => { value => undef,
                                help  => "The prefix of the output .fastb and .qualb files." },
     PAIRED                => { value => undef,
                                help  => "Whether files are paired or not. If paired, interleave files when merging." },
     INCLUDE_NON_PF_READS  => { value => "True",
                                help  => "Whether or not to include non-PF reads." },
     PHRED_64              => { value => "False",
                                help  => "True: fastq quals are encoded with PHRED 64. False: PHRED 33." },
     FORCE_PHRED           => { value => "False",
                                help  => "True: accepts specified PHRED encoding. False: Aborts is incorrect PHRED detected." },
     OVERWRITE             => { value => "False",
                                help  => "Whether or not to overwrite cache entries." },
     SAVE_COMPRESSED_FASTQ => { value => "False",
                                help  => "Whether or not to save compressed copies of fastq files." },
     SAVE_INTERMEDIATES    => { value => "False",
                                help  => "Whether or not to keep intermediate files." },
     REVERT_SAM            => { value => "False",
                                help  => "Whether to revert alignment information in SAM files." },
     REVERSE_READS         => { value => "False",
                                help  => "Whether to reverse reads or not." },
     TRIM_START            => { value => 0,
                                help  => "Trim all bases in a read before this one." }, 
     TRIM_END              => { value => 0,
                                help  => "Trim all bases in a read after this one ('0' don't trim)." }, 
     JAVA_MEM_GB           => { value => 8,
                                help  => "Memory to reserve for java (only for SAM/BAM conversion)." }, 
     DRY_RUN               => { value => "False",
                                help  => "If 'True' only say what would be done." },
     VERBOSE               => { value => "True" });

 
set_verbose($args{VERBOSE});

my $fn_in = $args{DATA_FILE};
my $head_fastb_out = $args{OUT_HEAD};


my $picard_tool_root = "java -Xmx" . $args{JAVA_MEM_GB} . "g -jar $args{PICARD_TOOLS_DIR}";


# ---- Look for previously computed files and deal with them appropriately.
clear_old(map "$head_fastb_out.$_", (".fastb", ".qualb"));


# ---- directory to store temporary files (different from $args{TMP_DIR})
my $dir_tmp = "$args{OUT_HEAD}.tmp";
mkdir $dir_tmp unless (-d $dir_tmp);


# =================
#        BAM 
# =================

if ($fn_in =~ /\.bam$/i) 
{
    abort(Tag() . "Must specify PAIRED = 1 with BAM files.")
        if (!$args{PAIRED});

    my $fn_BAM = $fn_in;

    # ---- Remove alignment information from SAM files (if there is any)
    my $fn_BAM_reverted = "$dir_tmp/reverted.bam";
    if ($args{REVERT_SAM}) {
        sam_revert($fn_BAM, $fn_BAM_reverted);
        $fn_BAM = $fn_BAM_reverted;
    }
    my @fns_fastq_tmp = sam_to_paired_fastqs($fn_BAM, $dir_tmp);
    my @heads_tmp = fastqs_to_fastbs(\@fns_fastq_tmp, $dir_tmp);
    fastbs_merge_trim_reverse(\@heads_tmp, $head_fastb_out, $dir_tmp, $args{PAIRED}, $fn_in);
}

# ===================
#        fastq 
# ===================

elsif ($fn_in =~ /\.f(ast)?q(\.gz)?$/i) 
{
    my @fns_fastq_in = ($fn_in);

    if ($fn_in =~ /[\*\?]/) {
        @fns_fastq_in = sort glob $fn_in;
        
        abort(Tag() . "No files found with glob pattern '$fn_in'.")
            if (@fns_fastq_in == 0);
    }
    my $n_errors = 0;
    foreach my $fn (@fns_fastq_in) {
        if (!-e $fn) {
            print_error(Tag() . "Can't find '$fn'.");
            $n_errors++;
        }
    }
    exit(1) unless ($n_errors == 0);
    
    my @heads_tmp = fastqs_to_fastbs(\@fns_fastq_in, $dir_tmp);
    fastbs_merge_trim_reverse(\@heads_tmp, $head_fastb_out, $dir_tmp, $args{PAIRED}, $fn_in);
}

# ===================
#        fasta 
# ===================

elsif ($fn_in =~ /\.f(ast)?a$/i) 
{
    my @fns_fasta_in = ($fn_in);

    if ($fn_in =~ /[\*\?]/) {
        @fns_fasta_in = sort glob $fn_in;
        abort(Tag() . "No files found with glob pattern '$fn_in'.")
            if (@fns_fasta_in == 0);
    }
    my $n_errors = 0;
    foreach my $fn (@fns_fasta_in) {
        if (!-e $fn) {
            print_error(Tag() . "Can't find '$fn'.");
            $n_errors++;
        }
    }
    exit(1) unless ($n_errors == 0);
    
    my @heads_tmp = fastas_to_fastbs(\@fns_fasta_in, $dir_tmp);
    fastbs_merge_trim_reverse(\@heads_tmp, $head_fastb_out, $dir_tmp, $args{PAIRED}, $fn_in);

}

# ===================
#        fastb
# ===================

elsif ($fn_in =~ /\.fastb$/i) # ====== fastb ======
{
    my @fns_fastb_in = ($fn_in);

    if ($fn_in =~ /[\*\?]/) {
        @fns_fastb_in = sort glob $fn_in;
        print_error(Tag() . "No files found with glob pattern '$fn_in'.")
            if (@fns_fastb_in == 0);
    }
    my $n_errors = 0;
    foreach my $fn (@fns_fastb_in) {
        if (!-e $fn) {
            print_error(Tag() . "Can't find '$fn'.");
            $n_errors++;
        }
    }
    exit(1) unless ($n_errors == 0);


    my @heads_fastb_in = map /^(.+).fastb/i, @fns_fastb_in;
    fastbs_merge_trim_reverse(\@heads_fastb_in, $head_fastb_out, $dir_tmp, $args{PAIRED}, $fn_in);
}


# ==================
#        what?
# ==================

else {
    abort(Tag() . "File '$fn_in' is neither a '.bam','.fastq', '.fasta', nor 'fastb' file.");
}
                    
                        







# ---- Clean up
if (!$args{SAVE_INTERMEDIATES} && !$args{DRY_RUN}) 
{
    unlink (glob "$dir_tmp/*");
    rmdir $dir_tmp or print_error("Can't remove directory '$dir_tmp'.") if (-d $dir_tmp);
}



# ---- Signal to modules downstream that conversion was successful
system "touch", "$args{OUT_HEAD}.done" if (-e "$head_fastb_out.fastb" && !$args{DRY_RUN});













# ==== Subroutine definitions.





sub clear_old 
{
    foreach my $fn (@_) {
        if (-e $fn) {
            print_warning("Found previous '$fn'.");
            if ($args{OVERWRITE}) {
                if ($args{DRY_RUN}) {
                    print_warning("Would erase '$fn' if DRY_RUN was false.");
                }
                else {
                    print_warning("Erasing '$fn'.");
                    unlink $fn;
                }
            }
            else {
                abort("Found previously computed '$fn'.");
            }
        }
    }
}




# ---- RevertSam ----
#
# Remove alignment information from the SAM/BAM file that shouldn't
# be known by a de novo assembler.
#
# Info removed:
# -- Re-ordering the reads by their reference alignment
# -- Re-calibrating the quality scores in accordance with reference alignments
#
# The alignments themselves are kept, but won't be passed into the fastq file.
#
# input:  group.bam
# output: $dir_tmp/reverted.bam
#
sub sam_revert
{
    my ($BAM_fn, $BAM_reverted_fn, $dir_tmp) = @_;
    
    #print_action(Tag() . "Reverting SAM alignment information.");
    
    my $cmd = ("$picard_tool_root/RevertSam.jar" .
               " I=$BAM_fn" .
               " O=$BAM_reverted_fn" . 
               " SO=queryname" .
               " REMOVE_DUPLICATE_INFORMATION=true" .
               " REMOVE_ALIGNMENT_INFORMATION=true" .
               " TMP_DIR=$args{TMP_DIR}");
    
    run_cmd($cmd, $args{VERBOSE}, $args{DRY_RUN});

    abort(Tag() . "RevertSam call failed to generate '$BAM_reverted_fn'.")
        unless (-e $BAM_reverted_fn);
}




# ---- SamToFastq ----
#
# Format conversion: SAM/BAM -> FASTQ 
# Two separate files: - one for the 1st read of each read pair, and one for the 2nd
#
# input: $fn_BAM
# output: $dir_tmp/{A,B}.fastq
#
sub sam_to_paired_fastqs 
{
    my ($fn_BAM, $dir_tmp) = @_;
    #print_action(Tag() . "Converting BAM to paired fastqs.");

    my @fns_fastq_tmp = map "$dir_tmp/$_.fastq", ("A", "B");

    my $cmd = ("$picard_tool_root/SamToFastq.jar".
               " I=$fn_BAM" .
               " INCLUDE_NON_PF_READS=" . bool_str($args{INCLUDE_NON_PF_READS}) . 
               " F=$fns_fastq_tmp[0]" . 
               " F2=$fns_fastq_tmp[1]" . 
               " TMP_DIR=$args{TMP_DIR}");
    
    run_cmd($cmd, $args{VERBOSE}, $args{DRY_RUN});
    
    foreach (@fns_fastq_tmp) {
        abort(Tag() . "SamToFastq call failed to generate '$_'.")
            unless (-e $_);
    }

    if ($args{SAVE_COMPRESSED_FASTQ} && !$args{DRY_RUN}) {
        print `gzip < $fns_fastq_tmp[0] > $args{OUT_HEAD}.A.fastq.gz`;
        print `gzip < $fns_fastq_tmp[1] > $args{OUT_HEAD}.B.fastq.gz`;
    }

    return @fns_fastq_tmp;
}






# ---- FastqToFastbQualb ----
#
# Format conversion: FASTQ -> FASTB/QUALB
#
# input: $fns_fastq
# output: $heads_fastb_tmp).{fastb,qualb}
#
sub fastqs_to_fastbs
{
    my ($fns_fastq, $dir_tmp) = @_;

    #print_action(Tag() . "Converting fastqs to fastbs/qualbs.");

    my $n = @$fns_fastq;
    my @heads_tmp = map "$dir_tmp/$_", (0..$n-1);
    for (my $i = 0; $i != $n; $i++) {
        run_cmd("$FindBin::Bin/FastqToFastbQualb" . 
                " FASTQ=$fns_fastq->[$i]" . 
                " OUT_HEAD=$heads_tmp[$i]" .
                " PHRED_64=" . bool_str($args{PHRED_64}) .
                " FORCE_PHRED=" . bool_str($args{FORCE_PHRED}) .
                " NH=" . bool_str(!$args{VERBOSE}),
                $args{VERBOSE}, $args{DRY_RUN});

        abort(Tag() . "FastqToFastbQualb call failed to generate '$heads_tmp[$i].fastb'.")
            unless (-e "$heads_tmp[$i].fastb");
    }
    return @heads_tmp; 
}




# ---- Fasta2Fastb ---- Quala2Qualb ----
#
# Convert fasta/quala to fastb/qualb
#
# input: \@fns_fasta
# output: $tmp_dir/$i.{fastb,qualb}
#
sub fastas_to_fastbs
{
    my ($fns_fasta, $dir_tmp) = @_;
    
    #print_action(Tag() . "Converting fastas to fastbs.");
    
    my $n = @$fns_fasta;
    my @heads_tmp = map "$dir_tmp/$_", (0..$n-1);
    for (my $i = 0; $i != $n; $i++) {
        
        run_cmd("$FindBin::Bin/Fasta2Fastb" .
                " IN=$fns_fasta->[$i]" .
                " OUT=$heads_tmp[$i].fastb" .
                " NAMES=False" .
                " NH=" . bool_str(!$args{VERBOSE}),
                $args{VERBOSE}, $args{DRY_RUN});

        abort(Tag() . "Fasta2Fastb call failed to generate '$heads_tmp[$i].fastb'.")
            unless (-e "$heads_tmp[$i].fastb");

        my $fn_quala = $fns_fasta->[$i];
        $fn_quala =~ s/\.fasta$/.quala/;
        $fn_quala =~ s/\.fa$/.qa/;
        if (-e $fn_quala) {
            #print_action(Tag() . "Converting quala to qualb.");
            
            run_cmd("$FindBin::Bin/Quala2Qualb" .
                    " IN=$fn_quala" .
                    " OUT=$heads_tmp[$i].qualb" .
                    " NH=" . bool_str(!$args{VERBOSE}),
                    $args{VERBOSE}, $args{DRY_RUN});

            abort(Tag() . "Quala2Qualb call failed to generate '$heads_tmp[$i].qualb'.")
                unless (-e "$heads_tmp[$i].qualb");
        }
        else {
            print_warning(Tag() . "Can't find quality score file '$fn_quala'.");
        }
    }
    return @heads_tmp;
}








# ---- MergePairedFastbs ---- if only 2 files
# ---- FastbMerge ---- otherwise
# ---- FastbQualbTrimReverse ----
#
# Merge separate fastbs into a single fastb.  Same for qualbs.
#
# input: $heads_fastb_in.{fastb,qualb}
# output: $head_fastb_out.{fastb,qualb}
#
sub fastbs_merge_trim_reverse
{
    my ($heads_in, $head_out, $dir_tmp, $paired, $fn_in) = @_;

    my $head_tmp = "$dir_tmp/all";

    # ---- Check pairing and number of files consistency
    #      Note: paired = 1 and @$heads_in = 1 is NOT a problem,
    #            hence the '> 2' comparison below. 

    abort(Tag() . "PAIRED = 1 but there are " . (scalar @$heads_in) . " files to merge (!= 2).")
        if ($paired && @$heads_in > 2);

    # ---- Merge two files in an interleaved way

    if ($paired && @$heads_in == 2) {
        
        run_cmd("$FindBin::Bin/MergePairedFastbs" . 
                " HEAD1_IN=$heads_in->[0]" . 
                " HEAD2_IN=$heads_in->[1]" . 
                " HEAD_OUT=$head_tmp" . 
                " NH=" . bool_str(!$args{VERBOSE}),
                $args{VERBOSE}, $args{DRY_RUN});

        abort(Tag() . "MergePairedFastbs call failed to generate '$head_tmp.fastb'.")
            unless (-e "$head_tmp.fastb");
        abort(Tag() . "MergePairedFastbs call failed to generate '$head_tmp.qualb'.")
            unless (-e "$head_tmp.qualb" || !-e "$heads_in->[0].qualb");

        #run_cmd("$FindBin::Bin/MergePairedFastbs" . 
        #        " IN1=$heads_in->[0]" . 
        #        " IN2=$heads_in->[1]" . 
        #        " OUT=$head_out" . 
        #        " REVERSE_READS=" . bool_str($args{REVERSE_READS}) . 
        #        " TRIM_START=$args{TRIM_START}" . 
        #        " TRIM_END=$args{TRIM_END}" .
        #        " NH=" . bool_str(!$args{VERBOSE}),
        #        $args{VERBOSE}, $args{DRY_RUN});
    }

    # ---- Merge files (if more than one) by concatenating them together

    else {
        if ($paired && @$heads_in == 1) {
            print_warning(Tag() . "'$fn_in': PAIRED GROUP WITH ONLY ONE DATA FILE! This is OK if and only if reads are INTERLEAVED!");
        }
        if (@$heads_in == 1) {
            $head_tmp = $heads_in->[0];
        } 
        else {
            my $heads_in_arg = globable_str_from_array_or_value(@$heads_in);    
            
            run_cmd("$FindBin::Bin/FastbMerge" . 
                    " HEADS_IN='$heads_in_arg'" . 
                    " HEAD_OUT=$head_tmp" . 
                    " NH=" . bool_str(!$args{VERBOSE}),
                    $args{VERBOSE}, $args{DRY_RUN});

            abort(Tag() . "FastbMerge call failed to generate '$head_tmp.fastb'.")
                unless (-e "$head_tmp.fastb");
            abort(Tag() . "FastbMerge call failed to generate '$head_tmp.qualb'.")
                unless (-e "$head_tmp.qualb" || !-e "$heads_in->[0].qualb");
        }
    }
    

    # ---- Apply the trimming and reversing

    run_cmd("$FindBin::Bin/FastbQualbTrimReverse" .
            " IN_HEAD=$head_tmp" .
            " OUT_HEAD=$head_out" .
            " TRIM_START=$args{TRIM_START}" . 
            " TRIM_END=$args{TRIM_END}" . 
            " REVERSE=" . bool_str($args{REVERSE_READS}) . 
            " NH=" . bool_str(!$args{VERBOSE}), 
            $args{VERBOSE}, $args{DRY_RUN});

    abort(Tag() . "FastbQualbTrimReverse call failed to generate '$head_out.fastb'.")
        unless (-e "$head_out.fastb");
    abort(Tag() . "FastbQualbTrimReverse call failed to generate '$head_out.qualb'.")
        unless (-e "$head_out.qualb" || !-e "$head_tmp.qualb");


}










sub run_cmd 
{
    my ($cmd, $verbose, $dry_run) = @_;

    abort(Tag() . "Undefined command to run.")
        unless (defined $cmd);

    my $exe = (split " ", $cmd)[0];

    abort(Tag() . "Can't find executable '$exe'.")
        unless (`which $exe`);

    if ($dry_run) {
        print "DRY_RUN: $cmd\n";
        return;
    }
  
    if ($verbose) {
        run_or_die($cmd); # stdout and stderr are kept
    }
    else {
        my $output = run_or_die($cmd); # stdout and stderr go into $output
    } 
}




