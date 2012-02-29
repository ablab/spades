#!/usr/bin/perl -w
###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2011) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################
#
# PrepareAllPathsInputs.pl 
#
# Take sequence data, described in 'in_groups.csv', belonging to libraries, described in
# 'in_libs.csv', and convert them to input fastb/qualb files for the ALLPATHS pipeline.
#
# This script is a wrapper for other scripts.  
#
# 2010-05   Filipe Ribeiro
#



use strict;
use FindBin;

# ---- Local libraries, from elsewhere in the BroadCRD repository.
use lib "$FindBin::Bin/../";
use lib "$FindBin::Bin/";

use PerlRunTime;    # run_or_die
use PerlUtils;      # ISO_date, abort
use ArachneArgs;    # Command-line ARG=VALUE parser


sub Tag { return ISO_date() . " (PAPI): "; }



# ---- CONTROL BEGINS HERE
#      Parse command-line options of the form KEY=value.
#      This function comes from the ArachneArgs.pm module.
my %args = getCommandArguments
    (PICARD_TOOLS_DIR      => { value => "",
                                help  => "The Picard tools directory." }, 
     DATA_DIR              => { value => undef,
                                help  => "Directory in the ALLPATHS pipeline for input reads." },
     TMP_DIR               => { value => "",
                                help  => "A temporary directory." },
     PLOIDY                => { value => "",
                                help  => "The ploidy." },
     IN_GROUPS_CSV         => { value => "in_groups.csv",
                                help  => "Comma-separated-value info on read groups to import." }, 
     IN_LIBS_CSV           => { value => "in_libs.csv",
                                help  => "Comma-separated-value info on libraries to import." }, 
     GENOME_SIZE           => { value => "",
                                help  => "OPTIONAL estimated total genome size (downsampling purposes only)." },
     FRAG_COVERAGE         => { value => "",
                                help  => "OPTIONAL fragment library desired coverage (downsampling purposes only; must specify GENOME_SIZE)." },
     JUMP_COVERAGE         => { value => "",
                                help  => "OPTIONAL jumping library desired coverage (downsampling purposes only; must specify GENOME_SIZE)." },
     LONG_JUMP_COVERAGE    => { value => "",
                                help  => "OPTIONAL long jumps library desired coverage (downsampling purposes only; must specify GENOME_SIZE)." },
     FRAG_FRAC             => { value => "",
                                help  => "OPTIONAL fragment library desired read fraction (downsampling purposes only)." },
     JUMP_FRAC             => { value => "",
                                help  => "OPTIONAL jumping library desired read fraction (downsampling purposes only)." },
     LONG_JUMP_FRAC        => { value => "",
                                help  => "OPTIONAL long jumps library desired read fraction (downsampling purposes only)." },
     LONG_JUMP_MIN_SIZE    => { value => 20000,
                                help  => "Threshold for long jumps." },
     LONG_READ_MIN_LEN     => { value => 500,
                                help  => "Threshold for long reads." },
     INCLUDE_NON_PF_READS  => { value => "True",
                                help  => "Whether or not to include non-PF reads." },
     PHRED_64              => { value => "False",
                                help  => "True: fastq quals are encoded with PHRED 64. False: PHRED 33." },
     FORCE_PHRED           => { value => "False",
                                help  => "True: accepts specified PHRED encoding. False: Aborts is incorrect PHRED detected." },
     OVERWRITE             => { value => "False",
                                help  => "Whether or not to overwrite cache entries." },
     SAVE_INTERMEDIATES    => { value => "False",
                                help  => "Whether or not to keep intermediate files." },
     HOSTS                 => { value => "", 
                                help  => "e.g. \"2,4.remote1,2.remote2\"." },
     JAVA_MEM_GB           => { value => 8,
                                help  => "Memory to reserve for java (only for SAM/BAM conversion)." }, 
     DRY_RUN               => { value => "False",
                                help  => "If 'True' only say what would be done." },
     VERBOSE               => { value => "True" });


set_verbose($args{VERBOSE});

my $DRY_RUN = $args{DRY_RUN};



# ---- Some argument error checking 
abort(Tag() . "Specified DATA_DIR '$args{DATA_DIR}' is not an absolute path.") 
    unless ($args{DATA_DIR} =~ /^\//);

abort(Tag() . "Specified DATA_DIR '$args{DATA_DIR}' does not exist.")
    unless (-d $args{DATA_DIR});

abort(Tag() . "IN_GROUPS_CSV database file '$args{IN_GROUPS_CSV}' does not exist.") 
    unless (-e $args{IN_GROUPS_CSV});

abort(Tag() . "IN_LIBS_CSV database file '$args{IN_LIBS_CSV}' does not exist.") 
    unless (-e $args{IN_LIBS_CSV});





# ---- Make sure *_COVERAGE, *_FRAC, GENOME_SIZE options make sense
{
    my $errors = 0;
    if ($args{FRAG_COVERAGE} ne "" && $args{FRAG_FRAC} ne "") { 
        print_error(Tag() .  "Can't specify both FRAG_COVERAGE and FRAG_FRAC.");
	$errors++;
    }
    if ($args{JUMP_COVERAGE} ne "" && $args{JUMP_FRAC} ne "") { 
        print_error(Tag() .  "Can't specify both JUMP_COVERAGE and JUMP_FRAC.");
	$errors++;
    }
    if ($args{LONG_JUMP_COVERAGE} ne "" && $args{LONG_JUMP_FRAC} ne "") { 
        print_error(Tag() .  "Can't specify both LONG_JUMP_COVERAGE and LONG_JUMP_FRAC.");
	$errors++;
    }
    if ($args{FRAG_COVERAGE} ne "" && !$args{GENOME_SIZE}) {
        print_error(Tag() . "Must specify GENOME_SIZE if specifying FRAG_COVERAGE.");
	$errors++;
    }
    if ($args{JUMP_COVERAGE} ne "" && !$args{GENOME_SIZE}) {
        print_error(Tag() . "Must specify GENOME_SIZE if specifying JUMP_COVERAGE.");
	$errors++;
    }
    if ($args{LONG_JUMP_COVERAGE} ne "" && !$args{GENOME_SIZE}) {
        print_error(Tag() . "Must specify GENOME_SIZE if specifying LONG_JUMP_COVERAGE.");
	$errors++;
    }
    abort(Tag() . "Found $errors errors. See above.") 
        if ($errors);
}



# ---- Check if cache directory exists and create it if it doesn't

my $CACHE_DIR = "$args{DATA_DIR}/read_cache";
if (!-d $CACHE_DIR) {
    print_action(Tag() . "Creating cache directory '$CACHE_DIR'");
    mkdir $CACHE_DIR unless ($DRY_RUN);
}
else {
    print_action(Tag() . "Using cache directory '$CACHE_DIR'");
}
abort(Tag() . "Can't find cache directory '$CACHE_DIR'.") 
    unless (-d $CACHE_DIR);



# ---- Tee the stdout and stderr
{
    my $dir_log = "$CACHE_DIR/log";
    mkdir $dir_log unless (-e $dir_log);
    my $prefix = `date +"\%Y-\%m-\%d.\%H:\%M:\%S"`;
    chomp $prefix;
    my $fn_out = "$dir_log/$prefix.PAPI.out";
    my $fn_err = "$dir_log/$prefix.PAPI.err";
    open TEMP, ">$fn_out"; PrintCommandHeader(*TEMP); close TEMP;
    open TEMP, ">$fn_err"; PrintCommandHeader(*TEMP); close TEMP;
    open STDOUT, "| tee -a $fn_out" or die "Can't tee stdout to '$fn_out': $!\n";
    open STDERR, "| tee -a $fn_err" or die "Can't tee stderr to '$fn_err': $!\n";
}







print_action(Tag() . "Running on '$ENV{PWD}'.");




# ---- First merge the libraries 
if (1) {
    my $module = "CacheLibs.pl";
    my $fw_args = join " ", (map "$_=$args{$_}", ("OVERWRITE",
                                                  "DRY_RUN", 
                                                  "VERBOSE"));
    
    #print_action(Tag() . "Running '$module' to merge input libraries with cached libraries.");
    run_or_die("$FindBin::Bin/$module" .
               " CACHE_DIR=$CACHE_DIR" .
               " ACTION=Add" .
               " IN_LIBS_CSV=$args{IN_LIBS_CSV}" .
               " $fw_args");
    #print "\n";
}





# ---- First convert data files to fastb/qualb files in the cache
if (1) {
    my $module = "CacheGroups.pl";
    my $fw_args = join " ", (map "$_=$args{$_}", ("IN_GROUPS_CSV", 
                                                  "INCLUDE_NON_PF_READS",
                                                  "PHRED_64",
                                                  "FORCE_PHRED",
                                                  "OVERWRITE",
                                                  "SAVE_INTERMEDIATES",
                                                  "TMP_DIR",
                                                  "HOSTS", 
                                                  "JAVA_MEM_GB",
                                                  "DRY_RUN", 
                                                  "VERBOSE"));


    #print_action(Tag() . "Running $module to import data files to a temporary cache.");
    run_or_die("$FindBin::Bin/$module" .
               " ACTION=Add" .
               " PICARD_TOOLS_DIR=$args{PICARD_TOOLS_DIR}" .
               " CACHE_DIR=$CACHE_DIR" .
               " $fw_args");
    #print "\n";
}





# ---- Now merge the cached reads
if (1) {
    my $module = "CacheToAllPathsInputs.pl";
    my $fw_args = join " ", (map "$_=$args{$_}", ("PLOIDY",
                                                  "IN_GROUPS_CSV", 
                                                  "DATA_DIR",
                                                  "GENOME_SIZE", 
                                                  "FRAG_COVERAGE",
                                                  "FRAG_FRAC",
                                                  "JUMP_COVERAGE",
                                                  "JUMP_FRAC",
                                                  "LONG_JUMP_COVERAGE",
                                                  "LONG_JUMP_FRAC",
						  "LONG_JUMP_MIN_SIZE",
						  "LONG_READ_MIN_LEN",
                                                  "DRY_RUN", 
                                                  "VERBOSE"));
    
    #print_action(Tag() . "Running '$module' to merge cached reads into ALLPATHS input files.");
    run_or_die("$FindBin::Bin/$module" .
               " CACHE_DIR=$CACHE_DIR" .
               " $fw_args");
    #print "\n";
}




# ---- Erase the cache

#print_action(Tag() . "Erasing the cache in '$CACHE_DIR'.");
#system("rm -fr $CACHE_DIR") 
#    unless ($DRY_RUN);




print_action(Tag() . "Done.\n");


