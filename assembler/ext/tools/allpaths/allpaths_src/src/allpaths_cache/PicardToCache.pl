#!/usr/bin/perl -w
###########################################################################
#       This software and its documentation are copyright (2009) by the   #
#   Broad Institute/Massachusetts Institute of Technology.  All rights    #
#   are reserved.  This software is supplied without any warranty or      #
#   guaranteed support whatsoever. Neither the Broad Institute nor MIT    #
#   can be responsible for its use, misuse, or functionality.             #
###########################################################################
#
# PicardToCache.pl 
#
# Take Picard lanes and add them to the cache. 
# Assumes libraries have been added.  If not, aborts.
#
#
#
# This script is a wrapper for other scripts.  
#
# 2010-07   Filipe Ribeiro     ribeiro@broadinstitute.org
#



use strict;
use FindBin;

# ---- Local libraries, from elsewhere in the BroadCRD repository.
use lib "$FindBin::Bin/../";
use lib "$FindBin::Bin/";

use PerlRunTime;    # run_or_die
use PerlUtils;      # ISO_date, abort
use ArachneArgs;    # Command-line ARG=VALUE parser

use NaifDB;         # database handler
use AllPathsCache;  # definitions of database fields


sub Tag { return ISO_date() . " (PTC): "; }



# ---- CONTROL BEGINS HERE
#      Parse command-line options of the form KEY=value.
#      This function comes from the ArachneArgs.pm module.
my %args = getCommandArguments
    (PICARD_TOOLS_DIR     => { value => "",
                               help  => "The Picard tools directory." }, 
     PICARD_DIR           => { value => "",
                               help  => "The Picard sequence data directory." }, 
     CACHE_DIR            => { value => "",
                               help  => "Optional directory for temporary files."},
     TMP_DIR              => { value => "/broad/shptmp/$ENV{USER}",
                               help  => "A temporary directory." },
     GROUPS               => { value => undef,
                               help  => "The groups to import, e.g. \"{613F0AAXX.{2,3},61NCCAAXX.{7,8}}\"." },
     INCLUDE_NON_PF_READS => { value => "True",
                               help  => "Whether or not to include non-PF reads." },
     SAVE_INTERMEDIATES   => { value => "False",
                               help  => "Whether or not to keep intermediate files." },
     OVERWRITE            => { value => "False",
                               help  => "Whether or not to overwrite cache entries." },
     DB_ONLY              => { value => "False", 
                               help  => "update the databases only without reconverting BAMs." },
     HOSTS                => { value => "", 
                               help  => "e.g. \"2,4.remote1,2.remote2\"." },
     DRY_RUN              => { value => "False",
                               help  => "If 'True' only say what would be done (NOT WORKING YET!)." },
     VERBOSE              => { value => "True" });


print_action("Running on '$ENV{PWD}'.\n\n");

   

# ---- Make sure PICARD_DIR is defined 

#my $PICARD_DIR = get_directory_from_arg_or_env(\%args, "PICARD_DIR");
my $PICARD_DIR = "/seq/picard";



# ---- Make sure PICARD_TOOLS_DIR is defined and the needed Picard modules are present.

my $PICARD_TOOLS_DIR = get_directory_from_arg_or_env(\%args, "PICARD_TOOLS_DIR");

my $n_errors = 0;
foreach my $module ("SamToFastq.jar", "RevertSam.jar") 
{
    if (! -e "$PICARD_TOOLS_DIR/$module") 
    {
	print_error("Can't find Picard module '$module' in '$PICARD_TOOLS_DIR'.");
	$n_errors++;
    }
}
abort("Found $n_errors errors.")
    if ($n_errors);





# ---- Make sure CACHE_DIR is defined, exists, and is a full path.

my $CACHE_DIR = get_directory_from_arg_or_env(\%args, "CACHE_DIR");
directory_is_absolute_or_die($CACHE_DIR);

mkdir("$CACHE_DIR/tmp") unless (-d "$CACHE_DIR/tmp");




# ---- Build in_groups_db from the specified Picard pipeline groups.
#      Also save in_groups_db to a .csv file. It's the input to CacheGroups.pl

print_action(Tag() . "Finding groups in '$PICARD_DIR'");

my ($in_groups_fn, $revert_sam_arg) = build_groups_db($args{GROUPS});






# ---- Call CacheGroups.pl to convert BAMs to fastb/qualb 

my $module = "CacheGroups.pl";
print_action(Tag() . "Running $module to convert BAM files to AllPaths cache:");
print "\n";

my $fw_args = join " ", (map "$_=$args{$_}", ("TMP_DIR", 
                                              "INCLUDE_NON_PF_READS",
                                              "SAVE_INTERMEDIATES",
                                              "OVERWRITE",
                                              "DB_ONLY",
                                              "HOSTS", 
                                              "DRY_RUN", 
                                              "VERBOSE"));

run_or_die("$FindBin::Bin/$module" .
           " ACTION=Add" .
           " PICARD_TOOLS_DIR=$PICARD_TOOLS_DIR" .
           " CACHE_DIR=$CACHE_DIR" .
           " IN_GROUPS_CSV=$in_groups_fn" .
           " REVERT_SAM='$revert_sam_arg'" .
           " $fw_args");
print "\n";











sub build_groups_db
{
    my @group_names = array_from_ref_or_value($_[0]);


    # ---- Open the cache library database from a .csv file.

    my $libs_db_fn = "$CACHE_DIR/libraries.csv";
    print_action(Tag() . "Loading cached libraries '$libs_db_fn'.");
    my $libs_db = new_naif_db()->from_csv_or_continue($libs_db_fn, $libs_db_field_validator);


    # ---- Cycle through all input groups.

    my $in_groups_db = new_naif_db($in_groups_db_field_validator);
    $n_errors = 0;
    foreach my $group_name (@group_names) 
    {
        my ($flowcell, $lane, $group_lib); 
        if ($group_name =~ /^(.+)\.(\d+)$/) {
            ($flowcell, $lane, $group_lib) = ($1, $2, "*");
        }
        elsif ($group_name =~ /^(.+)\.(\d+)\.(.+)$/) {
            ($flowcell, $lane, $group_lib) = ($1, $2, $3);
        }
        else
        {
            abort(Tag() . "Group '$group_name' is not of form '<flowcell>.<lane>[.<library>]'");
        }
         
        print_action(Tag() . "Group '$group_name':");
        # /seq/picard/<flowcell>/<date>/<lane>/<library>/<flowcell>.<lane>.*.bam
        my @bams = (glob("$PICARD_DIR/$flowcell/*/$lane/$group_lib/$flowcell.$lane.unmapped.bam"),
                    glob("$PICARD_DIR/$flowcell/*/$lane/$group_lib/$flowcell.$lane.aligned.duplicates_marked.bam"));

        if (@bams) 
        {
            # ---- Separate BAMs by libraries.

            my %libs_to_bams = ();
            foreach my $bam (@bams) 
            {
                my ($lib) = $bam =~ /$lane\/([^\/]+)\/$flowcell/;
                $libs_to_bams{$lib} = [] if (!exists $libs_to_bams{$lib});
                push @{$libs_to_bams{$lib}}, $bam;
            }
            my @libs = keys %libs_to_bams;
            my $n_libs = scalar @libs;

            foreach my $lib (@libs) 
            {
                # ---- Check if lib is defined in the cache.
                if ($libs_db->sub_db_select_count("library_name", $lib)) 
                {
                    print_subaction("Library '$lib' is defined in the cache library database.");
                    # ---- Choose which BAM to import.

                    my $file_name = choose_bam($libs_to_bams{$lib});
                    if ($file_name) 
                    {
                        # ---- some lanes and flowcells have several libraries in them
                        my $group = ($n_libs == 1) ? $group_name : "$group_name.$lib";

                        my $group_rec = { group_name   => $group,
                                          library_name => $lib,
                                          file_name    => $file_name };
                        $in_groups_db->record_add($group_rec);
                    }
                    else {
                        print_error("Can't choose bam files in library '$lib':");
                        map { print "     $_\n" } @{$libs_to_bams{$lib}};
                        $n_errors++;
                    }
                }
                else 
                {
                    print_error("Library '$lib' is not defined in the cache library database.");
                }
            }
        }
        else {
            print_error("Didn't find any file of type '$PICARD_DIR/$flowcell/*/$lane/*/$flowcell.$lane.{unmapped,aligned.duplicates_marked}.bam'");
            $n_errors++;
        }
    }
    abort(Tag() . "Found $n_errors errors.")
        if ($n_errors);


    abort(Tag() . "Didn't find any groups.") 
        unless ($in_groups_db->size() > 0);


    my $date_str = ISO_date();
    $date_str =~ s/[-\:\s]/_/g;
    my $in_groups_fn = sprintf "$CACHE_DIR/tmp/$date_str.rnd%04d.in_groups.csv", rand(10000);
    
    $in_groups_db->to_csv($in_groups_fn);


    # ---- Figure out on which files it is necessary to remove alignment information.
    
    my @revert_sam = map { (/aligned/) ? "True" : "False"; } @{$in_groups_db->{file_name}};

    my $revert_sam_arg = globable_str_from_array_or_value(@revert_sam);





    return ($in_groups_fn, $revert_sam_arg);
}












sub choose_bam
{
    my ($bams) = @_;

    my @unmapped_bams = grep /unmapped/, @$bams;
    my @aligned_bams = grep /aligned/, @$bams;

    my $file_name = "";
    if (@unmapped_bams >= 1) {
        if (@unmapped_bams > 1) {
            print_warning("Found more than one unmapped bam files. Choosing latest.");
            map { print "     $_\n"; } @unmapped_bams;
        }

        $file_name = latest_bam(@unmapped_bams);
    }
    elsif (@aligned_bams >= 1) {
        if (@aligned_bams > 1) {
            print_warning("Found more than one aligned bam files. Choosing latest.");
            map { print "     $_\n"; } @aligned_bams;
        }

        $file_name = latest_bam(@aligned_bams);
    }

    # ---- Checked if 'finished.txt' is present
    my $finished = $file_name;
    $finished =~ s/\/[^\/]+$/\/finished.txt/;
    if (!-e $finished) {
        print_error("Couldn't find '$finished'! Skipping the associated BAM file.");
        return "";
    }

    return $file_name;
}







sub latest_bam
{
    my @bams = @_;

    my @dates = map { (/picard\/[^\/]+\/([^\/]+)\//)[0] } @bams;

    my $i_last = (sort { $dates[$a] cmp $dates[$b] } (0..@bams-1))[-1];

    return $bams[$i_last];
}





sub read_params_txt
{
    my ($fn) = @_;

    open F, "<$fn" or abort("Can't open '$fn'.");
    my %p;
    while (my $l = <F>) {
        if ($l =~ /^([^=]+)\=(.+)$/) {
            my ($key, $val) = ($1, $2);
            $val =~ s/,//g;
            $p{$key} = $val;
        }
    }
    close F;
    return \%p;
}

