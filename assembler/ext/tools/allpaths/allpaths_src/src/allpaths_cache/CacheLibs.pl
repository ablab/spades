#!/usr/bin/perl -w
###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2011) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################
#
# CacheLibs.pl
#
# Manipulate the cache libraries.
#
# 2010-06   Filipe Ribeiro
#


use strict;
use FindBin;

# Local libraries, from elsewhere in the BroadCRD repository.

use lib "$FindBin::Bin/../";
use lib "$FindBin::Bin/";
use PerlRunTime;    # run_or_die
use PerlUtils;      # ISO_date, abort
use ArachneArgs;    # Command-line ARG=VALUE parser

use NaifDB;         # Filipe's database handler
use AllPathsCache;  # definitions of database fields



sub Tag { return ISO_date() . " (CL): "; }




# ---- CONTROL BEGINS HERE
#      Parse command-line options of the form KEY=value.
#      This function comes from the ArachneArgs.pm module.

my %args = getCommandArguments
    (CACHE_DIR            => { value => "",
                               help  => "Directory where reads in AllPaths format are cached."},
     ACTION               => { value => "List",
                               help  => "Action to perform: List, Add, Remove." },
     IN_LIBS_CSV          => { value => "",
                               help  => "Comma-separated-value  info on libraries to add." },
     OVERWRITE            => { value => "False",
                               help  => "Whether or not to overwrite cache entries." },
     PROJECT              => { value => "",
                               help  => "Optional project name of reads." },
     ORGANISM             => { value => "",
                               help  => "Optional organism name of reads." },
     LIBRARIES            => { value => "",
                               help  => "Optional libraries names of reads." },
     TYPE                 => { value => "",
                               help  => "Optional type ('frag' or 'jump') of reads to import." },
     DRY_RUN              => { value => "False",
                               help  => "If 'True' only say what would be done." },
     VERBOSE              => { value => "True" });
                               

set_verbose($args{VERBOSE});



# ---- Make sure cache dir is defined. 

my $CACHE_DIR = get_directory_from_arg_or_env(\%args, "CACHE_DIR");
directory_is_absolute_or_die($CACHE_DIR);





# ---- Do stuff

if ($args{ACTION} eq "List") 
{
    libraries_list(); 
}
elsif ($args{ACTION} eq "Add") 
{
    libraries_add();
}
elsif ($args{ACTION} eq "Remove") 
{
    libraries_remove();
}
else 
{
    abort(Tag() . "ACTION '$args{ACTION}' not in {'List', 'Add', 'Remove'}."); 
}













sub libraries_list
{
    # ---- Load the library database.
    my $libs_db_fn = "$CACHE_DIR/libraries.csv";
    my $libs_db = new_naif_db()->from_csv_or_die($libs_db_fn, $libs_db_field_validator);
    my $n_db_libs = $libs_db->size();
    print_action("$n_db_libs cached libraries loaded.");

    # ---- Fetch only the library records that match the various criteria.
    
    my $libs_tmp_db = $libs_db
        ->sub_db_select("library_name", $args{LIBRARIES})
        ->sub_db_select("project_name", $args{PROJECT})
        ->sub_db_select("organism_name", $args{ORGANISM});


    # ---- Select libraries with the specified type.
    
    my $libs_sub_db;
    if ($args{TYPE} eq "frag") {
        $libs_sub_db = $libs_tmp_db->sub_db_select_test("insert_size", sub { return ! $_[0]; });
    }
    elsif ($args{TYPE} eq "jump") {
        $libs_sub_db = $libs_tmp_db->sub_db_select_test("insert_size", sub { return $_[0] && $_[0] > 0; });
    }
    else {
        $libs_sub_db = $libs_tmp_db;
    }
    

    # ---- Print libraries to stdout ("-")
    if ($libs_sub_db->size() > 0) {
        $libs_sub_db->to_csv("-", $libs_db_field_validator);
    }
    else {
        print_error("No library matches the specified parameters.");
    }
    print "\n";


}



















sub libraries_add
{
    my $libs_db_fn = "$CACHE_DIR/libraries.csv";

    # ---- Open the input lib database from .csv file.
    my $in_libs_db_fn = $args{IN_LIBS_CSV};
    abort("Unspecified parameter IN_LIBS_CSV.") 
        unless ($in_libs_db_fn);

    my $in_libs_db = new_naif_db()->from_csv_or_die($in_libs_db_fn, $in_libs_db_field_validator);

    abort(Tag() . "There are repeated library names in your '$in_libs_db_fn'.")
        unless $in_libs_db->is_key("library_name");

    my $n_in_libs = $in_libs_db->size();
    print_action("$n_in_libs input libraries loaded from '$in_libs_db_fn':");
    $in_libs_db->to_csv("-", $in_libs_db_fields, "*"); 


    # ---- Merging libraries
    print_action(Tag() . "Merging input libraries with cache libraries.");

    my @fields_report = ("library", 
                         "organism",
                         "project",
                         "type", 
                         "insert_size",
                         "frag_size",
                         "pairing",
                         "genomic_range",                         
                         "status", "action");
    my $db_report = new_naif_db(\@fields_report);


    my $n_new_libs = 0;
    my $n_updated_libs = 0;
    for (my $i = 0; $i != $n_in_libs; $i++) 
    {
        my $rec_in_lib = $in_libs_db->record_get($i);

        my $rec_report = $db_report->record_new();
        $rec_report->{library}       = $rec_in_lib->{library_name};
        $rec_report->{organism}      = $rec_in_lib->{organism_name};
        $rec_report->{project}       = $rec_in_lib->{project_name};
        $rec_report->{type}          = $rec_in_lib->{type};
        $rec_report->{insert_size}   = ($rec_in_lib->{insert_size}) ?
            sprintf("%5s +- %4s", $rec_in_lib->{insert_size}, $rec_in_lib->{insert_stddev}) :
            "";
        $rec_report->{frag_size}     = ($rec_in_lib->{frag_size}) ?
            sprintf("%4s +- %3s", $rec_in_lib->{frag_size}, $rec_in_lib->{frag_stddev}) :
            "";
        $rec_report->{pairing}       = $rec_in_lib->{read_orientation};
        $rec_report->{genomic_range} = sprintf("[ %3s : %3s ]", 
                                               $rec_in_lib->{genomic_start},
                                               $rec_in_lib->{genomic_end});


        if (!$args{DRY_RUN}) {
            # $outcome -> { new, equal, updated }
            my $outcome = cache_db_rec_add($libs_db_fn, 
                                           $rec_in_lib, 
                                           "library",
                                           $args{OVERWRITE});

            if ($outcome->{new}) {
                $rec_report->{status} = "NEW";
                $rec_report->{action} = "ADD";
            }
            else {
                $rec_report->{status} = "OLD";
                $rec_report->{status} .= ($outcome->{equal} ? " EQUAL" : " DIFFERENT")
                    if (defined $outcome->{equal});
                $rec_report->{action} = ($outcome->{updated}) ? "OVERWRITE" : "SKIP";
            }

            if ($outcome->{new})     { $n_new_libs++; }
            if ($outcome->{updated}) { $n_updated_libs++; }
        }
        else {
            $rec_report->{status} = "DRY RUN";
            $rec_report->{action} = "SKIP";
        }

        $db_report->record_add($rec_report);
    }


    $db_report->to_csv("-", \@fields_report, "=");

    
    print_action(Tag() . "$n_new_libs new libraries added to library database.");
    print_action(Tag() . "$n_updated_libs libraries updated in library database.");
}

















sub libraries_remove
{
    # ---- Validate LIBRARIES.
    abort(Tag() . "Must specify LIBRARIES.")
        unless ($args{LIBRARIES} ne "");


    # ---- Load the library database.
    my $libs_db_fn = "$CACHE_DIR/libraries.csv";
    my $libs_db = new_naif_db()->from_csv_or_continue($libs_db_fn, $libs_db_field_validator);
    my $n_db_libs = $libs_db->size();
    print_action("$n_db_libs cached libraries loaded.");

    # ---- Remove libraries.
    
    my $n_removed = 0;
    my @names = array_from_ref_or_value($args{LIBRARIES});
    foreach my $name (@names) 
    {
        my $removed = cache_db_rec_erase($libs_db_fn, "library", $name);

        if ($removed) {
            print_action("Removed library '$name'.");
            $n_removed++;
        }
        else {
            print_error("Didn't find library '$name' to remove.");
        }
    }    
    
    print_subaction("$n_removed libraries removed.");
    return $n_removed;

}

