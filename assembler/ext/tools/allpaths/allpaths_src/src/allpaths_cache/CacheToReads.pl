#!/usr/bin/perl -w
###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2011) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################
#
# CacheToReads.pl
#
# Take a set of cached fastb/qualb files and merge them into a single set of  
# reads, ready for the ALLPATHS pipeline.
#
# 2010-05   Filipe Ribeiro
#


use strict;
use FindBin;

# Local libraries, from elsewhere in the BroadCRD repository.
use lib "$FindBin::Bin/../";
use lib "$FindBin::Bin/";

use PerlRunTime;    # run_or_die
use PerlUtils;      # ISO_date, abort
use ArachneArgs;    # Command-line ARG=VALUE parser

use NaifDB;         # database handler
use AllPathsCache;  # definitions of database fields


sub Tag { return ISO_date() . " (CTR): "; }


my $LONG_JUMP_SIZE = 20000;


# CONTROL BEGINS HERE
# Parse command-line options of the form KEY=value.
# This function comes from the ArachneArgs.pm module.

my %args = getCommandArguments
    (CACHE_DIR            => { value => "",
                               help  => "Directory where reads in AllPaths format are cached."},
     GROUPS               => { value => "",
                               help  => "Groups names of reads to import." },
     FRACTIONS            => { value => 1,
                               help  => "Optional fraction of all reads to import." },
     OUT_HEAD             => { value => undef,
                               help  => "The prefix for the input files to be generated." },
     DRY_RUN              => { value => "False",
                               help  => "If 'True' only say what would be done." },
     VERBOSE              => { value => "True" });


set_verbose($args{VERBOSE});


# ---- Make sure cache dir is defined. 

my $CACHE_DIR = $args{CACHE_DIR};
if ($CACHE_DIR eq "") 
{
    abort(Tag() . "You didn't specify CACHE_DIR nor set the environment variable \$ALLPATHS_CACHE_DIR.")
        unless defined $ENV{ALLPATHS_CACHE_DIR};
    $CACHE_DIR = $ENV{ALLPATHS_CACHE_DIR};
    abort(Tag() . "Specified CACHE_DIR '$CACHE_DIR' is not an absolute path.") 
        unless ($CACHE_DIR =~ /^\//);
}
abort(Tag() . "Can't find cache directory '$CACHE_DIR'.") 
    unless (-d $CACHE_DIR);


# ---- Make sure the path in OUT_HEAD is defined. 

abort(Tag() . "the path to OUT_HEAD '$args{OUT_HEAD}' is not a directory.") 
    if ($args{OUT_HEAD} =~ /^(.+)\/[^\/]+$/ && ! -d $1);



# ---- Open the cache databases from their .csv files.

my $groups_db_fn = "$CACHE_DIR/groups.csv";
#print_action(Tag() . "Loading cached groups '$groups_db_fn'.");
my $groups_db = new_naif_db()->from_csv_or_continue($groups_db_fn, $groups_db_field_validator);
my $n_db_groups = $groups_db->size();
#print_subaction("$n_db_groups groups already in the cache.");


my $libs_db_fn = "$CACHE_DIR/libraries.csv";
#print_action(Tag() . "Loading cached libraries '$libs_db_fn'.");
my $libs_db = new_naif_db()->from_csv_or_continue($libs_db_fn, $libs_db_field_validator);
my $n_db_libs = $libs_db->size();
#print_subaction("$n_db_libs libraries already in the cache.");





#print_action(Tag() . "Selecting groups to import.");


# ---- Build file records for all the groups 

my @arg_group_names = array_from_ref_or_value($args{GROUPS}); 
my $n_groups = scalar @arg_group_names;


my @arg_group_fractions = array_from_ref_or_value($args{FRACTIONS}); 

if (@arg_group_fractions == 1) {
    @arg_group_fractions = map $arg_group_fractions[0], (1..$n_groups);
} 
elsif (@arg_group_fractions < $n_groups) {
    abort(Tag() . "Number of FRACTIONS provided does not agree with number of GROUPS.");
}
my %fraction_from_group = ();
map $fraction_from_group{$arg_group_names[$_]} = $arg_group_fractions[$_], (0..$n_groups-1);


# ---- NOTE: the order in the sub_db is unpredictable!
my $groups_sub_db = $groups_db->sub_db_select("group_name", @arg_group_names);


my @group_names      = @{$groups_sub_db->{group_name}};
my @group_lib_names  = @{$groups_sub_db->{library_name}};
my @group_fractions  = map $fraction_from_group{$_}, @group_names;
my @group_read_sizes = @{$groups_sub_db->{mean_read_len}};

#print "groups: @group_names\n";
#print "libs: @group_lib_names\n";
#print "frac: @group_fractions\n";

# -- We want to sort over lib names to bring reads from the same library together
 
my @is = sort { $group_lib_names[$a] cmp $group_lib_names[$b] } (0..$n_groups-1);

# -- For each group gather some stats

#my @sorted_group_names      = map $group_names[$_], @is;
#my @sorted_group_lib_names  = map $group_lib_names[$_], @is;
#my @sorted_group_fractions  = map $group_fractions[$_], @is;
#my @sorted_group_read_sizes = map $group_read_sizes[$_], @is;
my @sorted_group_names      = @group_names[@is];
my @sorted_group_lib_names  = @group_lib_names[@is];
my @sorted_group_fractions  = @group_fractions[@is];
my @sorted_group_read_sizes = @group_read_sizes[@is];

#print "groups: @sorted_group_names\n";
#print "libs: @sorted_group_lib_names\n";
#print "frac: @sorted_group_fractions\n";


my @sorted_group_paired = ();
my @sorted_group_sizes = ();
my @sorted_group_stddevs = ();

my %organisms;
for (my $i = 0; $i != $n_groups; $i++) {
    my $group_lib_name = $sorted_group_lib_names[$i];
    my $library_rec = $libs_db->sub_db_select("library_name", $group_lib_name)->record_get(0);
    
    $organisms{$library_rec->{organism_name}}++;

    push @sorted_group_paired, $library_rec->{paired};
    if ($library_rec->{insert_size}) { # jumping library
        push @sorted_group_sizes,   int(0.5 + $library_rec->{insert_size});
        push @sorted_group_stddevs, int(0.5 + $library_rec->{insert_stddev});
    }
    elsif ($library_rec->{frag_size}) { # fragment library; insert size for fragment libraries should not be defined 
        push @sorted_group_sizes,   int(0.5 + $library_rec->{frag_size});
        push @sorted_group_stddevs, int(0.5 + $library_rec->{frag_stddev});
    }
    else { # just set them at zero
        push @sorted_group_sizes, 0;
        push @sorted_group_stddevs, 0;
    }
}

# ---- check if all groups are from the same organism (spelling sensitive)

if (scalar keys %organisms > 1) {
    print_warning("Trying to merge reads from different organisms!");
    foreach my $org (keys %organisms) {
        print_sub_action("$organisms{$org} $org");
    }
}




# ---- Build command-line arguments for CacheReadsMerge



my $group_heads_arg       = globable_str_from_array_or_value(@sorted_group_names);
my $group_fracs_arg       = globable_str_from_array_or_value(@sorted_group_fractions);
my $group_lib_names_arg   = globable_str_from_array_or_value(@sorted_group_lib_names);
my $group_lib_paired_arg  = globable_str_from_array_or_value(@sorted_group_paired);
my $group_lib_sizes_arg   = globable_str_from_array_or_value(@sorted_group_sizes);
my $group_lib_stddevs_arg = globable_str_from_array_or_value(@sorted_group_stddevs);
my $group_read_sizes_arg  = globable_str_from_array_or_value(@sorted_group_read_sizes);





# ---- Actually do the merging.

#print_action(Tag() . "Merging reads into AllPaths format.");

my $cmd = 
    "$FindBin::Bin/CacheReadsMerge" .
    " CACHE_DIR=$CACHE_DIR" .
    " GROUP_HEADS='$group_heads_arg'" .
    " GROUP_LIB_NAMES='$group_lib_names_arg'" .
    " GROUP_LIB_PAIRED='$group_lib_paired_arg'" .
    " GROUP_LIB_SIZES='$group_lib_sizes_arg'" .
    " GROUP_LIB_STDDEVS='$group_lib_stddevs_arg'" .
    " GROUP_READ_SIZES='$group_read_sizes_arg'" .
    " GROUP_FRACS='$group_fracs_arg'" .
    " OUT_HEAD=$args{OUT_HEAD}" .
    " NH=" . bool_str(!$args{VERBOSE});


if ($args{DRY_RUN}) {
    print "DRY_RUN: $cmd\n";
}
else {
    if ($args{VERBOSE}) {
        run_or_die($cmd); # stdout and stderr are kept
    }
    else {
        my $output = run_or_die($cmd); # stdout and stderr go into $output
    } 
}











