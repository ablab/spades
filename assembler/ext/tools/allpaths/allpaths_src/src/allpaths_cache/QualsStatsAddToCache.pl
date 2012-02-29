#!/usr/bin/perl -w
###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2011) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################
#
#


use strict;
use FindBin;
use POSIX ":sys_wait_h";     # for the WNOHANG signal to wait pid

# Local libraries, from elsewhere in the BroadCRD repository.
use lib "$FindBin::Bin/../";
use lib "$FindBin::Bin/";

use PerlRunTime;    # run_or_die
use PerlUtils;      # ISO_date, abort
use ArachneArgs;    # Command-line ARG=VALUE parser

use NaifDB;         # database handler
use AllPathsCache;  # definitions of database fields


sub Tag { return ISO_date() . " (CG): "; }

my %args = getCommandArguments
    (CACHE_DIR            => { value => undef,
                               help  => "Directory where reads in AllPaths format are cached."});


my $groups_db_fn = "$args{CACHE_DIR}/groups.csv";
my $groups_db = new_naif_db()->from_csv_or_continue($groups_db_fn);
my $n_db_groups = $groups_db->size();
print_action("$n_db_groups cached groups loaded.");


my @new_fields = ($groups_db->fields(), "mean_r1_qual", "mean_r2_qual");

my $new_groups_db = new_naif_db(\@new_fields);

for (my $i = 0; $i != $n_db_groups; $i++)
{    
    my $group_rec = $groups_db->record_get($i);
    my $group_name = $group_rec->{group_name};
    
    print "group = $group_name\n";

    my $stats_fn = "$args{CACHE_DIR}/qualb_stats/$group_name.out";

    $group_rec->{mean_r1_qual} = -1;
    $group_rec->{mean_r2_qual} = -1;
    if (-e $stats_fn) {    
	open FILE, $stats_fn;
	while (<FILE>) {
	    
	    if (/^\s+(\S+)\s+(\S+)\s+mean Q\s*$/) {
		$group_rec->{mean_r1_qual} = $1;
		$group_rec->{mean_r2_qual} = $2;
	    }
	}
	close FILE;
    }

    

    $new_groups_db->record_add($group_rec);
}


$new_groups_db->to_csv($groups_db_fn);

