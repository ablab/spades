#!/usr/bin/perl -w
###########################################################################
#       This software and its documentation are copyright (2009) by the   #
#   Broad Institute/Massachusetts Institute of Technology.  All rights    #
#   are reserved.  This software is supplied without any warranty or      #
#   guaranteed support whatsoever. Neither the Broad Institute nor MIT    #
#   can be responsible for its use, misuse, or functionality.             #
###########################################################################
#
# OldPicardToCache.pl
#
# Temporary tool to convert old Picard cache to the new cache.
#
# 2010-06   Filipe Ribeiro     ribeiro@broadinstitute.org
#


use strict;
use FindBin;

# Local libraries, from elsewhere in the BroadCRD repository.

use lib "$FindBin::Bin/../";
use lib "$FindBin::Bin/";

use PerlRunTime; # run_or_die
use PerlUtils;   # ISO_date, abort
use ArachneArgs; # Command-line ARG=VALUE parser

use NaifDB;      # Filipe's database handler
use AllPathsCache;  # definitions of database fields


sub Tag { return ISO_date() . " (OPTC): "; }




# ---- CONTROL BEGINS HERE
#      Parse command-line options of the form KEY=value.
#      This function comes from the ArachneArgs.pm module.

my %args = getCommandArguments(PICARD_DIR           => undef,
                               CACHE_DIR            => "",
                               DRY_RUN              => 1,
                               VERBOSE              => 1);
                               

my $PICARD_DIR = $args{PICARD_DIR};
abort(Tag() . "**** '$PICARD_DIR' doesn't exist.") 
    unless (-d $PICARD_DIR);



# ---- Make sure cache dir is defined 
my $CACHE_DIR = $args{CACHE_DIR};
if ($CACHE_DIR eq "") 
{
    abort(Tag() . "**** You didn't specify CACHE_DIR nor set the environment variable \$ALLPATHS_CACHE_DIR.")
        unless defined $ENV{ALLPATHS_CACHE_DIR};
    $CACHE_DIR = $ENV{ALLPATHS_CACHE_DIR};
}






# ---- Load the libraries database

my $libs_db_fn = "$CACHE_DIR/libraries.csv";
print Tag() . "Loading cache libraries '$libs_db_fn'.\n";

my $libs_db = new_naif_db()->from_csv_or_continue($libs_db_fn, $libs_db_field_validator);
my $n_libs = $libs_db->size();
printf("     %d libraries already in the cache.\n", $n_libs);





# ---- Load the groups database

my $groups_db_fn = "$CACHE_DIR/groups.csv";
print Tag() . "Loading cache groups '$groups_db_fn'.\n";

my $groups_db = new_naif_db()->from_csv_or_continue($groups_db_fn, $groups_db_field_validator);







# ---- load FIR.libraries.txt

my $FIR_fn = "FIR.libraries.txt";
print Tag() . "Loading old Picard library data in '$FIR_fn'.\n";


my $FIR_libs = read_FIR_libs($FIR_fn);








# ---- go through all the Picard directories and convert to cache
print Tag() . "Converting groups.\n";


my @dns = glob "$PICARD_DIR/?????.?";

foreach my $dn (@dns) {

    print "-" x 80 . "\n";
    print "'$dn'\n\n";

    
    
    my ($in_lib_rec, $in_group_rec) = 
        rec_import_from_picard($FIR_libs, $dn);


    my $n_libs_added = libraries_db_add_rec($libs_db, $in_lib_rec);
    if ($n_libs_added) {
        $libs_db->to_csv($libs_db_fn);
    }


    my $group_name = $in_group_rec->{group_name};
    
    foreach my $ext ("fastb", "qualb") {
        my $in_fn = "$dn/reads.$ext";
        my $out_fn = "$CACHE_DIR/$group_name.$ext";
        
        print "     copying '$in_fn' to '$out_fn'.\n";
        symlink $in_fn, $out_fn  unless -e $out_fn;
    }



    print "     adding group '$group_name' to groups database.\n";
    my $added = groups_db_add_rec($groups_db, 
                                 $in_group_rec, 
                                 "$CACHE_DIR/$group_name.fastb");
    
    if ($added) {
        $groups_db->to_csv($groups_db_fn);
    }
}



















# ----------------------
# Subroutine definitions
# ----------------------



sub rec_import_from_picard
{
    my ($picard_libs, $picard_dn) = @_;

    my $params = read_params_txt("$picard_dn/params.txt");

    my $lib_name = $params->{LIBRARY_NAME};

    my ($size, $sd) = (exists $picard_libs->{$lib_name}) ? @{$picard_libs->{$lib_name}} : (0, 0);

    my $group_name = "$params->{FLOWCELL_BARCODE}.$params->{LANE}";

    my %in_group_rec = ( library_name => $params->{LIBRARY_NAME},
                         group_name   => $group_name,
                         file_name    => find_bam($params->{ANALYSIS_DIR}, $group_name) );


    my ($orientation, $type, $i_size, $i_sd, $f_size, $f_sd) = ("", "", "", "", "", "");
    if ($size == 0) {
        $type = "unknown";
        $orientation = "unknown";
    }
    elsif ($size < 500) {
        $type = "fragment";
        $orientation = "inward";
        $f_size = $size;
        $f_sd = $sd;
    } 
    else {
        $type = "jumping";
        $orientation = "outward";
        $i_size = $size;
        $i_sd = $sd;
    }

    my $project_name = exists $params->{PROJECT} ? $params->{PROJECT} : "";

    my %in_lib_rec = ( library_name     => $params->{LIBRARY_NAME},
                       project_name     => $project_name,
                       organism_name    => $params->{ORGANISM},
                       paired           => 1,
                       frag_size        => $f_size,
                       frag_stddev      => $f_sd,
                       insert_size      => $i_size,
                       insert_stddev    => $i_sd,
                       type             => $type,
                       read_orientation => $orientation,
                       genomic_start    => "",
                       genomic_end      => "" );


    return (\%in_lib_rec, \%in_group_rec);
}



sub find_bam
{
    my ($dn, $gn) = @_;
    $dn =~ s/\/$//;
    foreach my $fn (glob "$dn/*bam") {
        return $fn if ($fn eq "$dn/$gn.unmapped.bam");
        return $fn if ($fn eq "$dn/$gn.aligned.duplicates_marked.bam");
    }
    return "$dn/$gn.....bam";
}



sub read_params_txt
{
    my ($fn) = @_;

    open F, "<$fn" or abort("**** Can't open '$fn'.");
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




sub read_FIR_libs
{
    my ($FIR_fn) = @_;

    $FIR_fn = "$PICARD_DIR/$FIR_fn" if (-e "$PICARD_DIR/$FIR_fn");
    $FIR_fn = "$PICARD_DIR/../$FIR_fn" if (-e "$PICARD_DIR/../$FIR_fn");


    open FILE, "<$FIR_fn" or abort("**** Can't open '$FIR_fn'");
    my %picard_libs;
    while (<FILE>) 
    {
        if (/^\s*([^\#]\S+)\s+(\d+)\s+(\d+)/) {
            $picard_libs{$1} = [$2, $3];
            print "$1 - $2 - $3\n";
        }
    }
    close FILE;
    return \%picard_libs;
}














































