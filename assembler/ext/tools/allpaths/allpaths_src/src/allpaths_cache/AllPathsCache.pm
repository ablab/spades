###########################################################################
#       This software and its documentation are copyright (2009) by the   #
#   Broad Institute/Massachusetts Institute of Technology.  All rights    #
#   are reserved.  This software is supplied without any warranty or      #
#   guaranteed support whatsoever. Neither the Broad Institute nor MIT    #
#   can be responsible for its use, misuse, or functionality.             #
###########################################################################
# 
# AllPathsCache.pm
#
# 2010-06-16    Filipe Ribeiro     ribeiro@broadinstitute.org
#

package AllPathsCache;

use Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(
$groups_db_field_validator
$in_groups_db_field_validator
$libs_db_field_validator
$in_libs_db_field_validator

$groups_db_fields
$in_groups_db_fields
$libs_db_fields
$in_libs_db_fields


&get_directory_from_arg_or_env
&directory_is_absolute_or_die

&cache_db_rec_add
&cache_db_rec_erase
);

use strict;   
use IO::File;       # ->new()

use PerlRunTime;    # run_or_die
use PerlUtils;      # ISO_date, abort

use NaifDB;         # Filipe's database handler



# ---- Define fields for the input file database

$AllPathsCache::in_groups_db_field_validator = {group_name       => "safe_string_key",
						file_name        => "string",
						library_name     => "string"};

$AllPathsCache::in_groups_db_fields = ["group_name",
				       "library_name",
				       "file_name"];


# ---- Define fields for the cache file database

$AllPathsCache::groups_db_field_validator = {group_id         => "number_key",
					     read_count       => "number",
					     min_read_len     => "number",
					     max_read_len     => "number",
					     mean_read_len    => "number",
					     mean_r1_qual     => "number",
					     mean_r2_qual     => "number",
					     import_date      => "string"};

#      Import the fields from the input db
map $AllPathsCache::groups_db_field_validator->{$_} = $AllPathsCache::in_groups_db_field_validator->{$_}, 
    keys %$AllPathsCache::in_groups_db_field_validator;


$AllPathsCache::groups_db_fields = ["group_name",
				    "group_id",
				    "library_name",
				    "read_count",
				    "min_read_len",
				    "mean_read_len",
				    "max_read_len",
				    "mean_r1_qual",
				    "mean_r2_qual",
				    "import_date",
				    "file_name"];





# ---- Define fields for the input library database

$AllPathsCache::in_libs_db_field_validator = {library_name      => "string_key",
					      project_name      => "string",
					      organism_name     => "string",
					      type              => "string",
					      paired            => "bool",
					      frag_size         => "number",
					      frag_stddev       => "number",
					      insert_size       => "number",
					      insert_stddev     => "number",
					      read_orientation  => { inward => 1, outward => 1, "" => 1, "unknown" => 1 }, 
					      genomic_start     => "number",
					      genomic_end       => "number"};


$AllPathsCache::in_libs_db_fields = ["library_name",
				     "project_name",
				     "organism_name",
				     "type",
				     "paired",
				     "frag_size",
				     "frag_stddev",
				     "insert_size",
				     "insert_stddev",
				     "read_orientation",
				     "genomic_start",
				     "genomic_end"];



# ---- Define fields for the cache library database

$AllPathsCache::libs_db_field_validator = {library_id        => "number_key",
					   import_date       => "string"};

#      Import the fields from the input db
map $AllPathsCache::libs_db_field_validator->{$_} = $AllPathsCache::in_libs_db_field_validator->{$_}, 
    keys %$AllPathsCache::in_libs_db_field_validator;


$AllPathsCache::libs_db_fields = ["library_name",
				  "library_id",
				  "project_name",
				  "organism_name",
				  "type",
				  "paired",
				  "frag_size",
				  "frag_stddev",
				  "insert_size",
				  "insert_stddev",
				  "read_orientation",
				  "genomic_start",
				  "genomic_end",
				  "import_date"];






sub get_directory_from_arg_or_env
{
    my ($args, $arg_name) = @_;
    my $dir = $args->{$arg_name};
    if ($dir eq "") 
    {
        abort("You didn't specify $arg_name nor set the environment variable \$ALLPATHS_$arg_name.")
            unless defined $ENV{"ALLPATHS_$arg_name"};
        $dir = $ENV{"ALLPATHS_$arg_name"};
    }
    abort("Directory '$dir' does not exist.")
        unless (-d $dir);
    return $dir;
}


sub directory_is_absolute_or_die
{
    my ($dir) = @_;
    abort("Directory '$dir' is not an absolute path.")
        unless ($dir =~ /^\//);
}












sub cache_db_rec_add
{
    my ($db_fn, $in_rec, $type, $overwrite) = @_;

    my $db_field_validator;
    my $db_fields;
    if    ($type eq "library") { 
	$db_field_validator = $::libs_db_field_validator;
	$db_fields          = $::libs_db_fields;
    }
    elsif ($type eq "group")   { 
	$db_field_validator = $::groups_db_field_validator;
	$db_fields          = $::groups_db_fields;
    }
    else  { abort("Invalid database type '$type'."); }

    # ---- Load the database.

    #print_action("Loading $type database '$db_fn'.");
    my $db = new_naif_db($db_field_validator);
    my $lock_fh = $db->from_csv_lock_exclusive($db_fn, $db_field_validator);

    my $n = $db->size();
    #print_subaction("$n entries in $type database.");
    #print "\n";


    # ---- Count number of matches in db
    my $name = $in_rec->{$type."_name"};

    my $n_matches = $db->sub_db_select_count($type."_name", $name);
    if ($n_matches > 1) {
        abort("Entry '$name' matches more than one record in $type database.");
    }

    

    my $outcome = { new     => undef,
                    equal   => undef,
                    updated => undef };


    if ($n_matches == 0)         # ---- New entry!
    {
        $outcome->{new}     = 1;
        $outcome->{updated} = 1;
    }
    else                         # ---- Entry exists already!
    {
        $outcome->{new}     = 0;
        $outcome->{updated} = (defined $overwrite ? $overwrite : 0);

        # ---- Check if entry is identical or different 

        my $rec = $db->sub_db_select($type."_name", $name)->record_get(0);

        $outcome->{equal} = 1;
        foreach my $key (keys %$in_rec) {
            if ($in_rec->{$key} ne $rec->{$key}) {
                printf("     %-20s -> %-20s  !=  %-20s (cache)\n", 
                       $key, $in_rec->{$key}, $rec->{$key});
                $outcome->{equal} = 0;
            }
        }
    }


    

    # ---- Update the database in case of new or overwrite
    
    if (!$outcome->{updated})   # didn't change -> close file.
    {
        close($lock_fh);   
    }
    else                        # changed -> write it out.
    {
        my $rec = {%$in_rec};
        $rec->{import_date} = ISO_date();

        if ($outcome->{new})   # new entry
        {
            $rec->{$type."_id"} = $db->next_4d_index($type."_id");
            $db->record_add($rec);
        }
        else                   # replacement entry 
        {
            my $i_rec = $db->record_index($type."_name", $name);
            $rec->{$type."_id"} = $db->{$type."_id"}[$i_rec];
            $db->record_set($i_rec, $rec); 
        }
        
        $db->to_csv_lock_exclusive($lock_fh, $db_fields);
    }


    # ---- Issue report.

    if (0) {
        $db->sub_db_select($type."_name", $name)->to_csv("-", $db_fields, "-#-");
        
        print "\n";
        if (!$outcome->{new}) {
            print_subaction("$type '$name' EXISTS and is " .
                            ($outcome->{equal} ? "IDENTICAL to" : "DIFFERENT from") .
                            " the one in $type database.");
        }
        if ($outcome->{updated}) {
            print_subaction("$type '$name' is " .
                            (($outcome->{new}) ? 
                             "NEW and was ADDED to" : 
                             "NOT NEW and REPLACED the one in") .
                            " $type database.");
        }
        print "\n";
    }
    return $outcome;
}




sub cache_db_rec_erase
{
    my ($db_fn, $type, $name) = @_;

    my $db_field_validator;
    my $db_fields;
    if    ($type eq "library") { 
	$db_field_validator = $::libs_db_field_validator;
	$db_fields          = $::libs_db_fields;
    }
    elsif ($type eq "group")   { 
	$db_field_validator = $::groups_db_field_validator;
	$db_fields          = $::groups_db_fields;
    }
    else  { abort("Invalid database type '$type'."); }

    # ---- Load the database.

    my $db = new_naif_db();
    my $lock_fh = $db->from_csv_lock_exclusive($db_fn, $db_field_validator);

    # ---- Count number of matches in the libs db
    my $n_matches = $db->sub_db_select_count($type."_name", $name);
    if ($n_matches > 1) {
        abort("Entry '$name' matches more than one record in $type database.");
    }

    if ($n_matches == 1) 
    {
        # ---- Update the database in case of new or overwrite
        my $i_rec = $db->record_index($type."_name", $name);
        $db->record_erase($i_rec);
        $db->to_csv_lock_exclusive($lock_fh, $db_fields);
    }
    close($lock_fh);
    
    return $n_matches;
}

























1;
