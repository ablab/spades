#!/usr/bin/perl -w
###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2011) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################
#
# CacheToAllPathsInputs.pl
#
# Import reads from the cache and setup an AllPaths run.
#
# 2010-07   Filipe Ribeiro
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


sub Tag { return ISO_date() . " (CTAPI): "; }



# CONTROL BEGINS HERE
# Parse command-line options of the form KEY=value.
# This function comes from the ArachneArgs.pm module.

my %args = getCommandArguments
    (CACHE_DIR            => { value => "",
                               help  => "Directory where reads in AllPaths format are cached."},
     DATA_DIR             => { value => undef,
                               help  => "full path to the ALLPATHS DATA directory." },
     PLOIDY               => { value => "",
                               help  => "The ploidy." },
     IN_GROUPS_CSV        => { value => "",
                               help  => "Comma-separated-value info on read groups to import." }, 
     GROUPS               => { value => "",
                               help  => "Groups names of reads to import." },
     GENOME_SIZE          => { value => "",
                               help  => "Optional estimated total genome size." },
     COVERAGES            => { value => "",
                               help  => "Group by group desired coverage (must specify GENOME_SIZE)." },
     FRACTIONS            => { value => "",
                               help  => "Group by group desired read_fraction." },
     FRAG_COVERAGE        => { value => "",
                               help  => "Fragment library desired coverage (must specify GENOME_SIZE)." },
     JUMP_COVERAGE        => { value => "",
                               help  => "Jumping library desired coverage (must specify GENOME_SIZE)." },
     LONG_JUMP_COVERAGE   => { value => "",
                               help  => "Long jumps (insert size > 10 kb) library desired coverage (must specify GENOME_SIZE)." },
     LONG_JUMP_MIN_SIZE   => { value => 20000,
                               help  => "Threshold for long jumps." },
     LONG_READ_MIN_LEN    => { value => 500,
                               help  => "Threshold for long reads." },
     FRAG_FRAC            => { value => "",
                               help  => "Fragment library desired read fraction." },
     JUMP_FRAC            => { value => "",
                               help  => "Jumping library desired read fraction." },
     LONG_JUMP_FRAC       => { value => "",
                               help  => "Long jumps library desired read fraction." },
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



# ---- Some argument error checking 

my $DATA_DIR = $args{DATA_DIR};
mkdir $DATA_DIR or abort(Tag() . "Can't create '$DATA_DIR'.")
    unless -d $DATA_DIR;


abort(Tag() . "Can't specify both IN_GROUPS_CSV and GROUPS.") 
    if ($args{IN_GROUPS_CSV} && $args{GROUPS});

abort(Tag() . "Must specify one of IN_GROUPS_CSV and GROUPS.") 
    if (!$args{IN_GROUPS_CSV} && !$args{GROUPS});

abort(Tag() . "Can't specify both GROUP_FRAC and other *_FRAC.")
    if ($args{FRACTIONS} && ($args{FRAG_FRAC} ne "" || 
                             $args{JUMP_FRAC} ne "" || 
                             $args{LONG_JUMP_FRAC} ne ""));

abort(Tag() . "Can't specify both GROUP_COVERAGE and other *_COVERAGE.")
    if ($args{COVERAGES} && ($args{FRAG_COVERAGE} ne "" || 
                             $args{JUMP_COVERAGE} ne "" || 
                             $args{LONG_JUMP_COVERAGE} ne ""));

abort(Tag() . "Must also specify GENOME_SIZE when specifying a COVERAGE.")
    if (($args{COVERAGES} || 
         $args{FRAG_COVERAGE} ne "" || 
         $args{JUMP_COVERAGE} ne "" || 
         $args{LONG_JUMP_COVERAGE} ne "") &&
         !$args{GENOME_SIZE});




# ---- Open the cache databases from their .csv files.

my $groups_db_fn = "$CACHE_DIR/groups.csv";
#print_action(Tag() . "Loading cached groups '$groups_db_fn'.");
my $groups_db = new_naif_db()->from_csv_or_continue($groups_db_fn, $groups_db_field_validator);
my $n_db_groups = $groups_db->size();
print_action(Tag() . "$n_db_groups groups loaded from the cache.");




my $libs_db_fn = "$CACHE_DIR/libraries.csv";
#print_action(Tag() . "Loading cached libraries '$libs_db_fn'.");
my $libs_db = new_naif_db()->from_csv_or_continue($libs_db_fn, $libs_db_field_validator);
my $n_db_libs = $libs_db->size();
print_action(Tag() . "$n_db_libs libraries loaded from the cache.");



# ---- Create database for reporting purposes
my @fields_report =     ("group", 
                         "library", 
                         "type", 
                         "prefix_out",
                         "reads_in",
                         "fraction",
                         "reads_out",
                         "coverage",
                         "phys_cov",
                         "pairing",
                         "library_size",
                         "overlap",
                         "min_mean_max_len",
                         "mean_Qs");
my $db_report = new_naif_db(\@fields_report);

my @fields_report_all = ("prefix_out",
                         "n_groups", 
                         "n_libraries", 
                         "reads_in",
                         "fraction",
                         "reads_out",
                         "coverage",
                         "phys_cov",
                         "overlap",
                         "min_mean_max_len",
                         "mean_Qs");
my $db_report_all = new_naif_db(\@fields_report);


# ---- Compile the various libraries
my $module = "CacheToReads.pl";
my @not_found = ();
foreach my $type ("frag", "jump", "long_jump", "long") 
{
    my ($groups_arg, $fractions_arg) = 
        build_args_for_type(\%args, $groups_db, $libs_db, $type, $db_report, $db_report_all);

    if ($groups_arg ne "") {    
        #print_action(Tag() . "Merge cached '$type' reads into ALLPATHS input files.");
        system("$FindBin::Bin/$module" .
               " CACHE_DIR=$CACHE_DIR" .
               " GROUPS='$groups_arg'" .
               " FRACTIONS='$fractions_arg'" .
               " OUT_HEAD=$DATA_DIR/". $type . "_reads_orig" .
               " DRY_RUN=$args{DRY_RUN}" .
               " VERBOSE=$args{VERBOSE}");
        print "\n";
    }
    else {
        push @not_found, $type;
    }
}

# ---- Print report

$db_report->to_csv("-", \@fields_report, "=");
$db_report_all->to_csv("-", \@fields_report_all, "=");




if (@not_found) {
    print "\n";
    print "==================== WARNINGS ====================\n";
    foreach my $type (@not_found) {
        print "\n";
        print("!!!! No '$type' cached read groups found.\n");
        if ($type eq "frag") {
            print("     YOU CAN'T RUN AND ASSEMBLY WITHOUT FRAGMENT READS!\n");
            print("     Remember, fragment libraries must have empty 'insert_size' and 'insert_stddev' in the 'in_libs.csv'.\n");
        }
        elsif ($type eq "jump") {
            print("     YOU CAN'T RUN AND ASSEMBLY WITHOUT JUMPING READS!\n");
            print("     Remember, jumping libraries must have 'insert_size' and 'insert_stddev' defined in the 'in_libs.csv'.\n");
        }
        elsif ($type eq "long_jump") {
            print("     Long jumping reads (typically 40 kb, < 1x coverage) are useful only for scaffolding of vertebrate size genomes, and are not required for an assembly.\n");
        }
        elsif ($type eq "long") {
            print("     Long reads (typically, unpaired, 1 kb, 50x coverage) are useful only for gap patching of relatively small genomes, and are not required for an assembly.\n");
        }
        print "\n";
    }
    print "==================================================\n\n";
}



# ---- Setup Ploidy

if ($args{PLOIDY}) {
    print_action(Tag() . "Creating '$DATA_DIR/ploidy' with PLOIDY = $args{PLOIDY}.");
    system "echo $args{PLOIDY} > $DATA_DIR/ploidy";
}
elsif (-e "$DATA_DIR/ploidy") {
    my $ploidy = `cat $DATA_DIR/ploidy`;
    chomp $ploidy;
    print_action(Tag() . "'$DATA_DIR/ploidy' already exists. PLOIDY = $ploidy.");
}
else {
    print_warning(Tag() . "PLOIDY not specified and the '$DATA_DIR/ploidy' file was not created.  You must create it yourself then.");
}




# ---- Check results

my @input_fns = qw/
    frag_reads_orig.fastb
    frag_reads_orig.qualb
    frag_reads_orig.pairs
    jump_reads_orig.fastb
    jump_reads_orig.qualb
    jump_reads_orig.pairs
    /;

my $n_errors = 0;

foreach my $input_fn (@input_fns) 
{
    my $fn = "$DATA_DIR/$input_fn";
    if (!-s $fn) {
	print_error("Can't find '$fn'.  You can't run an assembly without this file.");
	$n_errors++;
    }
}
print_error(Tag() . "Found $n_errors errors.") if ($n_errors); 

print_action(Tag() . "Done.\n");







# ---- Subroutine definitions.

sub groups_validate
{
    my ($db) = @_;
    my $n_recs = $db->size();
    
    my $valid = 1;
    for (my $i = 0; $i != $n_recs; $i++) {
        my $rec = $db->record_get($i);

        # ---- swtiching this off.  jumps can have both.
	if (0 && $rec->{insert_size} && $rec->{frag_size}) {
	    $valid = 0;
	    print_warning(Tag() . "Library $rec->{library_name} can't have both 'insert_size' and 'fragment_size' specified.");
	}
    }
    abort(Tag() . "Some libraries have inconsistent entries. See above.")
	unless $valid;
}



sub is_frag
{
    my ($rec) = @_;
    return ($rec->{paired} &&
            ! $rec->{insert_size});
}

sub is_jump
{
    my ($rec) = @_;
    my $paired      = $rec->{paired};
    my $insert_size = $rec->{insert_size};
    return ($paired && 
            $insert_size && 
            $insert_size > 0 && 
            $insert_size < $args{LONG_JUMP_MIN_SIZE});
}


sub is_long_jump
{
    my ($rec) = @_;
    my $paired      = $rec->{paired};
    my $insert_size = $rec->{insert_size};
    return ($paired && 
            $insert_size && 
            $insert_size >= $args{LONG_JUMP_MIN_SIZE});
}


sub is_long
{
    my ($rec) = @_;
    my $paired        = $rec->{paired};
    my $mean_read_len = $rec->{mean_read_len};
    return (!$paired && 
            $mean_read_len && 
            $mean_read_len >= $args{LONG_READ_MIN_LEN});
}








# ---- Fetch the group names in IN_GROUPS_CSV and GROUPS.

sub build_args_for_type
{
    my ($args, $groups_db, $libs_db, $type, $db_report, $db_report_all) = @_;
    
    
    # ---- Get ALL the input groups to consider
    my $group_names = arg_group_names($args);
    my $n_groups = scalar @$group_names;

    if (@$group_names == 0) {
        print_warning(Tag() . "No groups to add.");
        return (0, 0);
    }

    # ---- get hashes of fractions or coverages 
    my ($fracs, $covs) = arg_group_fractions_coverages($args, $group_names);

    # ---- Build a db of group and lib info
    my @fields = ("group_name", 
                  "library_name", 
                  "organism_name",
                  "type",
                  "paired",
                  "frag_size",
                  "frag_stddev",
                  "insert_size", 
                  "insert_stddev", 
                  "read_count",
                  "mean_read_len",
                  "min_read_len",
                  "max_read_len");
    

    my $db = group_lib_db($group_names, $groups_db, $libs_db, \@fields);
    groups_validate($db);


    my @out_group_names = ();
    my @out_group_fracs = ();

    # ---- Select sub groups of type $type

    my $sub_db = $db->new();
    for (my $i = 0; $i != $n_groups; $i++) {
        my $rec = $db->record_get($i);
        
        if (($type eq "frag"      && is_frag($rec)) ||       # fragment library
            ($type eq "jump"      && is_jump($rec)) ||       # jumping library
            ($type eq "long_jump" && is_long_jump($rec)) ||  # long jumping library
            ($type eq "long"      && is_long($rec)))         # long unpaired library
        {
            $sub_db->record_add($rec);
            push @out_group_names, $rec->{group_name};
        }
    }
    my $n_sub_groups = $sub_db->size();


    # ---- Compute fractions (messy)

    # ----  COVERAGES  ----

    if ($args{COVERAGES}) 
    {
        for (my $i = 0; $i != $n_sub_groups; $i++) {
            my $rec = $sub_db->record_get($i);
            my $group_name = $rec->{group_name};
            my $n_bases = $rec->{read_count} * $rec->{mean_read_len};
            abort(Tag() . "No bases found for group '$group_name'.") 
                unless ($n_bases > 0);
            my $frac = $covs->{$group_name} * $args{GENOME_SIZE} / $n_bases;
            push @out_group_fracs, $frac;
        }
    }
    # ----  {FRAG,JUMP,LONG_JUMP}_COVERAGE  ----
    
    elsif (($type eq "frag"      && $args{FRAG_COVERAGE}      ne "") ||
           ($type eq "jump"      && $args{JUMP_COVERAGE}      ne "") ||
           ($type eq "long_jump" && $args{LONG_JUMP_COVERAGE} ne ""))
    {
        my $coverage = (($type eq "frag") ? $args{FRAG_COVERAGE} :
                        ($type eq "jump") ? $args{JUMP_COVERAGE} : $args{LONG_JUMP_COVERAGE});
        my $n_bases = 0;
        for (my $i = 0; $i != $n_sub_groups; $i++) {
            my $rec = $sub_db->record_get($i);
            $n_bases += $rec->{read_count} * $rec->{mean_read_len};
        }
        abort(Tag() . "No bases found for type '$type'.") 
            unless ($n_bases > 0);

        my $frac = $coverage * $args{GENOME_SIZE} / $n_bases;
        @out_group_fracs = map $frac, (1..$n_sub_groups);
    }
    # ----  FRACTIONS  ----

    elsif ($args{FRACTIONS}) 
    {
        @out_group_fracs = map filter_frac($fracs->{$_}), @{$sub_db->{group_name}};
    }
    # ----  {FRAG,JUMP,LONG_JUMP}_FRAC  ----

    elsif (($type eq "frag"      && $args{FRAG_FRAC}      ne "") ||
           ($type eq "jump"      && $args{JUMP_FRAC}      ne "") ||
           ($type eq "long_jump" && $args{LONG_JUMP_FRAC} ne ""))
    {
        my $frac = filter_frac(($type eq "frag") ? $args{FRAG_FRAC} :
                               ($type eq "jump") ? $args{JUMP_FRAC} : $args{LONG_JUMP_FRAC});
        @out_group_fracs = map $frac, (1..$n_sub_groups);
    }
    # ----  ALL READS  ----

    else 
    {
        @out_group_fracs = map 1, (1..$n_sub_groups);
    }
    

    # ---- Validate fractions
    for (my $i = 0; $i != $n_sub_groups; $i++) {
        if ($out_group_fracs[$i] > 1) {
            my $group_name = $out_group_names[$i];
            print_warning("Fraction of reads to be merged for group '$group_name' is larger than 100%.");
            $out_group_fracs[$i] = 1;
        }
    }
    
    # ---- Build the arg strings
    my $groups_arg = globable_str_from_array_or_value(@out_group_names);
    my $fracs_arg  = globable_str_from_array_or_value(@out_group_fracs);


    # ---- Add report entries
    my %libraries = ();
    my $reads_in_all = 0;
    my $reads_out_all = 0;
    my $bases_in_all = 0;
    my $bases_physical_all = 0;
    my $min_read_len_all = 0;
    my $max_read_len_all = 0;
    my $Q1_all = 0;
    my $Q2_all = 0;
    my $n_overlap_all = 0;
    
    my $G = $args{GENOME_SIZE};
    for (my $i = 0; $i != $n_sub_groups; $i++) {
        my $group     = $out_group_names[$i];
        my $rec_group = $groups_db->sub_db_select("group_name", $group)->record_get(0);
        my $library   = $rec_group->{library_name};
        my $rec_lib   = $libs_db->sub_db_select("library_name", $library)->record_get(0);
        
        my $rec = $db_report->record_new();
        $rec->{group}            = $group;
        $rec->{library}          = $library;
        $rec->{type}             = $rec_lib->{type};
        $rec->{pairing}          = ($rec_lib->{paired} ? 
                                    $rec_lib->{read_orientation} 
                                    : "unpaired");
        
        $rec->{reads_in}         = $rec_group->{read_count};
        $rec->{reads_out}        = int($rec_group->{read_count} * $out_group_fracs[$i]);
        $rec->{fraction}         = sprintf("%5.1f %%", 100 * $out_group_fracs[$i]);

        $rec->{min_mean_max_len} = sprintf("[%4d %4d %4d]", 
                                           $rec_group->{min_read_len}, 
                                           $rec_group->{mean_read_len},
                                           $rec_group->{max_read_len});
        
        $rec->{coverage}         = ($G && $G > 0) ? 
            sprintf("%6.1fx", $rec->{reads_out} * $rec_group->{mean_read_len} / $G) : "-";

        my $bases_physical = (($rec_lib->{paired} ? 0.5 : 1.0) 
                              * $rec->{reads_out} * ($rec_lib->{insert_size} ? 
                                                     $rec_lib->{insert_size} :
                                                     $rec_lib->{frag_size} ?
                                                     $rec_lib->{frag_size} : 
                                                     $rec_group->{mean_read_len}));

        $rec->{phys_cov}         = ($G && $G > 0) ? sprintf("%6.1fx", $bases_physical / $G) : "-";
        
        $rec->{library_size}     = 
            ($rec_lib->{insert_size} ? 
             sprintf("%5d +- %4d", $rec_lib->{insert_size}, $rec_lib->{insert_stddev}) :
             ($rec_lib->{frag_size} ? 
              sprintf("%5d +- %4d", $rec_lib->{frag_size}, $rec_lib->{frag_stddev}) :
              ""));

        $rec->{mean_Qs}          = ($rec_group->{mean_r1_qual} || $rec_group->{mean_r2_qual} ?
                                    sprintf("(%2.0f %2.0f)", 
                                            $rec_group->{mean_r1_qual},
                                            $rec_group->{mean_r2_qual}) :
                                    "no quals");
        
        $rec->{prefix_out}       = $type . "_reads_orig";

        my $frac_overlap = ($type eq "frag" ? 
                            gaussian_prob_le($rec_lib->{frag_size}, 
                                             $rec_lib->{frag_stddev},
                                             2 * $rec_group->{mean_read_len})
                            : 0);

        $rec->{overlap} = ($type eq "frag" ? 
                           sprintf("%4.0f %%", 100 * $frac_overlap) : "");
        
        $db_report->record_add($rec);


        # ---- for global report

        $libraries{$library} = 1;
        $reads_in_all += $rec->{reads_in};
        $reads_out_all += $rec->{reads_out};
        $bases_in_all += $rec->{reads_in} * $rec_group->{mean_read_len};
        $bases_physical_all += $bases_physical;
        $min_read_len_all = $rec_group->{min_read_len}  
            if ($rec_group->{min_read_len} < $min_read_len_all || $min_read_len_all == 0);
        $max_read_len_all = $rec_group->{max_read_len}  
            if ($rec_group->{max_read_len} > $max_read_len_all);
        $Q1_all += $rec_group->{mean_r1_qual} ? $rec_group->{mean_r1_qual} * $rec->{reads_in} : 0;
        $Q2_all += $rec_group->{mean_r2_qual} ? $rec_group->{mean_r2_qual} * $rec->{reads_in} : 0;
        $n_overlap_all += $frac_overlap * $rec->{reads_in};
    }


    # ---- Global report
    
    if ($reads_in_all) {
        my $mean_read_len_all = $bases_in_all / $reads_in_all;
        my $mean_r1_qual_all = $Q1_all / $reads_in_all;
        my $mean_r2_qual_all = $Q2_all / $reads_in_all;
        
        
        my $rec_all = $db_report_all->record_new();
        $rec_all->{n_groups}    = $n_sub_groups;
        $rec_all->{n_libraries} = scalar keys %libraries;
        $rec_all->{prefix_out}  = $type . "_reads_orig";
        $rec_all->{reads_in}    = $reads_in_all;
        $rec_all->{reads_out}   = $reads_out_all;
        $rec_all->{fraction}    = sprintf("%5.1f %%", 100.0 * $reads_out_all / $reads_in_all);
        $rec_all->{coverage}    = ($G && $G > 0) ? 
            sprintf("%6.1fx", $reads_out_all * $mean_read_len_all / $G) : "-";
        $rec_all->{phys_cov}    = ($G && $G > 0) ?  sprintf("%6.1fx", $bases_physical_all / $G) : "-";
        
        $rec_all->{overlap}     = ($type eq "frag" ? sprintf("%4.0f %%", 100.0 * $n_overlap_all / $reads_in_all) : ""); 
        $rec_all->{min_mean_max_len} = sprintf("[%4d %4d %4d]", 
                                               $min_read_len_all,
                                               $mean_read_len_all,
                                               $max_read_len_all);
        
        $rec_all->{mean_Qs}     = ($mean_r1_qual_all || $mean_r2_qual_all ?
                                   sprintf("(%2.0f %2.0f)", $mean_r1_qual_all, $mean_r2_qual_all) :
                                   "no quals");     
        $db_report_all->record_add($rec_all);
    }
    
    return ($groups_arg, $fracs_arg);
}





sub arg_group_names
{
    my ($args) = @_;
    return [($args{GROUPS}) ?
            array_from_ref_or_value($args{GROUPS}) :
            @{new_naif_db()
                  ->from_csv_or_die($args{IN_GROUPS_CSV}, $in_groups_db_field_validator)
                  ->{group_name}}];
}






sub group_lib_db
{
    my ($group_names, $groups_db, $libs_db, $fields) = @_;
    
    my $groups_sub_db = $groups_db
        ->sub_db_select("group_name", @$group_names);

    my $db = new_naif_db($fields);
 
    my $n = $groups_sub_db->size();
    for (my $i = 0; $i < $n; $i++) 
    {
        my $group_rec =  $groups_sub_db->record_get($i);
        my $lib_name = $group_rec->{library_name};
        my $lib_rec = $libs_db->sub_db_select("library_name", $lib_name)->record_get(0);

        my %rec = ();
        map { $rec{$_} = (exists $group_rec->{$_}) ? $group_rec->{$_} : $lib_rec->{$_} } @$fields;
        
        $db->record_add(\%rec);
    }
    return $db;
}






sub arg_group_fractions_coverages
{
    my ($args, $group_names) = @_;

    # if GROUPS was not specified don't look at GROUPS_{FRAC,COVERAGE}
    return ("", "")  if (!$args{GROUPS}); 

    my $n_groups = scalar @$group_names;

    my $fracs = "";
    my $covs = "";

    if ($args{FRACTIONS}) {
        $fracs = {};
        my @arg_frac = array_from_ref_or_value($args{FRACTIONS});
        if (@arg_frac == 1) {
            map $fracs->{$_} = $arg_frac[0], @$group_names;
        }
        elsif (@arg_frac == $n_groups) {
            map $fracs->{$group_names->[$_]} = $arg_frac[$_], (0..$n_groups-1);
        }
        else {
            abort(Tag() . "Number of GROUPS != number of FRACTIONS.");
        }
    }

    if ($args{COVERAGES}) {
        $covs = {};
        my @arg_cov = array_from_ref_or_value($args{COVERAGES});
        if (@arg_cov == 1) {
            map $covs->{$_} = $arg_cov[0], @$group_names;
        }
        elsif (@arg_cov == $n_groups) {
            map $covs->{$group_names->[$_]} = $arg_cov[$_], (0..$n_groups-1);
        }
        else {
            abort(Tag() . "Number of GROUPS != number of COVERAGES.");
        }
    }

    return ($fracs, $covs);
}



# deal with fractions like 45% just as with 0.45
sub filter_frac
{
    my ($frac) = @_;
    return 0.01 * $1 if ($frac =~ /^(\S+)%$/);
    return $frac;
}





# see wikipedia for erf(x) approximation

sub gaussian_prob_le
{
    my ($mu, $sig, $x) = @_;
    return ($x > $mu ? 1.0 : $x == $mu ? 0.5 : 0.0) unless ($sig != 0);

    my $PI = 3.14159265358979323844;
    my $a = 8 * ($PI - 3) / (3 * $PI * (4 - $PI));

    my $z = ($x - $mu) / ($sig * sqrt(2.0));
    my $z2 = $z * $z;
    my $sgn = ($z == 0 ? 0 : ($z > 0 ? +1 : -1));
    my $erf = $sgn * sqrt(1.0 - exp(-$z2 * (4.0 / $PI + $a * $z2) / (1.0 + $a * $z2)));

    return 0.5 + 0.5 * $erf;

}




