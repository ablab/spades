#!/usr/bin/perl -w
###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2011) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################
#
# CacheGroups.pl
#
# Tool to manipulate the cache groups database. 
#
# 2010-07   Filipe Ribeiro
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




# CONTROL BEGINS HERE
# Parse command-line options of the form KEY=value.
# This function comes from the ArachneArgs.pm module.

my %args = getCommandArguments
    (CACHE_DIR             => { value => "",
                                help  => "Directory where reads in AllPaths format are cached."},
     ACTION                => { value => "List",
                                help  => "Action to perform: List, Add, Remove." },
#    add args:
     PICARD_TOOLS_DIR      => { value => "",
                                help  => "The Picard tools directory." }, 
     TMP_DIR               => { value => "",
                                help  => "A large temporary directory for Java modules." },
     IN_GROUPS_CSV         => { value => "in_groups.csv",
                                help  => "Comma-separated-value info on read groups to import." }, 
     INCLUDE_NON_PF_READS  => { value => "True",
                                help  => "Whether or not to include non-PF reads." },
     PHRED_64              => { value => "False",
                                help  => "True: fastq quals are encoded with PHRED 64. False: PHRED 33." },
     FORCE_PHRED           => { value => "False",
                                help  => "True: accepts specified PHRED encoding. False: Aborts is incorrect PHRED detected." },
     REVERT_SAM            => { value => "",
                                help  => "Whether to revert alignment information in SAM files." },
     OVERWRITE             => { value => "False",
                                help  => "Whether or not to overwrite cache entries." },
     DB_ONLY               => { value => "False", 
                                help  => "Update the databases only, without reconverting BAMs." },
     SAVE_COMPRESSED_FASTQ => { value => "False",
                                help  => "Whether or not to save compressed copies of fastq files." },
     SAVE_INTERMEDIATES    => { value => "False",
                                help  => "Whether or not to keep intermediate files." },
     HOSTS                 => { value => "", 
                                help  => "e.g. \"2,4.host1,2.host2\"." },
#    list + remove args:
     PROJECT               => { value => "",
                                help  => "Optional project name." },
     ORGANISM              => { value => "",
                                help  => "Optional organism name." },
     LIBRARIES             => { value => "",
                                help  => "Optional list of libraries names." },
     GROUPS                => { value => "",
                                help  => "Optional list of groups names." },                    
     TYPE                  => { value => "",
                                help  => "Optional type ('frag' or 'jump') of reads." },
     READ_LEN              => { value => 0,
                                help  => "Optional read length of reads." },
     DISPLAY_INPUTS        => { value => "False",
                                help  => "Display input .csv files that correspond to cache entries." },
#    other 
     JAVA_MEM_GB           => { value => 8,
                                help  => "Memory to reserve for java (only for SAM/BAM conversion)." }, 
     DRY_RUN               => { value => "False",
                                help  => "If 'True' only say what would be done." },
     VERBOSE               => { value => "True" });


set_verbose($args{VERBOSE});



# ---- Make sure CACHE_DIR is defined. 

my $CACHE_DIR = get_directory_from_arg_or_env(\%args, "CACHE_DIR");
directory_is_absolute_or_die($CACHE_DIR);

my $PICARD_TOOLS_DIR = $args{PICARD_TOOLS_DIR};



# ---- Do stuff

if ($args{ACTION} eq "List") 
{
    groups_list(); 
}
elsif ($args{ACTION} eq "Add") 
{
    groups_add();
}
elsif ($args{ACTION} eq "Remove") 
{
    groups_remove();
}
else 
{
    abort(Tag() . "ACTION '$args{ACTION}' not in {'List', 'Add', 'Remove'}."); 
}






sub groups_list
{

    # ---- Validate TYPE.

    my %types = ( frag => 1, jump => 1, "" => 1 );
    
    abort(Tag() . "Invalid TYPE. Should be one of ('" . ( join "', '", keys %types) . "').")
        unless exists $types{$args{TYPE}};




    # ---- Open the cache databases from their .csv files.
    
    my $libs_db_fn = "$CACHE_DIR/libraries.csv";
    my $libs_db = new_naif_db()->from_csv_or_continue($libs_db_fn, $libs_db_field_validator);
    my $n_db_libs = $libs_db->size();
    print_action("$n_db_libs cached libraries loaded.");
    
    
    my $groups_db_fn = "$CACHE_DIR/groups.csv";
    my $groups_db = new_naif_db()->from_csv_or_continue($groups_db_fn, $groups_db_field_validator);
    my $n_db_groups = $groups_db->size();
    print_action("$n_db_groups cached groups loaded.");
    
    


    my $libs_sub_db;
    my $groups_sub_db;
    
    if ($args{GROUPS}) 
    {
        # ---- Fetch group records for all the groups matching the various criteria.
        
        $groups_sub_db = $groups_db
            ->sub_db_select("group_name", $args{GROUPS});
        
        
        # ---- Fetch only the library records that match the various criteria.
        
        $libs_sub_db = $libs_db
            ->sub_db_select("library_name", $groups_sub_db->{library_name});
    }
    else 
    {
        # ---- Fetch only the library records that match the various criteria.
        
        my $libs_tmp_db = $libs_db
            ->sub_db_select("library_name", $args{LIBRARIES})
            ->sub_db_select("project_name", $args{PROJECT})
            ->sub_db_select("organism_name", $args{ORGANISM});
        
        
        
        # ---- Select libraries with the specified type.
        
        if ($args{TYPE} eq "frag") {
            $libs_sub_db = $libs_tmp_db->sub_db_select_test("insert_size", sub { return ! $_[0]; });
        }
        elsif ($args{TYPE} eq "jump") {
            $libs_sub_db = $libs_tmp_db->sub_db_select_test("insert_size", sub { return $_[0] && $_[0] > 0; });
        }
        else {
            $libs_sub_db = $libs_tmp_db;
        }
        
        abort(Tag() . "No library matches the specified parameters.")
            unless $libs_sub_db->size();
        
        # ---- Fetch group records for all the groups matching the various criteria.
        
        $groups_sub_db = $groups_db
            ->sub_db_select("library_name", $libs_sub_db->{library_name});
    }
    





    # ---- Print all the groups
    if ($args{DISPLAY_INPUTS})
    {
        my @in_libs_fields = ("library_name",
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
                              "genomic_end");

        my %all_libs = ();
        my $n_groups = $groups_sub_db->size();
        for (my $i = 0; $i < $n_groups; $i++) 
        {
            my $group_rec = $groups_sub_db->record_get($i);
            my $lib_name = $group_rec->{library_name};
            $all_libs{$lib_name} = 1;
        }


        my $in_libs_db = new_naif_db(\@in_libs_fields);

        foreach my $lib_name (keys %all_libs)
        {
            my $lib_rec = $libs_sub_db->sub_db_select("library_name", $lib_name)->record_get(0);
            my %rec = ();
            map { $rec{$_} = $lib_rec->{$_} } @in_libs_fields;
            $in_libs_db->record_add(\%rec);
        }
        print "\n";
        $in_libs_db->to_csv("-", \@in_libs_fields);



        my @in_groups_fields = ("group_name",
                                "library_name",
                                "file_name");

        my $in_groups_db = new_naif_db(\@in_groups_fields);

        for (my $i = 0; $i < $n_groups; $i++) 
        {
            my $group_rec = $groups_sub_db->record_get($i);
            my %rec = ();
            map { $rec{$_} = $group_rec->{$_} } ("group_name", "library_name");
            $rec{file_name} = "./$rec{group_name}.unmapped.bam";
            $in_groups_db->record_add(\%rec);
        }
        print "\n";
        $in_groups_db->to_csv("-", \@in_groups_fields);
    }
    else     
    {
        my @fields = ("import_date",
                      "proj",
                      "group", 
                      "gid", 
                      "library", 
                      "lid", 
                      "organism",
                      "type",
                      "pairing",
                      "insert_size", 
                      "frag_size", 
                      "n_reads",
                      "min_mean_max_len",
		      "mean_Qs");

        my $db = new_naif_db(\@fields);

        my $n = $groups_sub_db->size();
        for (my $i = 0; $i < $n; $i++) 
        {
            my $rec_group = $groups_sub_db->record_get($i);
            my $lib_name = $rec_group->{library_name};
            my $rec_lib = $libs_sub_db->sub_db_select("library_name", $lib_name)->record_get(0);

            if (!$args{READ_LEN} ||  
                ($args{READ_LEN} >= $rec_group->{min_read_len} &&
                 $args{READ_LEN} <= $rec_group->{max_read_len})) 
            {
                
                my $rec = $db->record_new();
                $rec->{import_date}      = $rec_group->{import_date};
                $rec->{proj}             = $rec_lib->{project_name};
                $rec->{group}            = $rec_group->{group_name};
                $rec->{gid}              = $rec_group->{group_id};
                $rec->{library}          = $rec_lib->{library_name};
                $rec->{lid}              = $rec_lib->{library_id};
                $rec->{organism}         = $rec_lib->{organism_name};
                $rec->{type}             = $rec_lib->{type};
                $rec->{pairing}          = $rec_lib->{read_orientation};
                $rec->{insert_size}      = ($rec_lib->{insert_size}) ?
                    sprintf("%5d +- %4d", $rec_lib->{insert_size}, $rec_lib->{insert_stddev}) :
                    "";
                $rec->{frag_size}        = ($rec_lib->{frag_size}) ?
                    sprintf("%4d +- %3d", $rec_lib->{frag_size}, $rec_lib->{frag_stddev}) :
                    "";
                $rec->{n_reads}          = $rec_group->{read_count};
                $rec->{min_mean_max_len} = sprintf("[%4d %4d %4d]", 
                                                   $rec_group->{min_read_len}, 
                                                   $rec_group->{mean_read_len},
                                                   $rec_group->{max_read_len});
                
                $rec->{mean_Qs}          = ($rec_group->{mean_r1_qual} || $rec_group->{mean_r2_qual} ?
                                            sprintf("(%2.0f %2.0f)", 
                                                    $rec_group->{mean_r1_qual},
                                                    $rec_group->{mean_r2_qual}) :
                                            "no quals");
                
                
                $db->record_add($rec);
            }
        }
        $db->to_csv("-", \@fields, "#");
    }

}













sub groups_remove
{
    # ---- Validate GROUPS.
    abort(Tag() . "Must specify GROUPS.")
        unless ($args{GROUPS} ne "");


    # ---- Load the group database.
    my $groups_db_fn = "$CACHE_DIR/groups.csv";
    my $groups_db = new_naif_db()->from_csv_or_continue($groups_db_fn, $groups_db_field_validator);
    my $n_db_groups = $groups_db->size();
    print_action("$n_db_groups cached groups loaded.");
    print "\n";

    # ---- Remove groups.
    
    my $n_removed = 0;
   
    my @names = array_from_ref_or_value($args{GROUPS});
    foreach my $name (@names) 
    {
        my $removed = cache_db_rec_erase($groups_db_fn, "group", $name);

        if ($removed) {
            print_action("Removed group '$name'.");
            print `rm -v $CACHE_DIR/$name.*`;
            $n_removed++;
        }
        else {
            print_error("Didn't find group '$name' to remove.");
        }
    }    
    
    print_action("$n_removed groups removed.");
    return $n_removed;
}











sub groups_add
{
    # ---- Check HOSTS argument.

    my $hosts = parse_hosts_arg($args{HOSTS});


    # ---- Open the input group database from .csv file.

    my $in_groups_db_fn = $args{IN_GROUPS_CSV};
    my $in_groups_db = new_naif_db()->from_csv_or_die($in_groups_db_fn, $in_groups_db_field_validator);

    abort(Tag() . "There are repeated group names in your '$in_groups_db_fn'.")
        unless $in_groups_db->is_key("group_name");

    my $n_in_groups = $in_groups_db->size();
    print_action("$n_in_groups input groups loaded.");
    print_action("$n_in_groups input groups loaded from '$in_groups_db_fn':");
    $in_groups_db->to_csv("-", $in_groups_db_fields, "*"); 








    # ---- Open the cache databases from their .csv files.

    my $groups_db_fn = "$CACHE_DIR/groups.csv";
    my $groups_db = new_naif_db()->from_csv_or_continue($groups_db_fn, $groups_db_field_validator);
    my $n_db_groups = $groups_db->size();
    print_action("$n_db_groups cached groups loaded.");
    

    my $libs_db_fn = "$CACHE_DIR/libraries.csv";
    my $libs_db = new_naif_db()->from_csv_or_continue($libs_db_fn, $libs_db_field_validator);
    my $n_db_libs = $libs_db->size();
    print_action("$n_db_libs cached libraries loaded.");



    # ---- Make sure all libraries in in_groups_db are defined in libs_db.

    if (1) {
        my $in_groups_lib_names_freq = $in_groups_db->values_frequencies("library_name");
        my $lib_names_freq = $libs_db->values_frequencies("library_name");
        my $n_errors = 0;
        foreach my $lib_name (keys %$in_groups_lib_names_freq)
        {
            if (! exists $lib_names_freq->{$lib_name}) {
                print_error(Tag() . "Library name '$lib_name' in '$in_groups_db_fn' undefined in '$libs_db_fn'.");
                $n_errors++;
            }
        }
        abort(Tag() . "$n_errors undefined libraries.")
            if ($n_errors);
    }


    # ---- Make sure if picard tools are available in case we're importing bam files

    if (1) {
        my @fns = @{$in_groups_db->{file_name}};
        my $do_bams = 0;
        map { $do_bams = 1 if (/.bam$/i); } (@fns);
        picard_tools_verify() if $do_bams;
    }



    # ---- Parse REVERT_SAM argument.

    my $revert_sam = parse_revert_sam_arg($args{REVERT_SAM}, $in_groups_db->size());


    # ---- database for reporting purposes

    my @fields_report = ("group", 
                         "library", 
                         "organism",
                         "project",
                         "type", 
                         "insert_size",
                         "frag_size",
                         "pairing",
                         "genomic_range",                         
                         "status", 
                         "revert_sam", 
                         "action", 
                         "notes");
    my $db_report = new_naif_db(\@fields_report);

    # ---- Loop over each data file in need of importing to the cache
    #print_action(Tag() . "Checking which groups to add to cache database.");

    my @recs_in_group = ();

    my $n_warnings = 0;
    my $separator = "-" x 80 . "\n";




    for (my $i = 0; $i != $n_in_groups; $i++)
    {    
        my $rec_in_group = $in_groups_db->record_get($i);

        my $file_name    = $rec_in_group->{file_name}; # the original data file
        my $group_name   = $rec_in_group->{group_name};     
        my $library_name = $rec_in_group->{library_name};
        my $rec_lib = $libs_db->sub_db_select("library_name", $library_name)->record_get(0); 

        my $rec_report = $db_report->record_new();

        $rec_report->{group}         = $group_name;
        $rec_report->{library}       = $library_name;
        $rec_report->{organism}      = $rec_lib->{organism_name};
        $rec_report->{project}       = $rec_lib->{project_name};
        $rec_report->{type}          = $rec_lib->{type};
        $rec_report->{insert_size}   = ($rec_lib->{insert_size}) ?
            sprintf("%5s +- %4s", $rec_lib->{insert_size}, $rec_lib->{insert_stddev}) :
            "";
        $rec_report->{frag_size}     = ($rec_lib->{frag_size}) ?
            sprintf("%4s +- %3s", $rec_lib->{frag_size}, $rec_lib->{frag_stddev}) :
            "";
        $rec_report->{pairing}       = ($rec_lib->{read_orientation} eq "outward" ? "out (rev)" : "in (!rev)");
        $rec_report->{genomic_range} = sprintf("[ %3s : %3s ]", 
                                               $rec_lib->{genomic_start},
                                               $rec_lib->{genomic_end});
        $rec_report->{revert_sam}    = ($revert_sam->[$i] ? "yes" : "no");

        #print $separator;
        #group_print($rec_in_group); 
        #print "\n";
        #library_print($rec_lib);
        #print "\n";
        

        # check if file exists
        my @file_names = glob $file_name;
        if (!-e $file_name && !@file_names)
        {
            print_error("Error: group '$group_name': Can't find '$file_name'.");
            $n_warnings++;
            $rec_report->{status} = "new";
            $rec_report->{action} = "skip";
            $rec_report->{notes} .= "(can't find file)";
        }
        else  
        {   
            # ---- Check if file already has an entry in groups_db

            if ((my $rec = $groups_db->sub_db_select("group_name", $group_name))->size())
            {
                if ($args{OVERWRITE}) {
                    print_action("Converting: group '$group_name' is already in the database and it WILL be REPLACED.");
                    push @recs_in_group, $rec_in_group;
                    $rec_report->{status} = "OLD";
                    $rec_report->{action} = "OVERWRITE";
                }
                else {
                    print_action("Skipping: group '$group_name' is already in the database.");
                    $rec_report->{status} = "OLD";
                    $rec_report->{action} = "SKIP";
                }

                if ($rec->{file_name}[0] ne $file_name) {
                    print_error("Warning: group '$group_name' exists in groups_db but with a DIFFERENT file name.");
                    $rec_report->{notes} .= "(!= file_name)";
                    $n_warnings++;
                }
                if ($rec->{library_name}[0] ne $library_name) {
                    print_error("Warning: group '$group_name' exists in groups_db but with a DIFFERENT library name.");
                    $rec_report->{notes} .= "(!= library_name)";
                    $n_warnings++;
                }

            }
            else {
                print_action("Converting: group '$group_name' will be added to the database.");
                push @recs_in_group, $rec_in_group;
                $rec_report->{status} = "NEW";
                $rec_report->{action} = "ADD";
            }
        }

        $db_report->record_add($rec_report);

    }
    print $separator;


    # ---- Import data files to cache in parallel by forking.

    my $n_errors = fork_ImportToCache($PICARD_TOOLS_DIR, $hosts, $revert_sam, 
                                      \@recs_in_group, $libs_db, $groups_db_fn, $db_report);



    # ---- Final outputs.

    $groups_db = new_naif_db()->from_csv_or_continue($groups_db_fn, $groups_db_field_validator);    
    print_action(Tag() . "Final size of groups database: ", $groups_db->size());


    
    $db_report->to_csv("-", \@fields_report, "=");

    print_error(Tag() . "There were $n_warnings warnings (see above).")
        if ($n_warnings);

    abort(Tag() . "$n_errors errors occurred while converting groups. Aborting.")
        if ($n_errors);
    
}

















# SUBROUTINES











# Run ConvertToFastbQualb by forking
#
# Paradigm is:
#
#  while (not done) {
#    if (my $pid = fork()) {        # parent
#      ...                          # record child info
#    }
#    else {                         # child
#      exec("ConvertToFastbQualb ...");  
#    }
#
#    # check, without waiting, if any child is zombie and reap it.
#    if ((my $pid = waitpid(-1, WNOHANG)) > 0) {
#      ...                          # check results
#    }


sub fork_ImportToCache
{
    my ($PICARD_TOOLS_DIR, $hosts, $revert_sam, 
        $recs_in_group, $libs_db, $groups_db_fn, $db_report) = @_;
    
    my $n_hosts = @$hosts;           # maximum number of hosts available
    my @hosts_used = map 0, @$hosts;
    my %children = ();
    
    my $n_todo = @$recs_in_group;
    my $n_started = 0;
    my $n_running = 0;

    my $n_errors = 0;

    my $LOG_DIR = "$CACHE_DIR/log";
    mkdir $LOG_DIR unless (-d $LOG_DIR);


    while ($n_started < $n_todo || $n_running) 
    {

        # ---- check if there is stuff that can be processed

        if ($n_started < $n_todo && $n_running < $n_hosts) 
        {



            # ---- Find first available host.
            my $i_host = 0; 
            $i_host++ while ($i_host != $n_hosts && $hosts_used[$i_host]);
            
            abort(Tag() . "Weird... Can't find a host to fork.  I should try to spoon.")
                unless ($i_host < $n_hosts);
            




            # ---- Fork

            if (my $pid = fork())     # Parent process
            {
                $hosts_used[$i_host] = 1;      # set usage flag
                $children{$pid} = {};
                $children{$pid}{i_host} = $i_host;
                $children{$pid}{i_rec}  = $n_started;
                $n_started++;
                $n_running++;
            }
            else                      # Forked child
            {
                my $host = $hosts->[$i_host];
                my $rec_in_group = $recs_in_group->[$n_started];
                my $group_name   = $rec_in_group->{group_name};
                my $library_name = $rec_in_group->{library_name};

                # Get library info for this read set.
                my $rec_lib       = $libs_db->sub_db_select("library_name", $library_name)->record_get(0);

                my $reverse_reads = ($rec_lib->{read_orientation} eq "outward" ? "True" : "False");
                my $trim_start    = ($rec_lib->{genomic_start} ? $rec_lib->{genomic_start} : 0);
                my $trim_end      = ($rec_lib->{genomic_end}   ? $rec_lib->{genomic_end}   : 0);
                my $paired        = $rec_lib->{paired};

                my $fastb_fn = "$CACHE_DIR/$group_name.fastb";
                my $qualb_fn = "$CACHE_DIR/$group_name.qualb";
                my $done_fn  = "$CACHE_DIR/$group_name.done";


                # Cleanup previous done files
                map { unlink $_ if (-e $_) } glob "$done_fn*";



                if ($args{DB_ONLY}) # Update groups_db assuming conversion has been done
                {
                    print_action(Tag() . "Skipping import for group '$group_name'.");
                    system "touch", $done_fn
                        if (-e $fastb_fn and -e $qualb_fn);
                }
                else # Convert the data file
                {
                    print_action(Tag() . "Importing group '$group_name'.");
                    my $data_fn = $rec_in_group->{file_name};
                    $data_fn = $ENV{PWD} . "/" . $data_fn unless ($data_fn =~ /^\//);
                
                    # picard tools temporary directory
                    my $TMP_DIR = $args{TMP_DIR};
                    if ($TMP_DIR eq "") {
                        $TMP_DIR = "$CACHE_DIR/picard_tmp";
                        mkdir $TMP_DIR
                            unless (-e $TMP_DIR);
                    }

                    # Call ConvertToFastbQualb. This is the time-consuming part.
                    my $cmd = ("$FindBin::Bin/ConvertToFastbQualb.pl" . 
                               " PICARD_TOOLS_DIR=$PICARD_TOOLS_DIR" .
                               " TMP_DIR=$TMP_DIR" .
                               " DATA_FILE=\"$data_fn\"" . 
                               " PAIRED=$paired" .
                               " OUT_HEAD=$CACHE_DIR/$group_name" .
                               " REVERSE_READS=$reverse_reads" .
                               " TRIM_START=$trim_start" .
                               " TRIM_END=$trim_end" .
                               " INCLUDE_NON_PF_READS=$args{INCLUDE_NON_PF_READS}" .
                               " PHRED_64=$args{PHRED_64}" .
                               " FORCE_PHRED=$args{FORCE_PHRED}" .
                               " REVERT_SAM=$revert_sam->[$n_started]" .
                               " OVERWRITE=$args{OVERWRITE}" .
                               " SAVE_COMPRESSED_FASTQ=$args{SAVE_COMPRESSED_FASTQ}" .
                               " SAVE_INTERMEDIATES=$args{SAVE_INTERMEDIATES}" .
                               " JAVA_MEM_GB=$args{JAVA_MEM_GB}" .
                               " DRY_RUN=$args{DRY_RUN}" .
                               " VERBOSE=$args{VERBOSE}" .
                               " | tee $LOG_DIR/$group_name.log");
                    run_or_die($host eq "" ? $cmd : "ssh $host -x '$cmd'");

                    abort(Tag() . "ConvertToFastbQualb.pl call failed to generate '$CACHE_DIR/$group_name.fastb'.")
                        unless (-e "$CACHE_DIR/$group_name.fastb");


                }

    
                exit 0 if ($args{DRY_RUN});

                
                if (-e $done_fn) { # $done_fn created by ConvertToFastbQualb.pl

                    #print_subaction("Computing stats on group '$group_name'.");
                    my $stats = compute_group_stats($CACHE_DIR, $group_name, $reverse_reads, $host);
                    if ($stats->{read_count} > 0) {

                        map $rec_in_group->{$_} = $stats->{$_}, keys %$stats;
                        
                        #print_action(Tag() . "Adding record for group '$group_name' to groups database.");
                        

                        # $outcome -> { new =>, equal =>, updated => };
                        my $outcome = cache_db_rec_add($groups_db_fn, 
                                                       $rec_in_group,
                                                       "group",
                                                       $args{OVERWRITE});
                        
                        # this hack signals the outcome to the parent process
                        map { system "touch", "$done_fn.$_" if $outcome->{$_} } keys %$outcome;
                    }
                    else {
                        print_error(Tag() . "'$fastb_fn' has no reads.");
                    }

                    exit 0;
                }
                else {
                    print_error(Tag() . "ConvertToFastbQualb.pl failed for group '$group_name'.");
                    exit 1;
                }
            }


        }
        else # All hosts are busy or there is nothing left to start... sit back and relax.
        {
            sleep(1);    
        }
            





        
        # ---- The parent tests for dead children and reaps them. 

        if ((my $pid = waitpid(-1, WNOHANG)) > 0)
        {
            $n_running--;
            $hosts_used[$children{$pid}{i_host}] = 0;  # clear up usage flag

            my $rec_in_group = $recs_in_group->[$children{$pid}{i_rec}];
            my $group_name = $rec_in_group->{group_name};

            my $i_rec = $db_report->record_index("group", $group_name);
            my $rec_report = $db_report->record_get($i_rec);
            
            if ($args{DRY_RUN}) {
                $rec_report->{notes} = "[DRY_RUN][db not updated]";
            }
            else {
                my $done_fn = "$CACHE_DIR/$group_name.done";
                if (-e $done_fn) # child succeded!
                {  
                    $rec_report->{notes} .= (-e "$done_fn.updated") ? "[db updated]" : "[db not updated]";
                    map { unlink $_ } glob "$done_fn*";
                }
                else # child failed!
                {  
                    $rec_report->{notes} .= "[** FAILED **][db not updated]";
                    $n_errors++;
                }
            }
            

            $db_report->record_set($i_rec, $rec_report);

        }            
    }  #  while ($n_started < $n_todo || $n_running) 

    return $n_errors;
}









# parse_hosts_args
#
# Expects a string of the form "2,4.host1,2.host2", which translates to
# 
#  3 processes on the localhost
#  4 processes on wga7
#  8 processes on crd8 
#
sub parse_hosts_arg
{
    my ($s)= @_;

    return [""] unless $s;
    
    my @hosts = ();

    my @s = split "\,", $s; 

    foreach (@s) {
       
        if (/^(\d+)$/) {
            push @hosts, map "", (1..$1);
        }
        elsif (/^(\d+)\.(\w+)$/) {
            push @hosts, map $2, (1..$1);
        }
        elsif (/^(\w+)$/) {
            push @hosts, $1;
        }
        else {
            abort(Tag() . "Invalid HOSTS '$s' parameter. Example: '2,4.host1,2.host2'");
        }
    }
    return \@hosts;
}






sub parse_revert_sam_arg
{
    my ($revert_sam_arg, $n_in_groups) = @_;

    my @revert_sam = array_from_ref_or_value($revert_sam_arg);

    if (scalar @revert_sam == 1 && $revert_sam[0] eq "") 
    {
        @revert_sam = map 0, (1..$n_in_groups); # no reversal will take place for any group
    }
    else {

        my $good = 1;
        foreach my $rs (@revert_sam) 
        {
            $good = 0 unless ($rs =~ /True|true|False|false|0|1/);
        }
        abort(Tag() . "REVERT_SAM must be a set (e.g. '{\"False\", \"True\", \"False\", ...}') with the same number of elements as groups.") 
            unless (@revert_sam == $n_in_groups and $good);
    }
    return \@revert_sam;
}


sub picard_tools_verify
{
    # ---- Make sure PICARD_TOOLS_DIR is defined and the needed Picard modules are present.
    
    if (!$args{DB_ONLY}) 
    {
        if ($PICARD_TOOLS_DIR eq "") 
        {
            abort(Tag() . "You didn't specify PICARD_TOOLS_DIR nor set the environment variable \$ALLPATHS_PICARD_TOOLS_DIR.")
                unless defined $ENV{ALLPATHS_PICARD_TOOLS_DIR};
            $PICARD_TOOLS_DIR = $ENV{ALLPATHS_PICARD_TOOLS_DIR};
        }
        
        
        my $n_errors = 0;
        my @modules = ("SamToFastq.jar");
        push @modules, "RevertSam.jar" if ($args{REVERT_SAM} ne "");
        
        foreach my $module (@modules) 
        {
            if (! -e "$PICARD_TOOLS_DIR/$module") 
            {
                print_error(Tag() . "Can't find Picard module '$module' in '$PICARD_TOOLS_DIR'.");
                $n_errors++;
            }
        }
        abort(Tag() . "Found $n_errors errors.")
            if ($n_errors);
    }

}







# Call QualbStats on a .qualb and return:
# mean quality score of 1st and 2nd reads 

sub compute_group_stats
{
    my ($cache_dir, $group_name, $rc, $host) = @_;

    my $nr = 0;
    my $nb = 0;
    my $len_min = 0;
    my $len_max = 0;
    my $len_mean = 0;
    my $q1 = 0;
    my $q2 = 0;
    my $q_max = 0;

    if (-e "$cache_dir/$group_name.qualb") {
        my $dn = "$cache_dir/qualb_stats/";
        mkdir $dn unless (-d $dn);
        
        # Run QualbStats.  This can be a little slow.
        my $cmd = ("$FindBin::Bin/QualbStats" .
                   " HEAD=$cache_dir/$group_name" .
                   " OUT_HEAD=$dn/$group_name" .
                   " RC=$rc" .
                   " TEE=$dn/$group_name.out");
        my @QualbStats_output = run_or_die($host eq "" ? $cmd : "ssh $host -x '$cmd'");
        
        # Parse the output of QualbStats and get the values we want.
        foreach (@QualbStats_output) {
            if (/^\s+(\S+)\s+reads\s*$/) {
                $nr = $1;
            }
            elsif (/^\s+(\S+)\s+total length\s*$/) {
                $nb = $1;
            }
            elsif (/^\s*(\S+)\s+minimum length\s*$/) {
                $len_min = $1;
            }
            elsif (/^\s*(\S+)\s+maximum length\s*$/) {
                $len_max = $1;
            }
            elsif (/^\s*(\S+)\s+mean length\s*$/) {
                $len_mean = $1;
            }
            elsif (/^\s*(\S+)\s+max Q\s*$/) {
                $q_max = $1;
            }
            elsif (/^\s*(\S+)\s+(\S+)\s+mean Q\s*$/) {
                ($q1, $q2) = ($1, $2);
            }
        }
        
        if ($q_max > 45) {
            print_warning(Tag() . "Maximum quality score = $q_max!");
            print_subaction("Are you sure the quality scores in the fastq are encoded PHRED " 
                            . (($args{PHRED_64}) ? "64" : "33") . "?\n");
            print_subaction("Please verify and try again (use option PHRED_64 as appropriate).\n");
        }  

    }
    else {
        # Run FastbStats.  This can be a little slow.
        my $cmd = ("$FindBin::Bin/FastbStats" .
                   " FASTB=$cache_dir/$group_name.fastb");
        my @FastbStats_output = run_or_die($host eq "" ? $cmd : "ssh $host -x '$cmd'");
        
        # Parse the output of QualbStats and get the values we want.
        foreach (@FastbStats_output) {
            if (/^\s+total objects count:\s+(\S+)\s*$/) {
                $nr = $1;
            }
            elsif (/^\s+total length:\s+(\S+)\s*$/) {
                $nb = $1;
            }
            elsif (/^\s+avg length:\s+(\S+)\s*$/) {
                $len_mean = $1;
            }
            elsif (/^\s+range:\s+\[(\d+),\s+(\d+)\]\s*$/) {
                $len_min = $1;
                $len_max = $2;
            }
        }
    }


    if ($nr == 0 || 
        $nb == 0 ||
        $len_min == 0 ||
        $len_max == 0 ||
        $len_mean == 0) 
    {
        print_error(Tag() . "Something went wrong parsing the output of QualbStats for '$group_name'.");
        exit(1);
    }
    
    return { read_count    => $nr,
             min_read_len  => $len_min,
             max_read_len  => $len_max,
             mean_read_len => $len_mean,
	     mean_r1_qual  => $q1,
	     mean_r2_qual  => $q2 };
}

