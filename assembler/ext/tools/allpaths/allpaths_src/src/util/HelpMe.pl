#!/usr/bin/perl -w
###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2011) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################
#
#  HelpMe.pl script that packages several logs together. 
#
#  2011-03   Filipe Ribeiro



use strict;
use FindBin;

# ---- Local libraries, from elsewhere in the BroadCRD repository.
use lib "$FindBin::Bin/../";
use lib "$FindBin::Bin/";

use PerlRunTime;    # run_or_die
use PerlUtils;      # ISO_date, abort
use ArachneArgs;    # Command-line ARG=VALUE parser


# ---- CONTROL BEGINS HERE
#      Parse command-line options of the form KEY=value.
#      This function comes from the ArachneArgs.pm module.
my %args = getCommandArguments
    (SUBDIR               => { value => undef,
                               help  => "The SUBDIR of the assembly." });

abort("Specified SUBDIR '$args{SUBDIR}' is not a directory.") 
    unless (-d $args{SUBDIR});

my $sub_dir = ($args{SUBDIR} =~ /^\//) ? $args{SUBDIR} : "$ENV{PWD}/$args{SUBDIR}";

my $log_dir = "make_log";

my $dir_info = find_base_dir($sub_dir, $log_dir) ||
    abort("Can't find '$log_dir' directory in any of the parent directories of '$sub_dir'. Aborting.");


my $full_ref_dir = $dir_info->{ref};
my $full_log_dir = "$full_ref_dir/$log_dir";


chdir $full_ref_dir;




# ---- Capturing directory listings
system "find . -type f -or -type l | xargs ls -lt > $full_log_dir/all_files.ls";



# ---- Capturing some environment settings
{
    system("free -m > $full_log_dir/free-m.out");
    system("df -h $full_ref_dir > $full_log_dir/df.out");
    system("dmesg > $full_log_dir/dmesg.out");
    system("cat /proc/sys/vm/overcommit_memory > $full_log_dir/overcommit_memory.out")
        if (-e "/proc/sys/vm/overcommit_memory");
    open(FILE, ">$full_log_dir/env.out");
    map { print FILE "$_\t$ENV{$_}\n"; } (keys %ENV);
    close FILE;
}


# ---- Collecting some files

{
    my $full_files_dir = "$full_log_dir/files";
    mkdir $full_files_dir unless (-d $full_files_dir);

    my @kspec_fns = `find . -name '*.kspec' | grep -v $log_dir`; # exclude files already in $log_dir
    chomp @kspec_fns;
    foreach my $fn (@kspec_fns) {
	my $dir = ($fn =~ /^(.+)\/[^\/]+$/)[0];
	system "mkdir -p $full_files_dir/$dir ; cp $fn $full_files_dir/$fn"; 
    }

    my @libs_ascii_dir = `find . -name '*libs.ascii' | grep -v $log_dir`; # exclude files already in $log_dir
    chomp @libs_ascii_dir;
    foreach (@libs_ascii_dir) {
	system "mkdir -p $full_files_dir/$_ ; cp $_/* $full_files_dir/$_";
    }

    my $full_cache_dir = "$full_files_dir/cache";
    mkdir $full_cache_dir unless (-d $full_cache_dir);

    foreach ("*.csv", "*.PAPI.out", "*.PAPI.err") {
	system "find . -name '$_' | grep -v $log_dir | xargs --replace=FN cp FN $full_cache_dir"; 
    } 

    my @qualb_stats_dir = `find . -name qualb_stats`;
    chomp @qualb_stats_dir;
    system "cp $qualb_stats_dir[0]/* $full_cache_dir" if (@qualb_stats_dir == 1);
    
}





# ---- Tar-ing everything together
chdir $full_ref_dir;

my $tar_fn = "help_me.tar.gz";

run_or_die("find . -name '*.report' | xargs tar cvzf $tar_fn ./$log_dir/");

print "\n";
print "Please attach the file '$full_ref_dir/$tar_fn' to your email to crdhelp\@broadinstitute.org.\n";
print "\n";






sub find_base_dir
{
    my ($dir, $base_dir) = @_;

    $dir = "$ENV{PWD}/$dir" unless $dir =~ /^\//;
    
    my %dirs = ( names => [],
                 ref   => undef);

    do {
        $dir =~ s/([^\/]+)\/?$//;
        push @{$dirs{names}}, $1;
    } until ($dir eq "/" or -d "$dir/$base_dir");

    $dirs{ref} = $dir;

    return ($dir eq "/") ? 0 : \%dirs;
}


