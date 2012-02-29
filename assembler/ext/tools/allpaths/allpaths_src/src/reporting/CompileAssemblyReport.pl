#!/usr/bin/perl
#///////////////////////////////////////////////////////////////////////////////
#//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
#//       This software and its documentation are copyright (2009) by the     //
#//   Broad Institute.  All rights are reserved.  This software is supplied   //
#//   without any warranty or guaranteed support whatsoever. The Broad        //
#//   Institute is not responsible for its use, misuse, or functionality.     //
#///////////////////////////////////////////////////////////////////////////////
#
# Parses different output files for key statistics.
#
# outputs report to '<SUBDIR>/assembly.report' 
#
# Dec 2009 - Filipe Ribeiro - creation; FindErrors
# Jan 2010 - Filipe Ribeiro - ScaffoldAccuracy, AssemblyCoverage, AssemblyAccuracy 
# Feb 2011 - Filipe Ribeiro - Major clean up.  Move to PERFSTAT parsing.

use strict;
use FindBin;
use POSIX ":sys_wait_h";     # for the WNOHANG signal to wait pid

# Local libraries, from elsewhere in the BroadCRD repository.
use lib "$FindBin::Bin/../";
use lib "$FindBin::Bin/";

use PerlRunTime;    # run_or_die
use PerlUtils;      # ISO_date, abort
use ArachneArgs;    # Command-line ARG=VALUE parser


my %args = getCommandArguments
    (SUBDIR => { value => undef,
                 help  => "SUBDIR of the ALLPATHS-LG pipeline."});


my $subdir = ($args{SUBDIR} =~ /^\//) ? $args{SUBDIR} : "$ENV{PWD}/$args{SUBDIR}";


die "**** invalid SUBDIR directory '$subdir'.\n" unless (-d $subdir);


my ($ass_stats, $wc_time) = assembly_stats_compile($subdir);
my $mem_cpu_stats = mem_cpu_stats_compile($subdir);



# ---- build report string

my $report = "";
$report .= assembly_stats_sprint($ass_stats);
$report .= mem_cpu_stats_sprint($mem_cpu_stats, $wc_time);


# ---- print report to stdout

print "\n" . $report . "\n";


# ---- write report to 'assembly.report'

my $report_fn = "$subdir/assembly.report";
open(FILE, ">$report_fn") or die "**** Can't open '$report_fn'.\n";
print FILE $report;
close(FILE);



# ---- build and write report for dexter


if (defined $ENV{CRDPERFSTATLOG}) {
    my $ds = dexter_sprint($ass_stats, $mem_cpu_stats, disk_usage_GB($subdir), $wc_time);
    my $fn = $ENV{CRDPERFSTATLOG};
    open(FILE, ">$fn") or die "**** Can't open '$fn'.\n";
    print FILE $ds;
    close(FILE);
}


# ---- write the environment to a report file

my $env_report_fn = "$subdir/environment.report";
open(FILE, ">$env_report_fn") or die "**** Can't open '$env_report_fn'.\n";
map { printf FILE "%20s %s\n", $_, $ENV{$_} } keys %ENV;
close(FILE);








sub assembly_stats_sprint
{
    my ($stats) = @_;
    my $s = "";
    foreach my $sts (@$stats) {
        if (@{$sts->{stats}} || @{$sts->{blocks}}) {
            $s .= report_header($sts->{module}, $sts->{target});
            $s .= "\n";
            foreach my $st (@{$sts->{blocks}}) {
                $s .= "$st->{description}:\n";
                $s .= $st->{block};
            }
            foreach my $st (@{$sts->{stats}}) {
                $s .= sprintf("%15s    %s\n", $st->{value}, $st->{description});
            }
            $s .= "\n";
        }
    }
    return $s;
}



sub mem_cpu_stats_sprint
{
    my ($stats, $wc_time) = @_;
    my $s = report_header("Memory and CPU usage");
    $s .= "\n";
    $s .= sprintf("%15d    available cpus\n", number_cpus());
    $s .= sprintf("%15.1f    GB of total available memory\n", total_memory_GB());
    $s .= sprintf("%15.1f    GB of available disk space\n", available_disk_GB());
    $s .= sprintf("%15.2f    hours of total elapsed time\n", $wc_time / 3600.0); 
    $s .= sprintf("%15.2f    hours of total per-module elapsed time\n", 
                  $stats->{etime} / 3600.0);
    $s .= sprintf("%15.2f    hours of total per-module user time\n", 
                  $stats->{utime} / 3600.0);
    $s .= sprintf("%15.2f    effective parallelization factor\n", 
                  ($stats->{etime} > 0) ? $stats->{utime} / $stats->{etime} : 0);
    $s .= sprintf("%15.2f    GB memory usage peak\n", 
                  $stats->{max_vmrss} / 1024**2);

    $s .= "\n";
    return $s;
}




sub dexter_sprint
{
    my ($a_stats, $mc_stats, $du_stat, $wc_time) = @_;
    my $s = "";
    foreach my $sts (@$a_stats) {
        foreach my $st (@{$sts->{stats}}) {
            $s .= sprintf("%s\t%s\t%s\n", $st->{short}, $st->{value}, $st->{description});
        }
    }
    $s .= sprintf("%s\t%.2f\t%s\n", 
                  "utime_h", $mc_stats->{utime} / 3600.0, "user time in hours");
    $s .= sprintf("%s\t%.2f\t%s\n", 
                  "sum_etime_h", $mc_stats->{etime} / 3600.0, "sum of elapsed times in hours");
    $s .= sprintf("%s\t%.2f\t%s\n", 
                  "etime_h", $wc_time / 3600.0, "elapsed time in hours");
    $s .= sprintf("%s\t%.2f\t%s\n", 
                  "parallel_factor", 
		  (($mc_stats->{etime} > 0) ? $mc_stats->{utime} / $mc_stats->{etime} : 0), 
		  "effective parallelization factor");
    $s .= sprintf("%s\t%.2f\t%s\n", 
                  "mem_gigabytes", $mc_stats->{max_vmrss} / 1024**2, "memory usage peak in GB");
    $s .= sprintf("%s\t%.2f\t%s\n", 
                  "du_gigabytes", $du_stat, "disk usage in GB");
    return $s;
}




















sub assembly_stats_compile
{
    my ($subdir) = @_;

    # ---- find all the '*.out.*' files 
    my @out_fns = (sort { (stat($a))[9] <=> (stat($b))[9]; } 
                   (find_lower_makeinfo_files($subdir, "out"), 
                    find_upper_makeinfo_files($subdir, "out")));

    die "**** Couldn't find any '*.out.*' files to compile assembly stats. Aborting.\n"
        unless (@out_fns);

    # ---- the full wallclock time
    my $wc_time = ((stat($out_fns[-1]))[9] - (stat($out_fns[0]))[9]);

    # ---- parse all the '*.out.*' files
    my @stats = ();

    foreach my $fn (@out_fns) {
        my ($target, $module) = $fn =~ /([^\/]+)\.out\.([A-Z]\S+)$/;
        if (defined $target and defined $module) {
            my $time = (stat($fn))[9];  # modification time
            
            my $s = { time => $time,
                      module => $module,
                      target => $target,
                      stats => [],
                      blocks => []};
            
            open FILE, "<$fn";
            while (<FILE>) {
                if (/PERFSTAT: BLOCK_START \[(.+)\]\s*$/) {
                    my $block_name = $1;
                    my $block_str = "";
                    my $line = <FILE>;
                    while ($line && $line !~ /PERFSTAT: BLOCK_STOP/) {
                        $block_str .= $line;
                        $line = <FILE>;
                    }
                    push @{$s->{blocks}}, { description => $block_name,
                                            block       => $block_str };
                }
                elsif (/PERFSTAT:\s*(\S[^\[]+\S)\s*\[(.+)\] = (\S+)\s*$/) {
                    push @{$s->{stats}}, { description => $1,
                                           short       => $2,
                                           value       => $3 };
                }
            }
            close FILE;
            push @stats, $s;
        }
    }

    return (\@stats, $wc_time);
}












sub mm_file_parse
{
    my ($fn) = @_;
    my ($mod_target, $mod_name) = $fn =~ /([^\/]+)\.mm\.(\S+)$/;
    
    my $t1 = (stat($fn))[9];   # the file's modification time
    
    my @line = split ' ', `grep ' Summary: ' $fn`;  # get memmonitor data
    
    return {time0   => $line[3] - $line[4],
            time1   => $line[3],
            etime   => $line[4], 
            detime  => $line[5], 
            utime   => $line[6], 
            dutime  => $line[7], 
            vmsize  => $line[8],
            vmrss   => $line[9],
            vmdata  => $line[10],
            vmstk   => $line[11],
            vmexe   => $line[12],
            vmlib   => $line[13],
            name    => $mod_name,
            target  => $mod_target};
}



sub global_mem_stats_compile
{
    my ($mods) = @_;
    my ($max_vmrss, $max_vmsize) = (0, 0);

    # ---- find the memory peak for the greediest module.
    if (1) {
        foreach my $mod (@$mods) {
            $max_vmrss = max ($mod->{vmrss}, $max_vmrss);
            $max_vmsize = max ($mod->{vmsize}, $max_vmsize);
        }  
    }
    else {
        # ---- find worst-case max memory usage 

        # create an array with all the start(0) and end(1) times pegged with 
        # the effects on memory usage (+ on start, - on end)

        my @ts = ((map { t => $_->{time0}, dvmrss =>   $_->{vmrss}, dvmsize =>   $_->{vmsize} }, @$mods),
                  (map { t => $_->{time1}, dvmrss => - $_->{vmrss}, dvmsize => - $_->{vmsize} }, @$mods));

        my $vmrss = 0;
        my $vmsize = 0;

        my @ts_sorted = sort { $a->{t} <=> $b->{t} } @ts;

        foreach my $t (@ts_sorted) 
        {
            $vmrss += $t->{dvmrss};
            $vmsize += $t->{dvmsize};

            $max_vmrss = max($max_vmrss, $vmrss);
            $max_vmsize = max($max_vmsize, $vmsize);
        }
    }
    return ($max_vmrss, $max_vmsize);
}




sub mem_cpu_stats_compile
{
    my ($subdir) = @_;

    # ---- find all the '*.mm.*' files
    my @mm_fns = (find_lower_makeinfo_files($subdir, "mm"), 
                  find_upper_makeinfo_files($subdir, "mm"));

    die "**** Couldn't find any '*.mm.*' files to compile memory and cpu stats. Aborting.\n"
        unless (@mm_fns);

    # ---- parse all the '*.mm.*' files
    my @mods = map mm_file_parse($_), @mm_fns;

    
    my ($max_vmrss, $max_vmsize) = global_mem_stats_compile(\@mods);
      

    # ---- compile results
    my $etime = 0;
    my $utime = 0;
    foreach (@mods) 
    { 
        $etime += $_->{etime};
        $utime += $_->{utime};
    }

    return { max_vmrss  => $max_vmrss,
             max_vmsize => $max_vmsize,
             etime      => $etime,
             utime      => $utime };
}













sub find_lower_makeinfo_files 
{
    my ($curdir, $label) = @_;
    $curdir =~ s/\/?$//;

    my @fns = (-d "$curdir/makeinfo") ? `find $curdir/makeinfo -name '*.$label.*'` : ();
    #print "lower $curdir\n";

    chomp @fns;

    foreach my $dir (`ls $curdir`) 
    {
        chomp $dir;
        if (-d "$curdir/$dir" and 
            $dir ne "seed" and 
            $dir ne "makeinfo") {
            push @fns, find_lower_makeinfo_files("$curdir/$dir", $label)
        }
    }
    return @fns;
}



sub find_upper_makeinfo_files 
{
    my ($curdir, $label) = @_;

    return () if (!$curdir or -d "$curdir/make_log");

    $curdir =~ s/\/[^\/]+\/?$//;   # get rid of top level directory

    my @fns = (-d "$curdir/makeinfo") ? `find $curdir/makeinfo -name '*.$label.*'` : ();

    #print "higher $curdir\n";
    chomp @fns;
    
    return (@fns, find_upper_makeinfo_files($curdir, $label));
}



sub disk_usage_GB
{
    my ($dir) = @_;
    
    while ($dir && !-d "$dir/make_log") {
	$dir =~ s/\/[^\/]+\/?$//;   # get rid of top level directory
    }
    return (-d "$dir/make_log") ? (split " ", `du -s $dir`)[0] / 1024.0**2 : 0;
}



sub total_memory_GB
{
    return (split " ", `head -1 /proc/meminfo`)[1] / 1024.0**2
        if (-e '/proc/meminfo');
    return 0;
}

sub available_disk_GB
{
    return (split " ", `df -k . | tail -1`)[2] / 1024.0**2;
}



sub number_cpus
{
    if (-e '/proc/cpuinfo') {
        my @a = `grep processor /proc/cpuinfo`;
        return scalar @a;
    }
    return 1;
}




sub report_header
{
    my ($module, $target) = @_;
    return "------------------ $module" . ((@_ >= 2) ? " -> $target\n" : "\n");
}



sub format_t_DHMS
{
    my ($t_s) = @_;
    my $s = "";
    
    my $m_s = 60;
    my $h_m = 60;
    my $d_h = 24;
    my $h_s = $h_m * $m_s;
    my $d_s = $d_h * $h_s;

    my $t_d = $t_s / $d_s;
    my $t_h = $t_s / $h_s;
    my $t_m = $t_s / $m_s;

    $s = sprintf("%dd %02dh", $t_d, $t_h % $d_h) if ($t_s >= $d_s);
    $s = sprintf("%2dh %02dm", $t_h, $t_m % $h_m) if ($t_s < $d_s);
    $s = sprintf("%2dm %02ds", $t_m, $t_s % $m_s) if ($t_s < $h_s);
    $s = sprintf("%ss", $t_s)                    if ($t_s < $m_s);

    return $s;
}





sub t_fmt
{
    my ($t) = @_;
    my @aux = localtime($t);
    $aux[5] += 1900; # dates start in 1900
    $aux[4]++; # month would otherwise be in [0..11]
    return sprintf "%4d.%02d.%02d-%02d:%02d:%02d", @aux[5,4,3,2,1,0];
}


sub max
{
    my ($a, $b) = @_;
    return ($a > $b) ? $a : $b;
}
