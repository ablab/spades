#!/usr/bin/perl

#my @children;

use Cwd 'abs_path';
use File::Basename;
my $script_dir = dirname(abs_path($0));

system("gcc $script_dir/pstree.c -o $script_dir/pstree");

my $child = fork;
exec "@ARGV 1>>rc_stdout.txt 2>>rc_stderr.txt" if $child == 0;
print "started...\n";

my $stat_proc = fork;
if ($stat_proc == 0) {
  exec "python $script_dir/stats_counter.py $child 1>>rc_stats.txt 2>>rc_stats.txt";
}

waitpid($child, 0);
waitpid($stat_proc, 0);
print "finished. Results are in rc_stats.txt\n";
