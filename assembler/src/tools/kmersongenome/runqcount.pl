#!/usr/bin/perl -w
use strict;
use warnings;
use Cwd;
use Getopt::Long;


if ($#ARGV < 3) {
	die "\nUsage: runqcount.pl min_k max_k output_prefix input_fasta";
}

my $min_k = $ARGV[0];
my $max_k = $ARGV[1];
my $pref = $ARGV[2];
my $inp = $ARGV[3];
my $mem = "50G";
my $indbytes = "7";
my $numthreads = "32";

my $jelly_command = "jellyfish";

my $result = GetOptions ("-s=s" => \$mem,
			 "-c=s" => \$indbytes,
			 "-t=s" => \$numthreads
);


sub my_exec {
	my $cmd = shift;
	print $cmd . "\n";
	#system $cmd;
}

for (my $k = $min_k; $k <= $max_k; $k += 2) {
	my_exec( "jellyfish count --both-strands -m $k -o $pref$k -c $indbytes -s $mem -t $numthreads $inp" );
}

die("All done.");

