#!/usr/bin/perl
use Getopt::Long;
use strict;
use POSIX;
use File::Basename;

sub getTime{
	my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
	my @weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
	my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
	my $year = 1900 + $yearOffset;
	my $theTime = "$hour:$minute:$second, $weekDays[$dayOfWeek] $months[$month] $dayOfMonth, $year";
	return $theTime; 

}

sub execCommand {
	my $retval;
	my ($command) = @_;
	print LOGFILE getTime() ."\t\t". $command . "\n";
	print "Exec at " . getTime() . " $command...\n";
	$retval  = `$command`;
	if ( $? == -1 ) {
		print "command failed: $!";
	}
	chomp $retval;
	print $retval . "\n";
	return $retval;
}

sub usage{
	die("Usage: \nMandatory options:\n\t--mode: can be all|count|split|cluster|center|reconstruct|validate|misc\n\t--k: kmer size\n\t--tau\nOptional options:\n\t--trimend : bases to chop off\n\t--qv_off: ascii offset for quality vals (usually 33 or 64 (default))\n\t--dir: alternate directory where to look for files\n\t--cutoff\n\t--mem : in gigs\n\t--trimstart : bases to chop off from the start (def 0)\n\t--ref : reference genome for validation (default is ecoli)\n\t--threads : number of threads (currently not used)\n\t--valk : k used in validation step (def is same as k)\n");
}

my $tau;
my $readLen;
my $freqCutoff = -1;
my $startTrimVal = 0;
my $endTrimVal = -1;
my $kmerval = -1;
my $valk;
my $refgen = "/home/pashadag/ec/data/ecoli";
my %options=();
my $i;
my $command;
my $mode ;
my $altDir = ".";
my $readsMaster = "files.txt";
my $qvoffset = 64;
my $mem = 0;
my $threads = 1;

my $message = "";
foreach (@ARGV) { $message .= "$_ " }

GetOptions ('mode=s' => \$mode, 'dir=s' => \$altDir, 'tau=i' => \$tau, 'threads=i' => \$threads, 'cutoff=f' => \$freqCutoff, 'trimstart=i' => \$startTrimVal, 'trimend=i' => \$endTrimVal, 'qv_off=i' => \$qvoffset, 'k=i' => \$kmerval, 'mem=f' => \$mem, 'ref=s' => \$refgen, 'valk=i' => \$valk ) or usage();


if (!defined($mode) or !defined($kmerval)) {
	usage();
}


open ( LOGFILE, ">>log.txt") or die ("$0 : failed to open log file for output: $!\n");
print LOGFILE getTime() . "\t hammer called with $message\n";

if ($mode eq "misc") {
	goto MISC;
}

my $meminfull;
if ($mem == 0) {
	$meminfull = 4000000000;
} else {
	$meminfull = $mem . "000000000";
}

my $base="reads";
my $preprocReads = $base . ".prepro";
my $kmerfile = "kmers";

if ($mode eq "all") {
	goto PREPRO;
} elsif ($mode eq "count") {
		goto COUNT;
} elsif ($mode eq "split") {
		goto SPLIT;
} elsif ($mode eq "split3") {
	goto SPLIT2;
} elsif ($mode eq "cluster") {
	goto CLUSTER;
} elsif ($mode eq "center") {
	goto CENTER;
} elsif ($mode eq "reconstruct") {
	goto RECONSTRUCT;
} elsif ($mode eq "validate") {
	goto VALIDATE;
} elsif ($mode eq "misc") {
	goto MISC;
} else {
	print STDERR "Unknown mode: $mode.\n";
	die();
}

################################
####   PREPRO           ########
####   Files used: original read files
####   Files created: reads.prepro, 
################################
PREPRO:

if (!(-e $readsMaster)) {
	my ($base,$path) = fileparse($readsMaster);  
	$readsMaster = $altDir . "/" . $base;
}


open ( READMASTERFILE, "<$readsMaster") or die("$0 : failed to open file $readsMaster: $!\n");
my $readFiles;
while (my $line = <READMASTERFILE>) {
	chomp $line;
	$readFiles .= "$line ";
}
close(READMASTERFILE);



$command = "cat $readFiles | ./preproc $startTrimVal $endTrimVal $qvoffset | grep -v 'N' > $preprocReads";
execCommand($command);


################################
####   COUNT            ########
####   Files used: reads.prepro #
####   Files created: kmers    #
################################
COUNT:
if (!(-e $preprocReads)) {
	my ($base,$path) = fileparse($preprocReads);  
	$preprocReads = $altDir . "/" . $base;
}

print $preprocReads . "\n";
my @lsarray = split(/ /, `ls -l $preprocReads`);
my $filesize = $lsarray[4];
print $filesize . "\n";
chomp $filesize;
$command = "cat $preprocReads | ./getmers $kmerval $filesize $meminfull | grep -v 'N' > $kmerfile";
execCommand($command);



################################
####   SPLIT            ########
####   Files used: kmers
####   Files created: reads.?  
################################
SPLIT:

if (!(-e $kmerfile)) {
	my ($base,$path) = fileparse($kmerfile);  
	$kmerfile = $altDir . "/" . $base;
}

my $sortOpt = "";
if ($mem != 0) {
	my $mem2 = floor( (1000 *$mem) / ($tau + 1) );
	$sortOpt = "-S $mem2" . "M";
}


open ( CMDFILE, ">cmd.txt") or die ("$0 : failed to open cmd.txt for output: $!\n");
for ($i = 0; $i < $tau + 1 ; $i++) {
	$command = "cat $kmerfile | ./splinter $i " . ($tau + 1) . " | LC_ALL=C  sort $sortOpt -k1,1 -T . > $base.$i";
	print CMDFILE "$command & \n";
}
print CMDFILE "wait\n";
close(CMDFILE);

$command = "cat cmd.txt | sh";
execCommand($command);

SPLIT2:

my $numbytes = `ls -l reads.? | awk '{print \$5}' | awk '{ s += \$1 } END{print s}'`;
chomp $numbytes;
$numbytes = int($numbytes / (4* $threads )) + 1000;
$command = "cat reads.? | split -a 2 -C $numbytes -d  - readsdup\.";
execCommand($command); 



################################
####   CLUSTER                 #
####  needed: reads.prepro     #
####  needed: reads.?          #
####  needed: kmers
####  output: reads.uf         #
################################

CLUSTER:
if (!(-e $kmerfile)) {
	my ($base,$path) = fileparse($kmerfile);  
	$kmerfile= $altDir . "/" . $base;
}

if (!(-e $preprocReads)) {
	my ($base,$path) = fileparse($preprocReads);  
	$preprocReads = $altDir . "/" . $base;
}

#my $readBase = $base;
my $readBase = "readsdup";
if (!(-e "$base.0")) {
	$readBase = "$altDir/readsdup";
} 

die if (!($tau >= 1));
my $threads2 = `ls -1 $readBase.* | wc -l`;
chomp $threads2;
$command = "./cluster $tau $kmerval $kmerfile $readBase $meminfull 0 $threads2 ";
execCommand($command);




################################
####   CENTER           ########
#    needed: uf file     #
#    output uf.corr file #
################################
CENTER:
my $uffile = "reads.uf";
if (!(-e $uffile)) {
	$uffile = $altDir . "/reads.uf";
}

$command = "./center $uffile $freqCutoff $threads";
execCommand($command);





################################
####   RECONSTRUCT      ########
#  req: uf.corr file  #
# out: fixed file#
################################
RECONSTRUCT:
my $sortOpt = "";
if ($mem != 0) {
	$sortOpt = "-S $mem" . "G";
}

my $ufcorfile = "reads.uf.corr";
if (!(-e "$ufcorfile.0")) {
	$ufcorfile = $altDir . "/reads.uf.corr";
}

$command = " cat $ufcorfile.? |  grep -v '^\$' | awk '{print \$5}' | grep -v rem | uniq  >reads.fixed.kmers";
execCommand($command);

#$command = "cat reads.uf.corr.* |  grep -v '^\$' | awk '{print \$3, \$4}' |  grep -v rem  | multiply 5 > reads.fixed.full5"; execCommand($command);


exit;

################################
####   VALIDATE         ########
#  req: hammer files        #
#  out: stats to the screen #
################################
VALIDATE:
if (!defined($valk)) {
	$valk = $kmerval;
}
my $ufcorfile = "reads.uf.corr";
if (!(-e "$ufcorfile.0")) {
	$ufcorfile = $altDir . "/reads.uf.corr";
}


$command = " cat $base.fixed.kmers| ./checkss  -k $valk -g $refgen > scores$valk & ";
execCommand($command);

#$command = "cat $base.fixed.kmers | ./checkss  -k 35 -g $refgen > scores35 & "; execCommand($command);

#$command = "cat $ufcorfile.* | grep ingleton | dbt map -g $refgen -i raw -l 2 | awk '{print $8, $10}' | sort | uniq -c | awk '{print $1}' | tr '\\n' ' ' > breakdowns_old & "; execCommand($command);

$command = "cat $ufcorfile.* | grep -v '^\$'  | grep -v change | dbt map -g ~/ec/data/ecoli -l 2 | dbt map -g ~/ec/data/ecoli -l 3 > mapped ";
execCommand($command);

$command = "cat  mapped | awk '{print \$8, \$10, \$11}' | sort | uniq -c > breakdown;  ";
execCommand($command);

$command = "cat  mapped | grep ingle | awk '{print \$4, \$10}' | awk '{ if (\$1 >= 10) \$1 = 10; print \$1, \$2 }' | sort -k2,2 -k1n,1 | uniq -c > distribmult &";
execCommand($command);

$command = "cat  mapped | grep ingle | awk '{print \$5, \$10}' | awk '{ if (\$1 >= 2) \$1 = 2; print int(\$1*10), \$2 }' | sort -k2,2 -k1n,1 | uniq -c > distribfreq &";
execCommand($command);

exit;

MISC:
$command = "cat reads.uf.corr.* |  grep -v '^\$' | awk '{print \$3, \$4}' |  grep -v rem  | multiply 5 > reads.fixed.full5";
execCommand($command);


exit;

