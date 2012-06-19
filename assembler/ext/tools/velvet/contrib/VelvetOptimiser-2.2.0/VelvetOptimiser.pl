#!/usr/bin/perl
#
#       VelvetOptimiser.pl
#
#       Copyright 2008, 2009, 2010 Simon Gladman <simon.gladman@csiro.au>
#
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.

#		Version 2.2.0

#
#   pragmas
#
use strict;
use warnings;
#
#   includes
#
use POSIX qw(strftime);
use FindBin;
use lib "$FindBin::Bin";
use threads;
use threads::shared;
use VelvetOpt::Assembly;
use VelvetOpt::hwrap;
use VelvetOpt::gwrap;
use VelvetOpt::Utils;
use Data::Dumper;
use Storable qw (freeze thaw);
use Getopt::Long;


#
#   global var decs
#

#Change the following integer when compiling Velvet with the MAXKMERLENGTH
#greater than 31 to the value you used.
my $maxhash;
my @hashvals;
my %assemblies : shared;
my %assembliesObjs;
my @Options;
my $readfile;
my $interested = 0;
my $verbose : shared;
my $hashs;
my $hashe;
my $hashstep;
my $amos;
my $vgoptions;
my $genomesize;
my @shortInserts;
my $logfile = "logfile.txt";
my $ass_num = 1;
my $categories;
my $prefix;
my $OUT;
my $logSem : shared;
our $num_threads;
my $current_threads : shared = 0;
my $opt_func;
my $opt_func2;
my $OptVersion = "2.2.0";
my $threadfailed : shared = 0;
my $finaldir;

#
#
#	main script
#
#
print STDERR "
****************************************************

           VelvetOptimiser.pl Version $OptVersion

            Simon Gladman - CSIRO 2009

****************************************************\n";

my $currfreemem = VelvetOpt::Utils::free_mem;

print STDERR "Number of CPUs available: " . VelvetOpt::Utils::num_cpu . "\n";
printf STDERR "Current free RAM: %.3fGB\n", $currfreemem;

#get the velveth and velvetg version numbers...
my $response = VelvetOpt::hwrap::_runVelveth(" ");
$response =~ /Version\s+(\d+\.\d+\.\d+)/s;
my $vhversion = $1;
unless ($vhversion){ die "Unable to find velveth, please ensure that the velvet executables are in your PATH.\n";}
$response =~ /CATEGORIES = (\d+)/;
$categories = $1;
unless($categories){ $categories = 2; }

$response =~ /MAXKMERLENGTH = (\d+)/;
$maxhash = $1;
unless($maxhash){ $maxhash = 31; }

#check the number of threads that velvet was compiled with (OMP_NUM_THREADS) if it is the OMP version
#then warn the user that -t will multiply that value and VO will use more CPUs than they think..
my $thread_per_job = $ENV{OMP_NUM_THREADS} || 1;

print STDERR "Velvet OMP compiler setting: $thread_per_job\n";

#get the options!
&setOptions();

if($prefix eq "auto"){
	$logfile = strftime("%d-%m-%Y-%H-%M-%S", localtime) . "_Logfile.txt";
} else {
	$logfile = $prefix . "_logfile.txt";
}

print "Logfile name: $logfile\n";

#open the logfile
open $OUT, ">$logfile" or die "Couldn't open $logfile for writing.\n$!\n";

#
#
#   Perform common tasks - write details to log file and screen, run velveth and vanilla velvetg
#
#

print STDERR "\nMemory use estimation only!  Script will terminate after showing results.\n\n" if($genomesize);

print STDERR "Velvet details:\n";
print STDERR "\tVelvet version: $vhversion\n";
print STDERR "\tCompiled categories: $categories\n" if $categories;
print STDERR "\tCompiled max kmer length: $maxhash\n" if $maxhash;
print STDERR "\tMaximum number of velvetinstances to run: $num_threads\n";

#let user know about parameters to run with.
print STDERR "Will run velvet optimiser with the following paramters:\n";
print STDERR "\tVelveth parameter string:\n\t\t$readfile\n";
print STDERR "\tVelveth start hash values:\t$hashs\n";
print STDERR "\tVelveth end hash value:\t\t$hashe\n";
print STDERR "\tVelveth hash step value:\t$hashstep\n\n";
if($vgoptions){
	print $OUT "\tUser specified velvetg options: $vgoptions\n";
}
if($amos){
    print STDERR "\tRead tracking for final assembly on.\n";
} else {
    print STDERR "\tRead tracking for final assembly off.\n";
}

#build the hashval array - steps too...
for(my $i = $hashs; $i <= $hashe; $i += $hashstep){
    push @hashvals, $i;
}
#check for $hashe in array..
my $max = $hashvals[$#hashvals];
if($max < $hashe){
	push @hashvals, $hashe;
}

if($genomesize){
	my $x = &estMemUse();
	printf STDERR "\nMemory use estimated to be: %.1fGB for $num_threads threads.\n\n", $x;
	if ($x < $currfreemem){
		print STDERR "You should have enough memory to complete this job. (Though this estimate is no guarantee..)\n";
		exit;
	}
	else {
		print STDERR "You probably won't have enough memory to run this job.\nTry decreasing the maximum number of threads used.\n(use the -t option to set max threads.)\n";
		exit;
	}
}


print $OUT strftime("%b %e %H:%M:%S", localtime), "\n";

#send run parameters to log file.
print $OUT "Will run velvet optimiser with the following paramters:\n";
print $OUT "\tVelveth parameter string:\n\t\t$readfile\n";
print $OUT "\tVelveth start hash values:\t$hashs\n";
print $OUT "\tVelveth end hash value:\t\t$hashe\n";
print $OUT "\tVelveth hash step value:\t\t$hashstep\n\n";
if($vgoptions){
	print $OUT "\tUser specified velvetg options: $vgoptions\n";
}
if($amos){
    print $OUT "\tRead tracking for final assembly on.\n";
} else {
    print $OUT "\tRead tracking for final assembly off.\n";
}

print STDERR strftime("%b %e %H:%M:%S", localtime), " Beginning velveth runs.\n";
print $OUT strftime("%b %e %H:%M:%S", localtime), "\n\n\tBeginning velveth runs.\n";

#now run velveth for all the hashvalues in a certain number of threads..
my @threads;
foreach my $hashval (@hashvals){
	while($current_threads >= $num_threads){
		sleep(2);
	}
	if($threadfailed){
		for my $thr (threads->list) {
			#print STDERR "Waiting for thread ",$thr->tid," to complete.\n";
			$thr->join;
		}
		die "Velveth failed to run! Must be a problem with file types, check by running velveth manually or by using -v option and reading the log file.\n";
	}	
	$threads[$ass_num] = threads->create(\&runVelveth, $readfile, $hashval, $vhversion, \$logSem, $ass_num);
	$ass_num ++;
	sleep(2);
}

for my $thr (threads->list) {
    #print STDERR "Waiting for thread ",$thr->tid," to complete.\n";
    $thr->join;
}

#now run velvetg for the all the hashvalues in a certain number of threads..
#first get velvetg's version number.

$response = VelvetOpt::gwrap::_runVelvetg(" ");
$response =~ /Version\s+(\d+\.\d+\.\d+)/s;
my $vgversion = $1;

print STDERR strftime("%b %e %H:%M:%S", localtime), " Finished velveth runs.\n";

print STDERR strftime("%b %e %H:%M:%S", localtime), " Beginning vanilla velvetg runs.\n";
print $OUT strftime("%b %e %H:%M:%S", localtime), "\n\n\tBeginning vanilla velvetg runs.\n";

foreach my $key (sort { $a <=> $b } keys %assemblies){
	while($current_threads >= $num_threads){
		sleep(2);
	}
	$threads[$ass_num] = threads->create(\&runVelvetg, $vgversion, \$logSem, $key);
	sleep(2);
}

for my $thr (threads->list) {
    #print STDERR "Waiting for thread ",$thr->tid," to complete.\n";
    $thr->join;
}


#now to thaw it all out..

foreach my $key(sort keys %assemblies){
	my $obj = bless thaw($assemblies{$key}), "VelvetOpt::Assembly";
	$assembliesObjs{$key} = $obj;
}


#find the best assembly...

#
#
#   Now perform a velvetg optimisation based upon the file types sent to velveth
#
#

#
#   get the best assembly so far...
#

my $bestId;
my $maxScore = -100;
my $asmscorenotneg = 1;

foreach my $key (keys %assembliesObjs){
	if(($assembliesObjs{$key}->{assmscore} != -1) && $asmscorenotneg){
    	if($assembliesObjs{$key}->{assmscore} > $maxScore){
        	$bestId = $key;
        	$maxScore = $assembliesObjs{$key}->{assmscore};
    	}
	}
	elsif($assembliesObjs{$key}->{n50} && $asmscorenotneg){
		if($assembliesObjs{$key}->{n50} > $maxScore){
			$bestId = $key;
			$maxScore = $assembliesObjs{$key}->{n50};
		}
	}
	else {
		$asmscorenotneg = 0;
		if($assembliesObjs{$key}->{totalbp} > $maxScore){
        	$bestId = $key;
        	$maxScore = $assembliesObjs{$key}->{totalbp};
    	}
	}
}
print "\n\nThe best assembly so far is:\n" if $interested;
print $assembliesObjs{$bestId}->toStringNoV() if $interested;

#   determine the optimisation route for the assembly based on the velveth parameter string.
my $optRoute = &getOptRoutine($readfile);

print STDERR strftime("%b %e %H:%M:%S", localtime), " Hash value of best assembly by assembly score: ". $assembliesObjs{$bestId}->{hashval} . "\n";

print $OUT strftime("%b %e %H:%M:%S", localtime), " Best assembly by assembly score - assembly id: $bestId\n";

print STDERR strftime("%b %e %H:%M:%S", localtime), " Optimisation routine chosen for best assembly: $optRoute\n";
print $OUT strftime("%b %e %H:%M:%S", localtime), " Optimisation routine chosen for best assembly: $optRoute\n";

#now send the best assembly so far to the appropriate optimisation routine...

if($optRoute eq "shortOpt"){
	
	&expCov($assembliesObjs{$bestId});
    &covCutoff($assembliesObjs{$bestId});

}
elsif($optRoute eq "shortLong"){

    &expCov($assembliesObjs{$bestId});
    &covCutoff($assembliesObjs{$bestId});

}
elsif($optRoute eq "longPaired"){
    &expCov($assembliesObjs{$bestId});
    &insLengthLong($assembliesObjs{$bestId});
    &covCutoff($assembliesObjs{$bestId});
}
elsif($optRoute eq "shortPaired"){
    &expCov($assembliesObjs{$bestId});
    &insLengthShort($assembliesObjs{$bestId});
    &covCutoff($assembliesObjs{$bestId});
}
elsif($optRoute eq "shortLongPaired"){
    &expCov($assembliesObjs{$bestId});
    &insLengthShort($assembliesObjs{$bestId});
    &insLengthLong($assembliesObjs{$bestId});
    &covCutoff($assembliesObjs{$bestId});
}
else{
    print STDERR "There was an error choosing an optimisation routine for this assembly.  Please change the velveth parameter string and try again.\n";
    print $OUT "There was an error choosing an optimisation routine for this assembly.  Please change the velveth parameter string and try again.\n";
}

#   once it comes back from the optimisation routines, we need to turn on read tracking and amos output if it was selected in the options.
#
#
#   The final assembly run!
#
#
if($amos){
    $assembliesObjs{$bestId}->{pstringg} .= " -amos_file yes -read_trkg yes";

    my $final = VelvetOpt::gwrap::objectVelvetg($assembliesObjs{$bestId});
    $assembliesObjs{$bestId}->getAssemblyDetails();
}

print STDERR strftime("%b %e %H:%M:%S", localtime), "\n\n\nFinal optimised assembly details:\n";
print $OUT strftime("%b %e %H:%M:%S", localtime), "\n\n\nFinal optimised assembly details:\n";
print STDERR $assembliesObjs{$bestId}->toStringNoV() if !$verbose;
print $OUT $assembliesObjs{$bestId}->toStringNoV() if !$verbose;
print STDERR $assembliesObjs{$bestId}->toString() if $verbose;
print $OUT $assembliesObjs{$bestId}->toString() if $verbose;
if($finaldir eq "."){
	print STDERR "\n\nAssembly output files are in the following directory:\n" . $assembliesObjs{$bestId}->{ass_dir} . "\n\n";
	print $OUT "\n\nAssembly output files are in the following directory:\n" . $assembliesObjs{$bestId}->{ass_dir} . "\n";
}
else {
	print STDERR "\n\nAssembly output files are in the following directory:\n" . $finaldir . "\n\n";
	print $OUT "\n\nAssembly output files are in the following directory:\n" . $finaldir . "\n";
}

#delete superfluous directories..
foreach my $key(keys %assemblies){
	unless($key == $bestId){ 
		my $dir = $assembliesObjs{$key}->{ass_dir};
		system('rm', '-r', '--preserve-root', $dir);
	} 
}
unless ($finaldir eq "."){
	rename $assembliesObjs{$bestId}->{ass_dir}, $finaldir;
	rename $logfile, "$finaldir/$logfile";
}

#
#
#	subroutines...
#
#
#----------------------------------------------------------------------

# Option setting routines

sub setOptions {
	use Getopt::Long;
	my $num_cpus = VelvetOpt::Utils::num_cpu;
	my $thmax = int($num_cpus/$thread_per_job);
	

	@Options = (
		{OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
		{OPT=>"v|verbose+", VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose logging, includes all velvet output in the logfile."},
		{OPT=>"s|hashs=i", VAR=>\$hashs, DEFAULT=>19, DESC=>"The starting (lower) hash value"}, 
		{OPT=>"e|hashe=i", VAR=>\$hashe, DEFAULT=>$maxhash, DESC=>"The end (higher) hash value"},
		{OPT=>"x|step=i", VAR=>\$hashstep, DEFAULT=>2, DESC=>"The step in hash search..  min 2, no odd numbers"},
		{OPT=>"f|velvethfiles=s", VAR=>\$readfile, DEFAULT=>0, DESC=>"The file section of the velveth command line."},
		{OPT=>"a|amosfile!", VAR=>\$amos, DEFAULT=>0, DESC=>"Turn on velvet's read tracking and amos file output."},
		{OPT=>"o|velvetgoptions=s", VAR=>\$vgoptions, DEFAULT=>'', DESC=>"Extra velvetg options to pass through.  eg. -long_mult_cutoff -max_coverage etc"},
		{OPT=>"t|threads=i", VAR=>\$num_threads, DEFAULT=>$thmax, DESC=>"The maximum number of simulataneous velvet instances to run."},
		{OPT=>"g|genomesize=f", VAR=>\$genomesize, DEFAULT=>0, DESC=>"The approximate size of the genome to be assembled in megabases.\n\t\t\tOnly used in memory use estimation. If not specified, memory use estimation\n\t\t\twill not occur. If memory use is estimated, the results are shown and then program exits."},
		{OPT=>"k|optFuncKmer=s", VAR=>\$opt_func, DEFAULT=>'n50', DESC=>"The optimisation function used for k-mer choice."},
		{OPT=>"c|optFuncCov=s", VAR=>\$opt_func2, DEFAULT=>'Lbp', DESC=>"The optimisation function used for cov_cutoff optimisation."},
		{OPT=>"p|prefix=s", VAR=>\$prefix, DEFAULT=>'auto', DESC=>"The prefix for the output filenames, the default is the date and time in the format DD-MM-YYYY-HH-MM_."},
		{OPT=>"d|dir_final=s", VAR=>\$finaldir, DEFAULT=>'.', DESC=>"The name of the directory to put the final output into."},
	);

	(@ARGV < 1) && (usage());

	&GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

	# Now setup default values.
	foreach (@Options) {
		if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
		${$_->{VAR}} = $_->{DEFAULT};
		}
	}
	
	print STDERR strftime("%b %e %H:%M:%S", localtime), " Starting to check input parameters.\n";
	
	unless($readfile){
		print STDERR "\tYou must supply the velveth parameter line in quotes. eg -f '-short .....'\n";
		&usage();
	}
	
    if($hashs > $maxhash){
        print STDERR "\tStart hash value too high.  New start hash value is $maxhash.\n";
        $hashs = $maxhash;
    }
    if(!&isOdd($hashs)){
        $hashs = $hashs - 1;
        print STDERR "\tStart hash value not odd.  Subtracting one. New start hash value = $hashs\n";
    }
    
    if(&isOdd($hashstep)){
		print STDERR "\tOld hash step value $hashstep\n";
		$hashstep --;
		print STDERR "\tHash search step value was odd, substracting one.  New hash step value = $hashstep\n";
	}
	if($hashstep < 2){
		$hashstep = 2;
		print STDERR "\tHash step set below minimum of 2.  New hash step value = 2\n";
	}
	if($hashstep > ($hashe - $hashs)){
		$hashstep = $hashe - $hashs;
		print STDERR "\tHash search step value was higher than start to end range.  Setting hash step to range. New hash step value = $hashstep\n";
	}
		
	if($hashe > $maxhash || $hashe < 1){
        print STDERR "\tEnd hash value not in workable range.  New end hash value is $maxhash.\n";
        $hashe = $maxhash;
    }
    if($hashe < $hashs){
        print STDERR "\tEnd hash value lower than start hash value.  New end hash value = $hashs.\n";
        $hashe = $hashs;
    }
    if(!&isOdd($hashe)){
        $hashe = $hashe - 1;
        print STDERR "\tEnd hash value not odd.  Subtracting one. New end hash value = $hashe\n";
    }
    
    if($num_threads > $thmax){
		print STDERR "\tWARNING: You have set the number of threads to use to a value greater than the number of available CPUs.\n";
		print STDERR "\tWARNING: This may be because of the velvet compile option for OMP.\n";
	}
	
	#check the velveth parameter string..
	my $vh_ok = VelvetOpt::hwrap::_checkVHString("check 21 $readfile", $categories);

	unless($vh_ok){ die "Please re-start with a corrected velveth parameter string." }
	
	print STDERR "\tVelveth parameter string OK.\n";
	
	#test if outdir exists...
	if(-d $finaldir && $finaldir ne "."){
		die "Output directory $finaldir already exists, please choose a different name and restart.\n";
	}

	print STDERR strftime("%b %e %H:%M:%S", localtime), " Finished checking input parameters.\n";
	
}

sub usage {
	print "Usage: $0 [options] -f 'velveth input line'\n";
	foreach (@Options) {
		printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
			defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
	}
	print "\nAdvanced!: Changing the optimisation function(s)\n";
	print VelvetOpt::Assembly::opt_func_toString;
	exit(1);
}
 
#----------------------------------------------------------------------


#
#	runVelveth
#

sub runVelveth{
	
	{
		lock($current_threads);
		$current_threads ++;
	}
	
	my $rf = shift;
	my $hv = shift;
	my $vv = shift;
	my $semRef = shift;
	my $anum = shift;
	my $assembly;
	
	print STDERR strftime("%b %e %H:%M:%S", localtime), "\t\tRunning velveth with hash value: $hv.\n";

    #make the velveth command line.
    my $vhline = $prefix . "_data_$hv $hv $rf";

    #make a new VelvetAssembly and store it in the %assemblies hash...
	$assembly = VelvetOpt::Assembly->new(ass_id => $anum, pstringh => $vhline, versionh =>$vv, assmfunc => $opt_func, assmfunc2 => $opt_func2);

    #run velveth on this assembly object
    my $vhresponse = VelvetOpt::hwrap::objectVelveth($assembly, $categories);

    unless($vhresponse){ die "Velveth didn't run on hash value of $hv.\n$!\n";}
	
	unless(-r ($prefix . "_data_$hv" . "/Roadmaps")){ 
		print STDERR "Velveth failed!  Response:\n$vhresponse\n";
		{
			lock ($threadfailed);
			$threadfailed = 1;
		}
	}
    
	#run the hashdetail generation routine.
    $vhresponse = $assembly->getHashingDetails();
    
	#print the objects to the log file...
	{
		lock($$semRef);
		print $OUT $assembly->toStringNoV() if !$verbose;
		print $OUT $assembly->toString() if $verbose;
	}
	
	{
		lock(%assemblies);
		my $ass_str = freeze($assembly);
		$assemblies{$anum} = $ass_str;
	}
	
	{
		lock($current_threads);
		$current_threads --;
	}
	print STDERR strftime("%b %e %H:%M:%S", localtime), "\t\tVelveth with hash value $hv finished.\n";
}

#
#	runVelvetg
#
sub runVelvetg{

	{
		lock($current_threads);
		$current_threads ++;
	}
	
	my $vv = shift;
	my $semRef = shift;
	my $anum = shift;
	my $assembly;
	
	#get back the object!
	$assembly = bless thaw($assemblies{$anum}), "VelvetOpt::Assembly";
	
	print STDERR strftime("%b %e %H:%M:%S", localtime), "\t\tRunning vanilla velvetg on hash value: " . $assembly->{hashval} . "\n";

	#make the velvetg commandline.
    my $vgline = $prefix . "_data_" . $assembly->{hashval};
	
	$vgline .= " $vgoptions";
	$vgline .= " -clean yes";

    #save the velvetg commandline in the assembly.
    $assembly->{pstringg} = $vgline;
	
	#save the velvetg version in the assembly.
	$assembly->{versiong} = $vv;

    #run velvetg
    my $vgresponse = VelvetOpt::gwrap::objectVelvetg($assembly);

    unless($vgresponse){ die "Velvetg didn't run on the directory $vgline.\n$!\n";}

    #run the assembly details routine..
    $assembly->getAssemblyDetails();

    #print the objects to the log file...
	{
		lock($$semRef);
		print $OUT $assembly->toStringNoV() if !$verbose;
		print $OUT $assembly->toString() if $verbose;
	}
	
	{
		lock(%assemblies);
		my $ass_str = freeze($assembly);
		$assemblies{$anum} = $ass_str;
	}
	
	{
		lock($current_threads);
		$current_threads --;
	}
	print STDERR strftime("%b %e %H:%M:%S", localtime), "\t\tVelvetg on hash value: " . $assembly->{hashval} . " finished.\n";
}

#
#   isOdd
#
sub isOdd {
    my $x = shift;
    if($x % 2 == 1){
        return 1;
    }
    else {
        return 0;
    }
}


#
#   getOptRoutine
#
sub getOptRoutine {

    my $readfile = shift;

    #   Choose the optimisation path depending on the types of read files sent to velveth
    #       For short only:                 shortOpt routine
    #       For short and long:             shortLong routine
    #       For short paired:               shortPaired routine
    #       For short and long paired:      longPaired routine
    #       For short paired and long:      shortPaired routine
    #       For short paired & long paired: shortlongPaired routine

    #look at velveth string ($readfile) and look for keywords from velvet manual...
    my $long = 0;
    my $longPaired = 0;
    my $shortPaired = 0;
    my $short = 0;

    #standard cases..
    if($readfile =~ /-short.? /) { $short = 1; }
    if($readfile =~ /-long /) { $long = 1; }
    if($readfile =~ /-shortPaired /) { $shortPaired = 1; }
    if($readfile =~ /-longPaired /) { $longPaired = 1; }

    #weird cases to cover the non-use of the short keyword (since its the default.)
    if(!($readfile =~ /(-short.? )|(-long )|(-shortPaired )|(-longPaired )/)) { $short = 1; } #if nothing is specified, assume short.
    if(!($readfile =~ /-short.? /) && ($readfile =~ /(-long )|(-longPaired )/)) { $short = 1; } #if long or longPaired is specified, also assum short since very unlikely to only have long...

    if($short && !($long || $longPaired || $shortPaired)){
        return "shortOpt";
    }
    elsif($short && $long && !($longPaired || $shortPaired)){
        return "shortLong";
    }
    elsif($short && $longPaired && !$shortPaired){
        return "longPaired";
    }
    elsif($short && $shortPaired && !$longPaired){
        return "shortPaired";
    }
    elsif($short && $shortPaired && $longPaired){
        return "shortLongPaired";
    }
    elsif($shortPaired && !$short && !$long && !$longPaired){
        return "shortPaired";
    }
    else {
        return "Unknown";
    }
}

#
#   covCutoff - the coverage cutoff optimisation routine.
#
sub covCutoff{

    my $ass = shift;
    #get the assembly score and set the current cutoff score.
    my $ass_score = $ass->{assmscore};
    print "In covCutOff and assembly score is: $ass_score..\n" if $interested;
	


	sub func {
		my $ass = shift;
		my $cutoff = shift;
		my $ass_score = $ass->{assmscore};
		my $ps = $ass->{pstringg};
        if($ps =~ /cov_cutoff/){
            $ps =~ s/cov_cutoff\s+\d+(\.\d+)?/cov_cutoff $cutoff/;
        }
        else {
            $ps .= " -cov_cutoff $cutoff";
        }
        $ass->{pstringg} = $ps;

        print STDERR strftime("%b %e %H:%M:%S", localtime);
		printf STDERR "\t\tSetting cov_cutoff to %.3f.\n", $cutoff;
        print $OUT strftime("%b %e %H:%M:%S", localtime);
		printf $OUT "\t\tSetting cov_cutoff to %.3f.\n", $cutoff;

        my $worked = VelvetOpt::gwrap::objectVelvetg($ass);
        if($worked){
            $ass->getAssemblyDetails();
        }
        else {
            die "Velvet Error in covCutoff!\n";
        }
        $ass_score = $ass->{assmscore};
		print $OUT $ass->toStringNoV();
		
		return $ass_score;
		
	}
	
	print STDERR strftime("%b %e %H:%M:%S", localtime), " Beginning coverage cutoff optimisation\n";
    print $OUT strftime("%b %e %H:%M:%S", localtime), " Beginning coverage cutoff optimisation\n";

	my $dir = $ass->{ass_dir};
    $dir .= "/stats.txt";
    #print "\tLooking for exp_cov in $dir\n";
    my $expCov = VelvetOpt::Utils::estExpCov($dir, $ass->{hashval});
	
	my $a = 0;
	my $b = 0.8 * $expCov;
	my $t = 0.618;
	my $c = $a + $t * ($b - $a);
	my $d = $b + $t * ($a - $b);
	my $fc = func($ass, $c);
	my $fd = func($ass, $d);

	my $iters = 1;
	
	printf STDERR "\t\tLooking for best cutoff score between %.3f and %.3f\n", $a, $b;
	printf $OUT "\t\tLooking for best cutoff score between %.3f and %.3f\n", $a, $b;
	
	while(abs($a -$b) > 1){
		if($fc > $fd){
			printf STDERR "\t\tMax cutoff lies between %.3f & %.3f\n", $d, $b;
			my $absdiff = abs($fc - $fd);
			print STDERR "\t\tfc = $fc\tfd = $fd\tabs diff = $absdiff\n";
			printf $OUT "\t\tMax cutoff lies between %.3f & %.3f\n", $d, $b;
			$a = $d;
			$d = $c;
			$fd = $fc;
			$c = $a + $t * ($b - $a);
			$fc = func($ass, $c);
		}
		else {
			printf STDERR "\t\tMax cutoff lies between %.3f & %.3f\n", $a, $c;
			my $absdiff = abs($fc - $fd);
			print STDERR "\t\tfc = $fc\tfd = $fd\tabs diff = $absdiff\n";
			printf $OUT "\t\tMax cutoff lies between %.3f & %.3f\n", $a, $c;
			$b = $c;
			$c = $d;
			$fc = $fd;
			$d = $b + $t * ($a - $b);
			$fd = func($ass, $d);
		}
		$iters ++;
	}

	printf STDERR "\t\tOptimum value of cutoff is %.2f\n", $b;
	print STDERR "\t\tTook $iters iterations\n";
	printf $OUT "\t\tOptimum value of cutoff is %.2f\n", $b;
	print $OUT "\t\tTook $iters iterations\n";

    return 1;

}

#
#   expCov - find the expected coverage for the assembly and run velvetg with that exp_cov.
#
sub expCov {

    print STDERR strftime("%b %e %H:%M:%S", localtime), " Looking for the expected coverage\n";
    print $OUT strftime("%b %e %H:%M:%S", localtime), " Looking for the expected coverage\n";

    my $ass = shift;

    #need to get the directory of the assembly and add "stats.txt" to it and then send it to
    #the histogram methods in SlugsUtils.pm...
    my $dir = $ass->{ass_dir};
    $dir .= "/stats.txt";
    my $expCov = VelvetOpt::Utils::estExpCov($dir, $ass->{hashval});

    print STDERR strftime("%b %e %H:%M:%S", localtime), "\t\tExpected coverage set to $expCov\n";
    print $OUT strftime("%b %e %H:%M:%S", localtime), "\t\tExpected coverage set to $expCov\n";

    #re-write the pstringg with the new velvetg command..
    my $vg = $ass->{pstringg};
    if($vg =~ /exp_cov/){
        $vg =~ s/exp_cov\s+\d+/exp_cov $expCov/;
    }
    else {
        $vg .= " -exp_cov $expCov";
    }

    $ass->{pstringg} = $vg;
    
    print $OUT $ass->toStringNoV();

}

#
#   insLengthLong - get the Long insert length and use it in the assembly..
#
sub insLengthLong {
    print STDERR strftime("%b %e %H:%M:%S", localtime), " Getting the long insert length\n";
    print $OUT strftime("%b %e %H:%M:%S", localtime), " Getting the long insert length\n";
    my $ass = shift;
    my $len = "auto";
    print STDERR strftime("%b %e %H:%M:%S", localtime), " Setting assembly long insert length $len\n";
    print $OUT strftime("%b %e %H:%M:%S", localtime), " Setting assembly long insert length $len\n";

    #re-write the pstringg with the new velvetg command..
    #my $vg = $ass->{pstringg};
    #if($vg =~ /ins_length_long/){
    #    $vg =~ s/ins_length_long\s+\d+/ins_length_long $len/;
    #}
    #else {
    #    $vg .= " -ins_length_long $len";
    #}
}

#
#   insLengthShort - get the short insert length and use it in the assembly..
#
sub insLengthShort {
    print STDERR strftime("%b %e %H:%M:%S", localtime), " Setting the short insert length\n";
    print $OUT strftime("%b %e %H:%M:%S", localtime), " Setting the short insert length\n";
    my $ass = shift;
	my $len = "auto";
    print STDERR strftime("%b %e %H:%M:%S", localtime), " Setting assembly short insert length(s) to $len\n";
    print $OUT strftime("%b %e %H:%M:%S", localtime), " Setting assembly short insert length(s) to $len\n";

    #re-write the pstringg with the new velvetg command..
    #my $vg = $ass->{pstringg};
    #if($vg =~ /ins_length /){
    #    $vg =~ s/ins_length\s+\d+/ins_length $len/;
    #}
    #else {
    #    $vg .= " -ins_length $len";
    #}
    #$ass->{pstringg} = $vg;
}


#
#	estMemUse - estimates the memory usage from 
#
sub estMemUse {
	
	my $max_runs = @hashvals;
	my $totmem = 0;
	#get the read lengths and the number of reads...
	#need the short read filenames...
	my ($rs, $nr) = VelvetOpt::Utils::getReadSizeNum($readfile);
	if ($max_runs > $num_threads){
		for(my $i = 0; $i < $num_threads; $i ++){
			$totmem += VelvetOpt::Utils::estVelvetMemUse($rs, $genomesize, $nr, $hashvals[$i]);
		}
	}
	else {
		foreach my $h (@hashvals){
			$totmem += VelvetOpt::Utils::estVelvetMemUse($rs, $genomesize, $nr, $h);
		}
	}
	return $totmem;
}
