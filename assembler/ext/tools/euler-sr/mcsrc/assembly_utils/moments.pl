#!/usr/bin/perl
# $Id: moments.perl,v 1.5 2001/01/05 22:36:44 doug Exp $
# http://www.bagley.org/~doug/shootout/

use strict;

my @nums = <STDIN>;
my $sum = 0;
foreach (@nums) { $sum += $_ }
my $n = scalar(@nums);
my $mean = $sum/$n;
my $average_deviation = 0;
my $standard_deviation = 0;
my $variance = 0;
my $skew = 0;
my $kurtosis = 0;
foreach (@nums) {    
    my $deviation = $_ - $mean;
    $average_deviation += abs($deviation);
    $variance += $deviation**2;
    $skew += $deviation**3;
    $kurtosis += $deviation**4;
}
$average_deviation /= $n;
$variance /= ($n - 1);
$standard_deviation = sqrt($variance);

if ($variance) {
    $skew /= ($n * $variance * $standard_deviation);
    $kurtosis = $kurtosis/($n * $variance * $variance) - 3.0;
}

# median
@nums = sort { $a <=> $b } @nums;
my $mid = int($n/2);
my $median = ($n % 2) ? $nums[$mid] : ($nums[$mid] + $nums[$mid-1])/2;

# mode
my $mode = 'NONE';
my $currCount = 1;
my $prevBestCount = 0;
for(my $i = 0; $i <= $#nums; $i++) {
    if ($i == $#nums or $nums[$i] != $nums[$i+1]) {
	if ($currCount > $prevBestCount) {
	    $prevBestCount = $currCount;
	    $mode = $nums[$i];
	}
	$currCount = 0;
    }
    $currCount++;
}

# n50 = weighted median
my $curr_sum = $sum;
my $j;
for ($j = $#nums; $j > 0 and $curr_sum > ($sum/2); $j--){
    $curr_sum -= $nums[$j];
}
my $n50 = $nums[$j];

# output
printf("n:                  %d\n", $n);
printf("max:                %d\n", $nums[$#nums]);
printf("mode:               %d\n", $mode);
printf("median:             %f\n", $median);
printf("mean:               %f\n", $mean);
printf("min:                %d\n", $nums[0]);
printf("n50:                %d\n", $n50);
printf("average_deviation:  %f\n", $average_deviation);
printf("standard_deviation: %f\n", $standard_deviation);
printf("variance:           %f\n", $variance);
printf("skew:               %f\n", $skew);
printf("kurtosis:           %f\n", $kurtosis);
