#!/usr/bin/env perl

# ****************** galign_with_fixed_gap.pl **********************
# ****************** Author: Chun-Nuan Chen **********************
# ****************** Last Modified: 10/8/2001 ********************

use strict;

##################### let the code start #######################
my @DD;
my $i;
my $j;
my $v=1;
my $arg0=$ARGV[0];
my $StringA;
my $StringB;
my @SeqA;
my @SeqB;
my $StrLenA;
my $StrLenB;
my $gap='-';
my @linklist;
my $result_number=0;
my $DEBUG=0;

my $msg1="Global alignment with fixed gap penalty\n";
my $msg2="Version 1.0 (Chun-Nuan Chen), 10/8/2001\n";
my $msg3="Usage:\n    $0 [-debug] sequence1 sequence2\n\n";
my $usage=($msg1).($msg2).($msg3);

die $usage if @ARGV<2 || $ARGV[0] =~ /^-h/;

if($arg0 =~ /^-d/ || $arg0 =~ /^-D/ ){
	$DEBUG=1;
	$StringA=$ARGV[1];
	$StringB=$ARGV[2];
}else{
	$StringA=$ARGV[0];
	$StringB=$ARGV[1];
}	       

$StrLenA=length $StringA;
$StrLenB=length $StringB;

push(@SeqA,$gap);
push(@SeqB,$gap);

for my $k (0..$StrLenA-1){
       push(@SeqA,substr($StringA,$k,1));
}
for my $k (0..$StrLenB-1){
       push(@SeqB,substr($StringB,$k,1));
}

print "The sequences to be globally aligned are:\n";
print "\t\SeqA: @SeqA\n";
print "\t\SeqB: @SeqB\n";

my $LenA=scalar @SeqA;
my $LenB=scalar @SeqB;

if($DEBUG){
	print "\$LenA is $LenA\n";
	print "\$LenB is $LenB\n";
}

print "\nThe dynamic programming matrix is: \n";
my $dash=('-')x(($LenB+1)*4);
print "$dash\n";
printf "%4s",' ';
foreach $j (0..$LenB-1){
	printf "%4s", $SeqB[$j];
}
print "\n";
#print "\t   @SeqB\n";
foreach $i (0..$LenA-1){
	printf "%4s", $SeqA[$i];
	$DD[$i]=[]; ## assign a reference to an anonymous array
	$linklist[$i]=[];## assign a reference to an anonymous array
	foreach $j (0..$LenB-1){
		$linklist[$i]->[$j]={};## assign a reference to an anonymous hash
		if($i==0){
	   		$DD[$i]->[$j]=$j*$v;
		        $linklist[$i]->[$j]{left}=1;	
		}
		elsif($j==0){
			$DD[$i]->[$j]=$i*$v;
			$linklist[$i]->[$j]{upper}=1;
		}
		else{
			$DD[$i]->[$j]=&mindist(
				$DD[$i-1]->[$j-1]+dist($SeqA[$i],$SeqB[$j]),
			   	$DD[$i-1]->[$j]+dist($SeqA[$i],$gap),
			  	$DD[$i]->[$j-1]+dist($gap,$SeqB[$j]),
				$i,
				$j
			);
       
		}
		printf "%4s",  @{$DD[$i]}[$j];
	}
	print " \n";
}

print "$dash\n";

my $BackSeqA=[];## anonymous arrary reference
my $BackSeqB=[];## anonymous arrary reference
&traceback($BackSeqA,$BackSeqB,$LenA-1,$LenB-1);

sub traceback(){
	if($DEBUG){
		print "\n############### beginning of traceback ###############\n";
	}
	##### convert array refereces to local-scoped array varibles \
	##### so that the recursive calls will not "remember" the array\
	##### values referenced. In another word, the referece will cause\
	##### problem when transferring the Seq values from 'diag' to 'upper'\
	##### or to 'left', and vice versa (see below)
	
	my @Back_SeqA_diag=@{$_[0]};
	my @Back_SeqA_upper=@{$_[0]};
	my @Back_SeqA_left=@{$_[0]};
	my @Back_SeqB_diag=@{$_[1]};
	my @Back_SeqB_upper=@{$_[1]};
	my @Back_SeqB_left=@{$_[1]};
	my $i=$_[2];
	my $j=$_[3];

	if($i>0 || $j>0){
		if($linklist[$i]->[$j]{diag}){
			unshift(@Back_SeqA_diag,$SeqA[$i]);
			unshift(@Back_SeqB_diag,$SeqB[$j]);
			if($DEBUG){
				print "[$i,$j] going diagonal\n";
				print "inside \@Back_SeqA_diag is @Back_SeqA_diag \n";
				print "inside \@\@Back_SeqB_diag is @{\@Back_SeqB_diag} \n";
			}
			&traceback(\@Back_SeqA_diag,\@Back_SeqB_diag,$i-1,$j-1);
		}
		if($linklist[$i]->[$j]{upper}){
			unshift(@Back_SeqA_upper,$SeqA[$i]);
			unshift(@Back_SeqB_upper,$gap);
			if($DEBUG){
				print "[$i,$j] going upper\n";
				print "inside \@Back_SeqA_upper is @Back_SeqA_upper \n";
				print "inside \@Back_SeqB_upper is @Back_SeqB_upper \n";
			}
			&traceback(\@Back_SeqA_upper,\@Back_SeqB_upper,$i-1,$j);
		}
		if($linklist[$i]->[$j]{left}){
			unshift(@Back_SeqA_left,$gap);
			unshift(@Back_SeqB_left,$SeqB[$j]);
			if($DEBUG){
				print "[$i,$j] going left\n";
				print "inside \@Back_SeqA_left is @Back_SeqA_left \n";
				print "inside \@Back_SeqB_left is @Back_SeqB_left \n";
			}
			&traceback(\@Back_SeqA_left,\@Back_SeqB_left,$i,$j-1 );
		}
		
	}
	else{	
		$result_number++;
		print "The No. $result_number alignment is:\n";
		print "\t\SeqA:  @{$_[0]} \n";
                print "\t      ";
                foreach my $i (0..(scalar @{$_[0]})-1){
                    if($_[0]->[$i] eq $_[1]->[$i]){
                            print " |";
                    }else{
                            print "  ";
                    }
                }
                print "\n";
		print "\t\SeqB:  @{$_[1]} \n\n";
	}
	if($DEBUG){
		print "################## end of traceback ##################\n";
	}
}

sub dist(){
	my $arg1=$_[0];
	my $arg2=$_[1];
	if($arg1 eq $arg2){
		return 0;
	}
	else{
		return 1;
	}
}

sub mindist(){
	my $arg1=$_[0];
	my $arg2=$_[1];
	my $arg3=$_[2];
	my $i=$_[3];
	my $j=$_[4];
	my $temp_min=&minvalue($arg1,$arg2,$arg3);
	if($arg1==$temp_min){
		$linklist[$i]->[$j]{diag}=1;
		#print "\$linklist[$i]-[$j]-{\"diag\"} is $linklist[$i]->[$j]{diag}\n";
		
	}
	if($arg2==$temp_min){
		$linklist[$i]->[$j]{upper}=1;
		#print "\$linklist[$i]-[$j]-{\"upper\"} is $linklist[$i]->[$j]{upper}\n";
	}
	if($arg3==$temp_min){
		$linklist[$i]->[$j]{left}=1;
		#print "\$linklist[$i]-[$j]-{\"left\"} is $linklist[$i]->[$j]{left}\n";
	}
	return $temp_min;
}		
		
sub minvalue(){
	my $arg1=$_[0];
	my $arg2=$_[1];
	my $arg3=$_[2];
	my $temp_min=$arg1;
	if($arg2 <$temp_min){
		$temp_min=$arg2;
	}
	if($arg3 < $temp_min){
		$temp_min=$arg3;
	}
	return $temp_min;
}
print "$dash\n";
1;
################# end of code #####################

