#       VelvetOpt::hwrap.pm
#
#       Copyright 2008 Simon Gladman <simon.gladman@csiro.au>
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
#
#		Version 1.1 - 14/07/2010 - Added support for changing input file types
#		Version 1.2 - 11/08/2010 - Changed velveth help parser for new velvet help format
#									Thanks to Alexie Papanicolaou - CSIRO for the patch.

package VelvetOpt::hwrap;

=head1 NAME

VelvetOpt::hwrap.pm - Velvet hashing program wrapper module.

=head1 AUTHOR

Simon Gladman, CSIRO, 2007, 2008.

=head1 LICENSE

Copyright 2008 Simon Gladman <simon.gladman@csiro.au>

       This program is free software; you can redistribute it and/or modify
       it under the terms of the GNU General Public License as published by
       the Free Software Foundation; either version 2 of the License, or
       (at your option) any later version.

       This program is distributed in the hope that it will be useful,
       but WITHOUT ANY WARRANTY; without even the implied warranty of
       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
       GNU General Public License for more details.

       You should have received a copy of the GNU General Public License
       along with this program; if not, write to the Free Software
       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
       MA 02110-1301, USA.

=head1 SYNOPSIS

    use VelvetOpt::hwrap;
    use VelvetOpt::Assembly;
    my $object = VelvetOpt::Assembly->new(
        timestamph => "23 November 2008 15:00:00",
        ass_id => "1",
        versionh => "0.7.04",
        pstringh => "test 21 -fasta test_reads.fna",
        ass_dir => "/home/gla048/Desktop/newVelvetOptimiser/data_1"
    );
    my $worked = VelvetOpt::hwrap::objectVelveth($object);
    if($worked){
        print $object->toString();
    }
    else {
        die "Error in velveth..\n" . $object->toString();
    }

=head1 DESCRIPTION

A wrapper module to run velveth on VelvetAssembly objects or on velveth
parameter strings. Also contains private methods to check velveth
parameter strings, run velveth and return results.

=head2 Uses

=over 8

=item strict

=item warnings

=item Carp

=item VelvetOpt::Assembly

=item POSIX qw(strftime)

=back

=head2 Private Fields

=over 8

=item interested

STDERR printing debug message toggle.  1 for on, 0 for off.

=back

=head2 Methods

=over 8

=item _runVelveth

Private method which runs velveth with the supplied velveth parameter string and returns velveth output messages as a string.

=item _checkVHString

Private method which checks for a correctly formatted velveth string.  Returns 1 or 0.

=item objectVelveth

Accepts a VelvetAssembly object and the number of categories velvet was compiled with, looks for the velveth parameter string it contains, checks it, sends it to _runVelveth, collects the results and stores them in the VelvetAssembly object.

=item stringVelveth

Accepts a velveth parameter string and the number of categories velvet was compiled with, checks it, sends it to _runVelveth and then collects and returns the velveth output messages.

=back

=cut

use warnings;
use strict;
use Carp;
use VelvetOpt::Assembly;
use POSIX qw(strftime);

my $interested = 0;

my @Fileformats;
my @Readtypes;
my $usage;
my $inited = 0;

sub init {
	#run a velveth to get its help lines..
	my $response = &_runVelveth(" ");
	
	$response =~ m/CATEGORIES = (\d+)/;
	my $cats = $1;
	unless($cats){$cats = 2;}
	
	$response =~ m/(File format options:(.*)Read type options)/s;
	my @t = split /\n/, $1;
	foreach(@t){
		#if(/\s+(-\S+)/){
		while(/\s+(-\S+)/g){
			push @Fileformats, $1;
		}
	}
	
	$response =~ m/(Read type options:(.*)Options:)/s;
	
	@t = ();
	@t = split /\n/, $1;
	foreach(@t){
		#if(/\s+(-\S+)/){
		while(/\s+(-\S+)/g){
			push @Readtypes, $1;
		}
	}
	
	for(my $i = 3; $i <= $cats; $i++){
		push @Readtypes, "-short$i";
		push @Readtypes, "-shortPaired$i";
	}
	
	$usage = "Incorrect velveth parameter string: Needs to be of the form\n{[-file_format][-read_type] filename}\n";
	$usage .= "Where:\n\tFile format options:\n";
	foreach(@Fileformats){
		$usage .= "\t$_\n";
	}
	$usage .= "Read type options:\n";
	foreach(@Readtypes){
		$usage .= "\t$_\n";
	}
	$usage .= "\nThere can be more than one filename specified as long as its a different type.\nStopping run\n";
	
	$inited = 1;
}

sub _runVelveth {
	#unless($inited){ &init(); }
    my $cmdline = shift;
    my $output = "";
    print STDERR "About to run velveth!\n" if $interested;
    $output = `velveth $cmdline`;
    $output .= "\nTimestamp: " . strftime("%b %e %Y %H:%M:%S", localtime) . "\n";
    return $output;
}

sub _checkVHString {
    unless($inited){ &init(); }
	my $line = shift;
	my $cats = shift;
	
	
	
	my %fileform = ();
    my %readform = ();
	
	foreach(@Fileformats){ $fileform{$_} = 1;}
    foreach(@Readtypes){ $readform{$_} = 1;}

    my @l = split /\s+/, $line;

    #first check for a directory name as the first parameter...
    my $dir = shift @l;
    if(!($dir =~ /\w+/) || ($dir =~ /^\-/)){
        carp "**** $line\n\tNo directory name specified as first parameter in velveth string. Internal error!\n$usage";
        return 0;
    }
    #print "VH Check passed directory..\n";
    my $hash = shift @l;
    unless($hash =~ /^\d+$/){
        carp "**** $line\n\tHash value in velveth string not a number. Internal error!\n$usage";
        return 0;
    }

    #print "VH check passed hash value..\n";

    my $i = 0;
    my $ok = 1;
    foreach(@l){
        if(/^-/){
            #s/-//;
            if(!$fileform{$_} && !$readform{$_}){
                carp "**** $line\n\tIncorrect fileformat or readformat specified.\n\t$_ is an invalid velveth switch.\n$usage";
                return 0;
            }
            elsif($fileform{$_}){
                if(($i + 1) > $#l){
                    carp "$line\n\tNo filename supplied after file format type $l[$i].\n$usage";
                    return 0;
                }
                if($readform{$l[$i+1]}){
                    if(($i+2) > $#l){
                        carp "$line\n\tNo filename supplied after read format type $l[$i+1].\n$usage";
                        return 0;
                    }
                    if(-e $l[$i+2]){
                        $ok = 1;
                    }
                    else{
                        carp "**** $line\n\tVelveth filename " . $l[$i+2] . " doesn't exist.\n$usage";
                        return 0;
                    }
                }
                elsif (-e $l[$i+1]){
                    $ok = 1;
                }
                else {
                   carp "**** $line\n\tVelveth filename " . $l[$i+1] . " doesn't exist.$usage\n";
                    return 0;
                }
            }
            elsif($readform{$_}){
                if(($i + 1) > $#l){
                    carp "$line\n\tNo filename supplied after read format type $l[$i].\n$usage";
                    return 0;
                }
                if($fileform{$l[$i+1]}){
                    if(($i+2) > $#l){
                        carp "$line\n\tNo filename supplied after file format type $l[$i+1].\n$usage";
                        return 0;
                    }
                    if(-e $l[$i+2]){
                        $ok = 1;
                    }
                    else{
                        carp "**** $line\n\tVelveth filename " . $l[$i+2] . " doesn't exist.\n$usage";
                        return 0;
                    }
                }
                elsif (-e $l[$i+1]){
                    $ok = 1;
                }
                else {
                    carp "**** $line\n\tVelveth filename " . $l[$i+1] ." doesn't exist.\n$usage";
                    return 0;
                }
            }
        }
        elsif(!-e $_){
            carp "**** $line\n\tVelveth filename $_ doesn't exist.\n$usage";
            return 0;
        }
        $i ++;
    }
    if($ok){
        return 1;
    }
}

sub objectVelveth {
    unless($inited){ &init(); }
    my $va = shift;
	my $cats = shift;
    my $cmdline = $va->{pstringh};
    if(_checkVHString($cmdline, $cats)){
        $va->{velvethout} = _runVelveth($cmdline);
        my @t = split /\n/, $va->{velvethout};
        $t[$#t] =~ s/Timestamp:\s+//;
        $va->{timestamph} = $t[$#t];
        return 1;
    }
    else {
        $va->{velvethout} = "Formatting errors in velveth parameter string.$usage";
        return 0;
    }
}

sub stringVelveth {
	unless($inited){ &init(); }
    my $cmdline = shift;
	my $cats = shift;
    if(_checkVHString($cmdline,$cats)){
        return _runVelveth($cmdline);
    }
    else {
        return "Formatting errors in velveth parameter string.$usage";
    }
}

1;
