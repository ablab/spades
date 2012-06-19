#       VelvetOpt::gwrap.pm
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
package VelvetOpt::gwrap;

=head1 NAME

VelvetOpt::gwrap.pm - Velvet graphing and assembly program wrapper module.

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

    use VelvetOpt::gwrap;
    use VelvetOpt::Assembly;
    my $object = VelvetOpt::Assembly->new(
        timestamph => "23 November 2008 15:00:00",
        ass_id => "1",
        versiong => "0.7.19",
        pstringg => "test",
        ass_dir => "/home/gla048/Desktop/newVelvetOptimiser/test"
    );
    my $worked = VelvetOpt::gwrap::objectVelvetg($object);
    if($worked){
        print $object->toString();
    }
    else {
        die "Error in velvetg..\n" . $object->toString();
    }

=head1 DESCRIPTION

A wrapper module to run velvetg on VelvetAssembly objects or on velvetg
parameter strings. Also contains private methods to check velvetg
parameter strings, run velvetg and return results.

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

=item _runVelvetg

Private method which runs velvetg with the supplied velvetg parameter string and returns velvetg output messages as a string.

=item _checkVGString

Private method which checks for a correctly formatted velvetg parameter string.  Returns 1 or 0.

=item objectVelvetg

Accepts a VelvetAssembly object, looks for the velvetg parameter string it contains, checks it, sends it to _runVelvetg, collects the results and stores them in the VelvetAssembly object.

=item stringVelvetg

Accepts a velvetg parameter string, checks it, sends it to _runVelvetg and then collects and returns the velvetg output messages.

=back

=cut

use warnings;
use strict;
use Carp;
use VelvetOpt::Assembly;
use POSIX qw(strftime);

my $interested = 0;

sub _runVelvetg {
    my $cmdline = shift;
    my $output = "";
    print STDERR "About to run velvetg!\n" if $interested;
    $output = `velvetg $cmdline`;
    $output .= "\nTimestamp: " . strftime("%b %e %Y %H:%M:%S", localtime) . "\n";
    return $output;
}

sub _checkVGString {
    return 1;
}

sub objectVelvetg {
    my $va = shift;
    my $cmdline = $va->{pstringg};
    if(_checkVGString($cmdline)){
        $va->{velvetgout} = _runVelvetg($cmdline);
        my @t = split /\n/, $va->{velvetgout};
        $t[$#t] =~ s/Timestamp:\s+//;
        $va->{timestampg} = $t[$#t];
        return 1;
    }
    else {
        $va->{velvetgout} = "Formatting errors in velvetg parameter string.";
        return 0;
    }
}

sub stringVelvetg {
    my $cmdline = shift;
    if(_checkVGString($cmdline)){
        return _runVelvetg($cmdline);
    }
    else {
        return "Formatting errors in velvetg parameter string.";
    }
}

1;
