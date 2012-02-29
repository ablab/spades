###########################################################################
#       This software and its documentation are copyright (2009) by the   #
#   Broad Institute/Massachusetts Institute of Technology.  All rights    #
#   are reserved.  This software is supplied without any warranty or      #
#   guaranteed support whatsoever. Neither the Broad Institute nor MIT    #
#   can be responsible for its use, misuse, or functionality.             #
###########################################################################
# 
# PerlUtils.pm
#
# A repository for useful functions
#
#
# 2010-07-07    Filipe Ribeiro     ribeiro@broadinstitute.org
#

package PerlUtils;

use Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(
&set_verbose
&ISO_date
&abort
&print_verbose
&print_error
&print_warning
&print_action
&print_subaction
&bool_str
&array_from_ref_or_value
&globable_str_from_array_or_value
&str_from_hash_ref
);


use strict;   

$PerlUtils::VERBOSE = 1;

sub set_verbose { $PerlUtils::VERBOSE = (defined $_[0] ? $_[0] : 1); };

sub ISO_date { my ($s, $m, $h, $D, $M, $Y) = localtime();
               return sprintf("%4d-%02d-%02d %02d:%02d:%02d", 1900+$Y, 1+$M, $D, $h, $m, $s); }


sub print_action { if ($PerlUtils::VERBOSE) { print "---- @_\n" ;} }
sub print_subaction { if ($PerlUtils::VERBOSE) { print "     @_\n" ;} }
sub print_error { print "**** @_\n"; }
sub print_warning { print "!!!! @_\n"; }


sub abort { print_error(@_); die "\n"; }



sub bool_str { return $_[0] ? "True" : "False"; }


sub array_from_ref_or_value 
{
    return () unless defined $_[0];
    return @{$_[0]} if ($_[0] =~ /ARRAY/);
    return $_[0];
}

sub globable_str_from_array_or_value
{
    return "" unless @_;
    if (ref $_[0] eq "ARRAY") {
        if (scalar @{$_[0]} == 1) {
            return ${$_[0]}[0];
        }
        return "{" . join(",", @{$_[0]}) . "}";
    }
    if (scalar @_ == 1) {
        return $_[0];
    }
    return "{" . join (",", @_) . "}";
}


sub str_from_hash_ref
{
    if (ref $_[0] eq "HASH") 
    {
        return "( " . join(", ", map("'$_' => '$_[0]->{$_}'", keys %{$_[0]})) . " )";
    }
    return undef;
}






1;
