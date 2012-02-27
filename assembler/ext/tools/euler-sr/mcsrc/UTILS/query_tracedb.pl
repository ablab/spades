#!/usr/bin/perl -w
use strict;
use LWP::UserAgent;
use HTTP::Request::Common 'POST';

$ENV{'LANG'}='C';
$ENV{'LC_ALL'}='C';

my $query = join ' ', @ARGV;
$query = 'help' if $query =~ /^(\-h|\-\-help|\-)$/;
$query = join('', <STDIN>) if ! $query;

my $req = POST 'http://www.ncbi.nlm.nih.gov/Traces/trace.cgi?cmd=raw', [query=>$query];
my $res =  LWP::UserAgent->new->request($req, sub { print $_[0] });
die "Couldn't connect to TRACE server\n" if ! $res->is_success;

# After copying this script, make it executable:
# chmod +x query_tracedb
#
#
# Then apply
#
# ./query_tracedb usage
#
# command to see the help screen with usage examples
