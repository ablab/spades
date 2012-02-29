package ArachneArgs;

use strict;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK);
use Term::ANSIColor qw(:constants);

require Exporter;

@ISA = qw(Exporter AutoLoader);

@EXPORT = qw( 
getCommandArguments 
printCommandArguments 
moduleArguments 
PrintCommandHeader 
);




sub _getFormattingCharacterMap 
{
    my %fcmap = (
	'section' => '',
	'arg' => '',
	'default' => '',
	'end' => ''
	);

    if (defined($ENV{'ARACHNE_PRETTY_HELP'})) {
	if ($ENV{'ARACHNE_PRETTY_HELP'} eq 'Color') {
	    $fcmap{'section'} = (RED . BOLD);
	    $fcmap{'arg'} = (MAGENTA . BOLD);
	    $fcmap{'default'} = BLUE;
	    $fcmap{'end'} = RESET;
	} elsif ($ENV{'ARACHNE_PRETTY_HELP'} eq 'Bold') {
	    $fcmap{'section'} = BOLD;
	    $fcmap{'arg'} = BOLD;
	    $fcmap{'arg'} = BOLD;
	    $fcmap{'default'} = "";
	    $fcmap{'end'} = RESET;
	}
    }

    return %fcmap;
}




sub _usage 
{
    my ($requiredArgsRef, $helpRef) = @_;
    my %requiredArgs = %$requiredArgsRef;
    my %help = (defined($helpRef)) ? %$helpRef : ();
    my %optionalArgs;
    my %fcmap = &_getFormattingCharacterMap();

    print "\n$fcmap{'section'}Usage: $0 arg1=value1 arg2=value2 ...$fcmap{'end'}\n\n";

    print "$fcmap{'section'}Required arguments:$fcmap{'end'}\n\n";

    foreach my $key (sort { $a cmp $b } keys(%requiredArgs)) {
	next if ($key =~ /_postprocess/ || $key =~ /_preprocess/);
	if (defined($requiredArgs{$key})) { $optionalArgs{$key} = $requiredArgs{$key}; }
	else {
	    print "$fcmap{'arg'}$key$fcmap{'end'}\n";

	    if (defined($help{$key})) {
		print "  $help{$key}\n";
	    }
	}
    }
    print "\n";

    return unless keys(%optionalArgs);

    print "$fcmap{'section'}Optional arguments:$fcmap{'end'}\n\n";

    foreach my $key (sort { $a cmp $b } keys(%optionalArgs)) {
	if (defined($requiredArgs{$key})) {
	    print("$fcmap{'arg'}$key$fcmap{'end'} $fcmap{'default'}default: " . 
		  ((ref($requiredArgs{$key}) eq 'ARRAY') ? 
		   "\"{" . join(",", @{$requiredArgs{$key}}) . "}\"" : 
		   $requiredArgs{$key}) . 
		  "$fcmap{'end'}\n");

	    if (defined($help{$key})) {
		print "  $help{$key}\n";
	    }
	}
    }
    print "\n";
}

# Parse the command-line arguments in Arachne style 
# (including the ability to have whitespace between an equal sign and the parameter.
sub getCommandArguments 
{
    my %requiredArgs = @_;
    my %help;

    # Clean up our required arguments
    foreach my $key (keys(%requiredArgs)) {
	if (ref($requiredArgs{$key}) eq 'HASH') {
	    $help{$key} = ${$requiredArgs{$key}}{'help'};
	    $requiredArgs{$key} = ${$requiredArgs{$key}}{'value'};
	}
	if (defined($requiredArgs{$key})) {
	    $requiredArgs{$key} =~ s/[\r\n]//g;
	    $requiredArgs{$key} = 1 if (lc($requiredArgs{$key}) eq "true");
	    $requiredArgs{$key} = 0 if (lc($requiredArgs{$key}) eq "false");
	}
    }

    # Set our required argument defaults.
    my %args = (
	'NO_HEADER' => 0,
	%requiredArgs,
	);

    if (defined($requiredArgs{'NH'}) && $requiredArgs{'NH'} =~ /(True|true|1)/) {
	$args{'NO_HEADER'} = 1;
	delete($args{'NH'});
	delete($requiredArgs{'NH'});
    }

    # Print usage and exit if we're not supplied with any arguments.
    if ($#ARGV == -1) {
	if (defined($requiredArgs{'_usage'})) 
	{ 
	    &{$requiredArgs{'_usage'}}(\%requiredArgs, \%help); 
	}
	else 
	{ 
	    &_usage(\%requiredArgs, \%help); 
	}
	exit(-1);
    }

    # Clean up the command-line arguments so that we can accept arguments with spaces 
    # and things like 'KEY= VALUE' in addition to the normal 'KEY=VALUE'.
    for (my $i = 0; $i <= $#ARGV; $i++) {
	my $arg = $ARGV[$i];
	if ($arg =~ /\w+=$/) {
	    until (($i+1) > $#ARGV || $ARGV[$i+1] =~ /\w+=/) { $arg .= $ARGV[++$i]; }
	}

	if ($arg =~ /(NO_HEADER|NH)=(\s)?(True|true|1)/) # Turn off automatic banner
	{
	    $args{'NO_HEADER'} = 1;
	} 
	elsif ($arg =~ /(.+)=(.+)/) 
	{
	    my ($key, $value) = ($1, $2);

	    # Store arguments that are of no interest to us in a separate variable.  
	    # This makes it convenient to allow certain arguments to pass through to 
	    # another script by simply appending this extra argument to its command-line.
	    if (!exists($requiredArgs{$key})) {
		$args{'_extra'} .= " $key=\"$value\"";
	    }

	    # Parse boolean values
	    if ($value eq 'True'  || $value eq 'true') 
	    { 
		$args{$key} = 1; 
	    } 
	    elsif ($value eq 'False' || $value eq 'false') 
	    { 
		$args{$key} = 0; 
	    }
	    # Parse array values
	    elsif ($value =~ /{(.+)}/) 
	    { 
                $value =~ s/\s//g;
		my @values = glob($value);
		$args{$key} = \@values;
	    } 
	    # Parse a regular ol' KEY=VALUE pair
	    else 
	    {		
		$args{$key} = $value;
	    }
	} 
	elsif ($arg =~ /(.+)=$/) # Parse a KEY=VALUE pair where VALUE is empty
	{
	    $args{$1} = "";
	} 
	elsif ($arg =~ /-(h|help|\?)/) # Print help
	{
	    if (defined($requiredArgs{'_usage'})) 
	    { 
		&{$requiredArgs{'_usage'}}(\%requiredArgs, \%help); 
	    }
	    else 
	    { 
		&_usage(\%requiredArgs, \%help); 
	    }

	    exit(-1);
	}
    }

    # Pre-process arguments
    if (defined($requiredArgs{'_preprocess'})) 
    {
	&{$requiredArgs{'_preprocess'}}(\%args);
	delete($args{'_preprocess'});
    }

    # Print the header box with info about the command
    &PrintCommandHeader() unless ($args{'NO_HEADER'});

    # Did the user forget any arguments?
    my $missingArgs = 0;
    foreach my $requiredArg (keys(%requiredArgs)) 
    {
	if ($requiredArg !~ /_(pre|post)process/ && !defined($args{$requiredArg})) 
	{
	    print "$requiredArg must be supplied.\n";
	    $missingArgs = 1;
	}
    }

    if ($missingArgs) { die ("Error: some required arguments were not supplied.\n"); }

    # Post-process arguments
    if (defined($requiredArgs{'_postprocess'})) {
	&{$requiredArgs{'_postprocess'}}(\%args);
	delete($args{'_postprocess'});
    }

    # We're all good!
    return %args;
}









# Print our command arguments
sub printCommandArguments 
{
    my ($argsref) = @_;
    my %args = %$argsref;

    foreach my $key (sort { $a cmp $b } keys(%args)) {
	if (ref($args{$key}) eq 'ARRAY') {
	    print "$key => {" . join(",", @{$args{$key}}) . "}\n";
	} else {
	    print "$key => $args{$key}\n";
	}
    }
}



# Returns a hash (by reference) with the required command-line arguments
# for the named C++ module.
# The trick is to run the module with no arguments, prompting the usage message
# (defined in system/ParsedArgs.cc), and then capture and parse this message.
sub moduleArguments($) 
{
    my ($module_name) = @_;
    my %args;
    my ($key, $value);

    # Temporarily setenv ARACHNE_PRETTY_HELP to "Bold" - to help with parsing
    my $temp = $ENV{'ARACHNE_PRETTY_HELP'};
    $ENV{'ARACHNE_PRETTY_HELP'} = "Bold";

    # Escape character - appears in output when ARACHNE_PRETTY_HELP="Bold"
    my $esc = chr(27);
    my $optional = 0;

    open (FH, "$module_name |");

    # Process each line of output into a command-line argument
    foreach my $line (<FH>) {

	$optional = 1 if ($line =~ /Optional arguments/);

	# Match line against the specific stdout format given in ParsedArgs
	next unless ($line =~ /$esc\[01m(.+?)$esc\[0m(.+)/);
	$key = $1;

	# If an argument is optional, but no value is specified in the usage
	# message, this means that the default value is in fact an empty string
	$value = $optional ? "" : undef;

	# Look for a default value in this line
	if ($2 =~ /default\: (.+)$/) {
	    $value = $1;
	}

	$args{$key} = $value;
    }
    close FH;
    $ENV{'ARACHNE_PRETTY_HELP'} = $temp;



    return \%args;
}




# Print the fancy header box with info about the command,
# including arguments supplied to it
# This parallels the function PrintTheCommandPretty in system/ParsedArgs.cc
sub PrintCommandHeader 
{

    my ($fh, $prefix, $thickbar) = (*STDOUT, ''); # default values
    $fh     = $_[0] if ($_[0]); # filehandle to print the header to
    $prefix = $_[1] if ($_[1]); # string to be prepended to every line of the header
    $thickbar = $_[2] if ($_[2]); # print a thicker version of the bar ('=')
    my $width = 80;

    my @stat = stat $0;
    my $mtime = localtime($stat[9]);
    my $now   = localtime();
    # drop the day of the week
    $mtime = substr( $mtime, 4 );
    $now   = substr( $now  , 4 );

    # truncate host name to 4 letters, e.g. 'crd1'
    my $host = (exists $ENV{HOST} ? substr($ENV{HOST}, 0, 4) : "localhost");

    my $bar = $thickbar ? '='x$width : '-'x$width;
    my $timestamp = "$now run on $host (pid=$$), last modified $mtime";

    my $command = $0;
    $command .= ' ' while length $command < $width - 1;
    $command .= "\\";

    # Fill @args_parsed with lines of (parsed) info about the args
    my @args = @ARGV;
    my @args_parsed = ();
    my $line = '    ';
    while (@args) {
	my $arg = shift @args;

        $arg = "$1\"$2\"" if ($arg =~ /^(.+=\s*)(\{.+\})$/);

	# Start a new line, if necessary
	if (length($line) + 1 + length($arg) >= $width - 1 &&
	    $line ne '') {
	    $line .= ' ' while length $line < $width - 1;
	    $line .= "\\";
	    push @args_parsed, $line;
	    $line = '    ';
	}

	$line .= "$arg ";
    }
    if ($line) {
	push @args_parsed, $line;
    }


    # We have prepared the output lines;
    # now, prepend each line with the prefix, and append a newline
    map {$_ = "$prefix$_\n"} ($bar, $timestamp, $command, @args_parsed);

    # Print lines to filehandle
    print $fh $bar;
    print $fh $timestamp;
    print $fh $command;
    print $fh @args_parsed;
    print $fh $bar;
    print $fh "\n";
}


1;
