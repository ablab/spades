###########################################################################
#       This software and its documentation are copyright (2009) by the   #
#   Broad Institute/Massachusetts Institute of Technology.  All rights    #
#   are reserved.  This software is supplied without any warranty or      #
#   guaranteed support whatsoever. Neither the Broad Institute nor MIT    #
#   can be responsible for its use, misuse, or functionality.             #
###########################################################################
# 
# NaifDB.pm
#
#  A NaifDB is a hash reference.
#  Each hash key is a field, and each hash value is a reference to an array.
#  All the arrays have the same size by construction.
#
#  $db->{name}   = ["Jerry",    "Elaine",             "George",    "Cosmo"]
#  $db->{job}    = ["Comedian", "Personal Assistant", "Architect", "Unemployed"] 
#  $db->{so}     = ["Dolores",  "Putty",              "Susan",     "Mona" ]
#  $db->{gender} = ["male",     "female",             "male",      "male" ] 
#  $db->{height} = ["175",      "160",                "160",       "185"  ]
#
# 2010-05-26    Filipe Ribeiro     ribeiro@broadinstitute.org
#

package NaifDB;

use Fcntl qw(:flock); # LOCK_EX, LOCK_SH

use strict;

use vars qw(@ISA @EXPORT);

use Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(
&new_naif_db
);



# ------------ list of defined subs

#sub new_naif_db
#sub new 
#sub DESTROY
#sub size
#sub fields
#sub record_add
#sub record_get
#sub record_set
#sub record_erase
#sub sub_db_select
#sub sub_db_select_test
#sub sub_db_select_count
#sub values_frequencies
#sub record_index
#sub record_indices
#sub is_key
#sub field_hash
#sub next_index
#sub next_4d_index
#sub from_csv_or_continue
#sub from_csv_or_die
#sub from_csv
#sub to_csv
#sub validate
#sub _csv_split_line
#sub _csv_compound_line
#sub _fields_validate
#sub _file_backup



# Exported routine to create a new NaifDB.
# 
#  my $db = new_naif_db();
#
sub new_naif_db { return __PACKAGE__->new(@_); }




# Contructor
# If called as $db->new(), creates a new NaifDB with the same fields.
#
#  my $new_db = $db->new()
#  my $new_db = NaifDB->new()
#
sub new 
{ 
    my $that = shift;
    my $class = ref($that) || $that;
    my $self = {};

    if (ref($that)) {
        map $self->{$_} = [], (keys %$that);
    }
    elsif (@_) {
        if (ref $_[0] eq "HASH") {
            map $self->{$_} = [], keys %{$_[0]};
        }
        elsif (ref $_[0] eq "ARRAY") {
            map $self->{$_} = [], @{$_[0]};
        }
        else {
            die "**** Don't know what to do with '$_[0]'.";
        }
    }

    bless $self, $class;
    return $self;
}




# Destructor
sub DESTROY {}



# Returns the number of records in the db.
#
#  my $size = $db->size()
#  
sub size
{
    my $self = shift;
    return 0 unless %$self;
    return scalar @{$self->{(keys %$self)[0]}};
}



# Returns an array with all the fields.
#
#  my @fields = $db->fields();
#
# You get: @fields = ( "name", "job", "so", "gender" )
#
sub fields
{
    my $self = shift;
    return () unless %$self;
    return keys %$self;
}


# Returns an empty record.
#
#  my $rec = $db->record_new();
#
sub record_new
{
    my $self = shift;
    my %rec = ();
    map $rec{$_} = "", keys %$self;
    return \%rec;
}


# Adds a record to the data base.
#
#  $db->record_add($record);
#
sub record_add
{
    my $self = shift;
    my ($record) = @_;

    

    if ($self->size()) 
    {
        foreach (keys %$self) 
        {
            die "**** missing key '$_' in hash being pushed back"
                unless exists $record->{$_};
        }
        foreach (keys %$record) 
        {
            die "**** extra key '$_' in hash being pushed back"
                unless exists $self->{$_};
        }
        
        map { push @{$self->{$_}}, $record->{$_} } (keys %$record);
    }
    else 
    {
        map { $self->{$_} = [ $record->{$_} ] } (keys %$record);
    }
}



# Retrives a record at the specified index.
#
#  my $rec = $db->record_get(2);
#
# You get: $rec = { name => "George", job => "Architect", so => "Susan", gender => "male" }
#
sub record_get
{
    my $self = shift;
    my ($i) = @_;
    
    if ($i > $self->size()) {
        my $n = $self->size();
        print STDERR "**** index $i out of range ($n)\n";
        exit(1);
    }

    my %record = ();
    map $record{$_} = $self->{$_}[$i], keys %$self;

    return \%record;
}



# Sets a record at the specified index.
#
#  $db->record_set($i, $rec);
#
sub record_set
{
    my $self = shift;
    my ($i, $record) = @_;
    die "**** record not defined"
        unless defined $record;
    die "**** record out of range"
        unless ($i < $self->size());

    foreach (keys %$self) 
    {
        die "**** missing key '$_' in input record."
            unless exists $record->{$_};
    }
    foreach (keys %$record) 
    {
        die "**** extra key '$_' in input record."
            unless exists $self->{$_};
    }

    map { $self->{$_}[$i] = $record->{$_} } (keys %$record);

    return $record;
}


# Erases a record at the specified index.
#
#  $db->record_erase($i);
#
sub record_erase
{
    my $self = shift;
    my ($i) = @_;
    die "**** record out of range"
        unless ($i < $self->size());

    map { splice @$_, $i, 1 } (values %$self);
}






# Retrieves a sub database with all the records with matching keys.
#
#  my $sub_db = $db->sub_db_select("gender", "male");
#
# You get: $sub_db = { name => ["Jerry", "George", "Cosmo" ], job => ... }
#
#  my $sub_db = $db->sub_db_select("so", "Putty", "Mona");
#
# You get: $sub_db = { name => ["Elaine", "Cosmo" ], job => ... }
#
sub sub_db_select
{
    my $self = shift;
    my $n = $self->size();
    return $self if ($n == 0);

    my ($key, @selected_key_values) = @_; 
    die "**** invalid key '$key'" 
        unless (defined $key and exists $self->{$key});

    # no key values, return all records
    return $self 
        unless (@selected_key_values && $selected_key_values[0]); 
    
    # make a hash of selected_key_values because is easily searchable with 'exists'
    my %selected_key_values = ();
    foreach my $kv (@selected_key_values) {
        if (ref($kv)) { # assume it's an array reference
            map $selected_key_values{$_} = 1, @$kv;
        }
        else { # it's a scalar
            $selected_key_values{$kv} = 1;
        }
    }

    my @keys = keys %$self;
    
    my $sub_db = $self->new();
    for (my $i = 0; $i != $n; $i++) 
    {
        if (exists $selected_key_values{$self->{$key}[$i]}) 
        {
            map { push @{$sub_db->{$_}}, $self->{$_}[$i] } @keys;
        }
    }
    return $sub_db;
}





# Retrieves a sub database with all the records that satisfy a test
#
#  my $sub_db = $db->sub_db_select_test("height", { return $_[0] < 170; });
#
# You get: $sub_db = { name => ["George", "Elaine" ], job => ... }
#
sub sub_db_select_test
{
    my $self = shift;
    my $n = $self->size();
    return $self if ($n == 0);

    my ($key, $test) = @_; 
    die "**** invalid key '$key'" 
        unless (defined $key and exists $self->{$key});

    # no test, return all records
    return $self 
        unless (defined $test);
    
    my @keys = keys %$self;
    
    my $sub_db = $self->new();
    for (my $i = 0; $i != $n; $i++) 
    {
        if ($test->($self->{$key}[$i]))
        {
            map { push @{$sub_db->{$_}}, $self->{$_}[$i] } @keys;
        }
    }
    return $sub_db;
}







# Counts the number of records with matching keys.
#
#  my $n = $db->sub_db_select_count("gender", "female", "male");
#
# You get: $n = 4
# 
sub sub_db_select_count
{
    my $self = shift;
    my $n = $self->size();
    return 0 if ($n == 0);

    my ($key, @selected_key_values) = @_; 
    die "**** invalid key '$key'"
        unless (defined $key and exists $self->{$key});
    
    # no key values, return size
    return $self->size() 
        unless (@selected_key_values && $selected_key_values[0]);

    # make a hash of selected_key_values because is easily searchable with 'exists'
    my %selected_key_values = ();
    foreach my $kv (@selected_key_values) {
        if (ref($kv)) { # assume it's an array reference
            map $selected_key_values{$_} = 1, @$kv;
        }
        else { # it's a scalar
            $selected_key_values{$kv} = 1;
        }
    }

    my @keys = keys %$self;
    
    my $count = 0;
    for (my $i = 0; $i != $n; $i++) {
        $count++ if (exists $selected_key_values{$self->{$key}[$i]});
    }
    return $count;
}









# Gathers the frequencies of all key values for a specific key.
#
#  my $val_freq = $db->values_frequencies("gender");
#
# You get: $val_freq = { male => 3, female => 1 }
# 
sub values_frequencies
{
    my $self = shift;
    return {} unless $self->size();

    my ($key) = @_; 

    my %values_freq = ();
    
    my $selected_key_values = $self->{$key};
    
    foreach my $value (@$selected_key_values) {
        if (exists $values_freq{$value}) {
            $values_freq{$value} ++;
        }
        else {
            $values_freq{$value} = 0;
        }
    }
    return \%values_freq;
}


# Verify if a field is a key (i.e., there is only one record for each field value)
#
sub is_key
{
    my $self = shift;
    return 1 unless $self->size();

    my ($key) = @_; 

    my %value_found = ();
    
    my $n = $self->size();

    for (my $i = 0; $i != $n; $i++) {
        
        my $value = $self->{$key}[$i];
        return 0 if (exists $value_found{$value});
        
        $value_found{$value} = 1;
    }
    return 1;
}


# Given a field assumed to be a key (i.e., there is only one record for each field value)
# return a hash reference with the record indices.
#
#  my $rec_i = $db->record_indices("name");
#
# You get: $rec_i = { "Elaine" => 2, "Cosmo" => 3, "Jerry" => 0, "George" => 1 }
# 
sub record_indices
{
    my $self = shift;
    return {} unless $self->size();

    my ($key) = @_; 

    my %rec_ind = ();
    
    my $n = $self->size();

    for (my $i = 0; $i != $n; $i++) {
        
        my $value = $self->{$key}[$i];
        die "**** Field '$key' refers to more than one record"
            if (exists $rec_ind{$value});
        
        $rec_ind{$value} = $i;
    }
    return \%rec_ind;
}





# Given a field assumed to be a key (i.e., there is only one record for each field value)
# return a hash reference with the record indices.
#
#  my $rec_i = $db->record_index("name", "Elaine");
#
# You get: $rec_i = 2;
# 
sub record_index
{
    my $self = shift;
    return {} unless $self->size();

    my ($key, $value) = @_; 

    my $rec_ind = -1;
    
    my $n = $self->size();

    for (my $i = 0; $i != $n; $i++) {
        
        my $v = $self->{$key}[$i];
        if ($v eq $value) {
            die "**** Field '$key' refers to more than one record"
                if ($rec_ind >= 0);
            $rec_ind = $i;
        }
    }
    return $rec_ind;
}





# Given a field assumed to be a key (i.e., there is only one record for each field value)
# return a hash reference with the values of the second field.
#
#  my $rec_field = $db->field_hash("name", "gender");
# 
# You get: $rec_field = { "Elaine" => "female", "Cosmo" => "male", 
#                         "Jerry" => "male", "George" => "male" }
#
sub field_hash
{
    my $self = shift;
    return {} unless $self->size();
    
    my ($key, $field) = @_; 
    
    my %rec_field = ();
    
    my $n = $self->size();
    
    for (my $i = 0; $i != $n; $i++) {
        
        my $value = $self->{$key}[$i];
        die "**** Field '$key' refers to more than one record"
            if (exists $rec_field{$value});
        
        $rec_field{$value} = $field;
    }
    return \%rec_field;
}
    






sub next_index
{
    my $self = shift;
    my ($field) = @_;

    die "**** You must specify a field"
        unless (defined $field);

    die "**** '$field' is not a field"
        unless (exists $self->{$field});

    my $index = 1 + ($self->size() ? $self->{$field}[-1] : 0);

    die "**** No room in database for new entry (limit is 10^4).\n" 
        if ( $index >= 1e4 );

    return $index;
}
 

sub next_4d_index
{
    my $self = shift;
    my ($field) = @_;
    return sprintf("%04d", $self->next_index($field));
}





# Reads a db from a comma-separated-value file.
# Returns an empty db if the file is not found.
#
#  $db->from_csv_or_continue("mydb.csv");
#
sub from_csv_or_continue
{
    my $self = shift;
    my ($fn, $valid_fields) = @_;

    return $self = new_naif_db($valid_fields) unless (-s $fn);
    
    return $self->from_csv($fn, $valid_fields);
}






# Reads a db from a comma-separated-value file.
# Dies if file is not found.
#
#  $db->from_csv_or_die("mydb.csv");
#
sub from_csv_or_die
{
    my $self = shift;
    my ($fn, $valid_fields) = @_;

    die "**** File '$fn' does not exist." unless (-s $fn);

    return $self->from_csv($fn, $valid_fields);
}







# Reads a db from a comma-separated-value file.
#
#  $db->from_csv("mydb.csv");
#
#  $db->from_csv("mydb.csv", { name   => "string", 
#                              job    => "string", 
#                              so     => "string", 
#                              gender => { "male" => 1, "female" => 1 },
#                              height => "number" } );
#
sub from_csv
{
    my $self = shift;
    my ($fn, $valid_fields) = @_;

    my $file;
    open($file, "<$fn") or die "**** can't open '$fn' for reading";
    flock($file, LOCK_SH);  # get a shared lock before reading (prevents writing)
    my @lines = <$file>;
    close($file);

    return $self->from_csv_lines(\@lines, $valid_fields);
}








# Writes a db to a comma-separated-value file
#
#  $db->to_csv("mydb.csv", [\@fields]);
#
sub to_csv
{
    my $self = shift;
    my ($fn, $fields, $separator) = @_;

    my $file;
    if ($fn eq "-") {
        open($file, ">&STDOUT") or die "**** can't open '$fn' for writing";
    }
    else {
        _file_backup($fn);
        open($file, ">$fn") or die "**** can't open '$fn' for writing";
        flock($file, LOCK_EX); # get an exclusive lock before writing
    }        
    
    print $file @{$self->to_csv_lines($fields, $separator)};
    close($file);
    return $self;
}







# The next to functions should be used in conjunction when updating an existing database.
# 
#  from_csv_lock_exclusive() 
#    opens a file and reads it, returning the open file handle with an exclusive flock.
#
#  to_csv_lock_exclusive()
#    writes the database to the file handle and closes the file.
#
sub from_csv_lock_exclusive
{
    my $self = shift;
    my ($fn, $valid_fields) = @_;

    _file_backup($fn);
    my $fh;
    system("touch", "$fn");
    open($fh, "+<$fn") or die "**** can't open '$fn' for updating";
    flock($fh, LOCK_EX);  # get an exclusive lock before reading and later writing
    my @lines = <$fh>;    # real entire file

    $self->from_csv_lines(\@lines, $valid_fields);
    return $fh;           # return the file handle for later updating
}

sub to_csv_lock_exclusive
{
    my $self = shift;
    my ($fh, $fields, $separator) = @_;
    
    seek($fh, 0, 0);                                       # go to the begining
    print $fh @{$self->to_csv_lines($fields, $separator)}; # update data
    truncate($fh, tell($fh));                              # truncate at the current position (end)
    close($fh);                                            # close, which releases lock
    return $self;
}










# Reads a db from a reference to a set of comma-separated-value strings.
#
#  $db->from_csv_lines($lines);
#
#  $db->from_csv_lines($lines, { name   => "string", 
#                                job    => "string", 
#                                so     => "string", 
#                                gender => { "male" => 1, "female" => 1 },
#                                height => "number" } );
#
sub from_csv_lines
{
    my $self = shift;
    my ($lines, $valid_fields) = @_;

    my $nl = scalar @$lines;

    my $il = 0;
    my $in_fields = [];

    # ---- Obtain field names from first line
    while (!@$in_fields && $il != $nl) { 
        my $line = $lines->[$il++];
        $in_fields = _csv_split_line($line);
        #print "@$in_fields\n";
    }
    my $n_fields = @$in_fields;
    
    # ---- Initialize the database fields
    map $self->{$_} = [], @$in_fields;

    
    # ---- Add entries to the database 
    while ($il != $nl) {
        my $line = $lines->[$il++];
        my $vals = _csv_split_line($line);
    
        if (@$vals) {
            die "**** line $il has wrong number of values" 
                if (@$vals != $n_fields);

            map { push @{$self->{$in_fields->[$_]}}, $vals->[$_] } (0..$n_fields-1);
        }
    }

    $self->validate($valid_fields) if ($valid_fields);

    return $self;
}




sub to_csv_lines
{
    my $self = shift;
    my ($fields, $sep) = @_;
    my @fields = (defined $fields && ref $fields eq "ARRAY") ? @$fields : keys %$self;

    # ---- Figure out the character length of each field
    my @lens = ();
    my $len = 0;
    foreach my $k (@fields) {
        my $l = length($k);
        map { $l = length($_) if ($_ && length($_) > $l) } @{$self->{$k}}; 
        push @lens, $l;
        $len += $l;
    }
    $len += 2 * (@lens - 1);
    my $separator = (defined $sep) ? $sep x (1 + ($len - 1) / length($sep)) . "\n" : "";

    my @lines = ();


    push @lines, $separator if (defined $sep);
    
    push @lines, (_csv_compound_line(\@fields, \@lens));
    my $n = $self->size();
    for (my $i = 0; $i != $n; $i++) 
    {
        push @lines, _csv_compound_line([ map $self->{$_}[$i], @fields ], \@lens);
    }

    push @lines, $separator if (defined $separator);


    return \@lines;
}






sub validate
{    
    my ($self, $fields) = @_;

    # ---- validate field NAMES

    _fields_validate($self, $fields)
	if (defined $fields);



    # ---- validate field VALUES

    my $nl = $self->size();
    
    my $n_errors = 0;
    foreach my $field (keys %$self) 
    {
	my $type = $fields->{$field};
	for (my $il = 0; $il != $nl; $il++) 
	{
	    my $v = $self->{$field}[$il];
	    if ($type =~ "^string") {
		# nothing to do 		
	    } 
	    elsif ($type =~ "^safe_string") {
                if ($v !~ /^[\w\.\-\+=]+$/) {
		    $n_errors++;
		    printf STDERR "**** line '$il', field '$field': '$v' not a safe string (only 'azAZ09.-+=').\n";
                }
	    } 
	    elsif ($type =~ "^number") {
		if ($v ne "" and
                    $v !~ /^[-\+]?\d*\.?\d+([eE][-\+]?\d+)?$/) 
                {
		    $n_errors++;
		    printf STDERR "**** line '$il', field '$field': '$v' is not a $type.\n";
		}
	    }
	    elsif ($type eq "bool") {
		if (lc $v ne "true" and 
		    lc $v ne "false" and
		    $v ne "1" and 
		    $v ne "0") 
		{
		    $n_errors++;
		    printf STDERR "**** line '$il', field '$field': '$v' is not a $type.\n";
		}
	    }
	    elsif (ref($type)) {
		my $vs = join ", ", (map "'$_'", keys %$type);
		if (! exists $type->{$v}) 
		{
		    $n_errors++;
		    printf STDERR "**** line '$il', field '$field': '$v' is not one of ($vs).\n";
		}
	    }
	    else {
                $n_errors++;
                printf STDERR "**** Don't know how to handle field type '$field'.\n";
	    }
	}
        
    }
    
    
    # ---- validate fields that are KEYS
    
    foreach my $field (keys %$fields) 
    {
        if ($field =~ /_key$/) {
            my %freqs = ();
            
            foreach my $val (@{$self->{$field}}) 
            {
                $freqs{$val} = 0 if (!exists $freqs{$val});
                $freqs{$val}++;
            }
            
            
            foreach my $val (keys %freqs) 
            {
                if ($freqs{$val} > 1) 
                {
                    $n_errors++;
                    printf STDERR "**** value '$val' for key field '$field' shows up $freqs{$val} times which is more than once.\n";
                }
            }
        }
    }



    die "**** Found $n_errors errors while validating database entries.\n" 
        if $n_errors;

}








# ---------------------------------
#  internal functions; not methods
# ---------------------------------

sub _csv_split_line
{
    my ($line) = @_;
    if ($line && $line !~ /^\s*$/) {
        my @vs = split ",", $line;
        map s/^\s+//, @vs;
        map s/\s+$//, @vs;
        map s/\s+/ /g, @vs;
        return \@vs;
    }
    return [];
}


sub _csv_compound_line
{
    my ($vals, $lens) = @_;
    #print "@$vals\n";

    # ABSOLUTELY no commas in values 
    map s/,//g, @$vals;

    my $s = "";
    if (defined $lens) 
    {
        return (join ", ", map sprintf("%$lens->[$_]s", $vals->[$_]), (0..@$vals-1)) . "\n";
        #return sprintf("$s%$lens->[-1]s\n", $vals->[-1]);
    }
    else 
    {
        map $s .= "$vals->[$_], ", (0..@$vals-2);
        return "$s".$vals->[-1]."\n";
    }
}



sub _fields_validate
{
    my ($in_fs, $fs) = @_;

    my $n_errors = 0;

    foreach my $f (keys %$fs) {
        if (! exists $in_fs->{$f}) {
            print "**** Can't find field '$f' in database.\n";
            $n_errors++;
        }
    }

    foreach my $f (keys %$in_fs) {
        if (! exists $fs->{$f}) {
            print "**** Found extra field '$f' in database.\n";
            $n_errors++;
        }
    }

    die "**** Found $n_errors errors while validating database fields.\n" 
        if $n_errors;
}




sub _file_backup
{
    my ($full_fn) = @_;

    if (-s $full_fn) 
    {
        my ($dn, $fn) = ($full_fn =~ /\// ? $full_fn =~ /^(.+)\/([^\/]+)$/ : (".", $full_fn));
        $dn .= "/backup";
        
        #print "dn = $dn\n";
        
        mkdir $dn unless (-d $dn);
        
        my ($s, $m, $h, $D, $M, $Y) = localtime();
        my $full_fn2 = sprintf("$dn/%4d-%02d-%02d.%02d-%02d-%02d.%04d.$fn", 
                               1900+$Y, 1+$M, $D, $h, $m, $s, rand(10000));

        system("cp", "$full_fn", "$full_fn2");
    }
}



1;
