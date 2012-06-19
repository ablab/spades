#       VelvetOpt::Assembly.pm
#
#       Copyright 2008,2009 Simon Gladman <simon.gladman@csiro.au>
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

#	Version 2.1.4
#
#	Changes for 2.0.1
#	*Bug fix in CalcAssemblyScore.  Now returns 0 if there is no calculable score instead of crashing.
#
#	Changes for 2.1.0
#	*Added 2 stage optimisation functions for optimising kmer size and cov_cutoff differently if required.
#
#	Changes for 2.1.1
#	*Allowed for non-word characters in prefix names.  (. - etc.)  Still no spaces allowed in prefix name or any filenames.
#
#	Changes for 2.1.2
#	*Now warns nicely of optimisation function returning undef or 0. Suggests you choose and alternative.
#
#	Changes for 2.1.3
#	*Now prints the velvet calculated insert sizes and standard deviations in the Assembly summaries (both log files and screen).
#
#	Changes for 2.1.4
#	*Fixed a bug where newer versions of velvet would cause the paired end library stats not to be displayed.

package VelvetOpt::Assembly;

=head1 NAME

VelvetOpt::Assembly.pm - Velvet assembly container class.

=head1 AUTHOR

Simon Gladman, CSIRO, 2007, 2008.

=head1 LICENSE

Copyright 2008, 2009 Simon Gladman <simon.gladman@csiro.au>

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

    use VelvetOpt::Assembly;
    my $object = VelvetOpt::Assembly->new(
        timestamph => "23 November 2008 15:00:00",
        ass_id => "1",
        versionh => "0.7.04",
        ass_dir => "/home/gla048/Desktop/newVelvetOptimiser/data_1"
    );
    print $object->toString();

=head1 DESCRIPTION

A container class to hold the results of a Velvet assembly.  Includes timestamps,
version information, parameter strings and assembly output metrics.

Version 2.1.4

=head2 Uses

=over 8

=item strict

=item warnings

=item Carp

=back

=head2 Fields

=over 8

=item assmscore

The assembly score metric for this object

=item timstamph

The timestamp of the start of the velveth run for this assembly

=item timestampg

The date and time of the end of the velvetg run.

=item ass_id

The assembly id number.  Sequential for all the runs for this optimisation.

=item versionh

The version number of velveth used in this assembly

=item versiong

The version number of velvetg used in this assembly

=item readfilename

The name of the file containing all the reads (or a qw of them if more than one...)

=item pstringh

The velveth parameter string used in this assembly

=item pstringg

The velvetg parameter string used in this assembly

=item ass_dir

The assembly directory path (full)

=item hashval

The hash value used for this assembly

=item rmapfs

The roadmap file size

=item sequences

The total number of sequences in the input files

=item nconts

The number of contigs in the final assembly

=item totalbp

The total number of bases in the contigs

=item n50

The n50 of the assembly

=item maxlength

The length of the longest contig in the assembly

=item maxcont

The size of the largest contig in the assembly

=item nconts1k

The number of contigs greater than 1k in size

=item totalbp1k

the sum of the length of contigs > 1k in size

=item velvethout

The velveth output

=item velvetgout

The velvetg output

=back

=head2 Methods

=over 8

=item new

Returns a new VelvetAssembly object.

=item accessor methods

Accessor methods for all fields.

=item calcAssemblyScore

Calculates the assembly score of the object (after velvetg has been run.) and stores it in self.

=item getHashingDetails

Gets the details of the outputs from the velveth run and stores it in self.

=item getAssemblyDetails

Gets the details of the outputs from the velvetg run and stores it in self.

=item toString

Returns a string representation of the object's contents.

=item toStringNoV

Returns a string representation of the object's contents without the velvet outputs which are large.

=item opt_func_toString

Returns the usage of the optimisation function.

=back

=cut

use strict;
use Carp;
use warnings;
#use base "Storable";
use Cwd;
use Bio::SeqIO;

my $interested = 0;



#constructor
sub new {
    my $class = shift;
    my $self = {@_};
    bless ($self, $class);
    return $self;
}

#optimisation function options...
my %f_opts;
	$f_opts{'ncon'}->{'intname'} = 'nconts';
	$f_opts{'ncon'}->{'desc'} = "The total number of contigs";
	$f_opts{'n50'}->{'intname'} = 'n50';
	$f_opts{'n50'}->{'desc'} = "The n50";
	$f_opts{'max'}->{'intname'} = 'maxlength';
	$f_opts{'max'}->{'desc'} = "The length of the longest contig";
	$f_opts{'Lcon'}->{'intname'} = 'nconts1k';
	$f_opts{'Lcon'}->{'desc'} = "The number of large contigs";
	$f_opts{'tbp'}->{'intname'} = 'totalbp';
	$f_opts{'tbp'}->{'desc'} = "The total number of basepairs in contigs";
	$f_opts{'Lbp'}->{'intname'} = 'totalbp1k';
	$f_opts{'Lbp'}->{'desc'} = "The total number of base pairs in large contigs";
	$f_opts{'LNbp'}->{'intname'} = 'numNs1k';
	$f_opts{'LNbp'}->{'desc'} = "The total number of Ns in large contigs";

#accessor methods
sub assmscore{ $_[0]->{assmscore}=$_[1] if defined $_[1]; $_[0]->{assmscore}}
sub timestamph{ $_[0]->{timestamph}=$_[1] if defined $_[1]; $_[0]->{timestamph}}
sub timestampg{ $_[0]->{timestampg}=$_[1] if defined $_[1]; $_[0]->{timestampg}}
sub ass_id{ $_[0]->{ass_id}=$_[1] if defined $_[1]; $_[0]->{ass_id}}
sub versionh{ $_[0]->{versionh}=$_[1] if defined $_[1]; $_[0]->{versionh}}
sub versiong{ $_[0]->{versiong}=$_[1] if defined $_[1]; $_[0]->{versiong}}
sub readfilename{ $_[0]->{readfilename}=$_[1] if defined $_[1]; $_[0]->{readfilename}}
sub pstringh{ $_[0]->{pstringh}=$_[1] if defined $_[1]; $_[0]->{pstringh}}
sub pstringg{ $_[0]->{pstringg}=$_[1] if defined $_[1]; $_[0]->{pstringg}}
sub ass_dir{ $_[0]->{ass_dir}=$_[1] if defined $_[1]; $_[0]->{ass_dir}}
sub hashval{ $_[0]->{hashval}=$_[1] if defined $_[1]; $_[0]->{hashval}}
sub rmapfs{ $_[0]->{rmapfs}=$_[1] if defined $_[1]; $_[0]->{rmapfs}}
sub nconts{ $_[0]->{nconts}=$_[1] if defined $_[1]; $_[0]->{nconts}}
sub n50{ $_[0]->{n50}=$_[1] if defined $_[1]; $_[0]->{n50}}
sub maxlength{ $_[0]->{maxlength}=$_[1] if defined $_[1]; $_[0]->{maxlength}}
sub nconts1k{ $_[0]->{nconts1k}=$_[1] if defined $_[1]; $_[0]->{nconts1k}}
sub totalbp{ $_[0]->{totalbp}=$_[1] if defined $_[1]; $_[0]->{totalbp}}
sub totalbp1k{ $_[0]->{totalbp1k}=$_[1] if defined $_[1]; $_[0]->{totalbp1k}}
sub numNs1k{ $_[0]->{numNs1k}=$_[1] if defined $_[1]; $_[0]->{numNs1k}}
sub velvethout{ $_[0]->{velvethout}=$_[1] if defined $_[1]; $_[0]->{velvethout}}
sub velvetgout{ $_[0]->{velvetgout}=$_[1] if defined $_[1]; $_[0]->{velvetgout}}
sub sequences{ $_[0]->{sequences}=$_[1] if defined $_[1]; $_[0]->{sequences}}
sub assmfunc{ $_[0]->{assmfunc}=$_[1] if defined $_[1]; $_[0]->{assmfunc}}
sub assmfunc2{ $_[0]->{assmfunc2}=$_[1] if defined $_[1]; $_[0]->{assmfunc2}}

#assemblyScoreCalculator
sub calcAssemblyScore {
    use Safe;
	
	my $self = shift;
	my $func = shift;
	
	my $cpt = new Safe;
	
	#Basic variable IO and traversal
	$cpt->permit(qw(null scalar const padany lineseq leaveeval rv2sv rv2hv helem hslice each values keys exists delete rv2cv));
	#Comparators
	$cpt->permit(qw(lt i_lt gt i_gt le i_le ge i_ge eq i_eq ne i_ne ncmp i_ncmp slt sgt sle sge seq sne scmp));
	#Base math
	$cpt->permit(qw(preinc i_preinc predec i_predec postinc i_postinc postdec i_postdec int hex oct abs pow multiply i_multiply divide i_divide modulo i_modulo add i_add subtract i_subtract));
	#Binary math
	$cpt->permit(qw(left_shift right_shift bit_and bit_xor bit_or negate i_negate not complement));
	#Regex
	$cpt->permit(qw(match split qr));
	#Conditionals
	$cpt->permit(qw(cond_expr flip flop andassign orassign and or xor));
	#Advanced math
	$cpt->permit(qw(atan2 sin cos exp log sqrt rand srand));

	foreach my $key (keys %f_opts){
		print "\nkey: $key\tintname: ", $f_opts{$key}->{'intname'}, "\n" if $interested;
		
		$func =~ s/\b$key\b/$self->{$f_opts{$key}->{'intname'}}/g;
	}
		
	my $r = $cpt->reval($func);
	warn $@ if $@;
	$self->{assmscore} = $r;
	unless($r =~ /^\d+/){ 
		warn "Optimisation function did not return a single float.\nOptimisation function was not evaluatable.\nOptfunc: $func";
		warn "Setting assembly score to 0\n"; 
		$self->{assmscore} = 0;
	}
	if($r == 0){
		print STDERR "**********\n";
		print STDERR "Warning: Assembly score for assembly_id " . $self->{ass_id} .  " is 0\n";
		print STDERR "You may want to consider choosing a different optimisation variable or function.\n";
		print STDERR "Current optimisation functions are ", $self->{assmfunc}, " for k value and ", $self->{assmfunc2}, " for cov_cutoff\n";
		print STDERR "**********\n";
	}
	return 1;
}

#getHashingDetails
sub getHashingDetails {
    my $self = shift;
    unless(!$self->timestamph || !$self->pstringh){
        my $programPath = cwd;
        $self->pstringh =~ /^(\S+)\s+(\d+)\s+(.*)$/;
        $self->{ass_dir} = $programPath . "/" . $1;
        $self->{rmapfs} = -s $self->ass_dir . "/Roadmaps";
        $self->{hashval} = $2;
        $self->{readfilename} = $3;
        my @t = split /\n/, $self->velvethout;
        foreach(@t){
            if(/^(\d+).*total\.$/){
                $self->{sequences} = $1;
                last;
            }
        }
        return 1;
    }
    return 0;
}

#getAssemblyDetails
sub getAssemblyDetails {
    my $self = shift;
	my $file = $self->ass_dir . "/contigs.fa";
    unless(!(-e $file)){
		
		my $all = &contigStats($file,1);
		my $large = &contigStats($file,1000);
		
		$self->{nconts} = defined $all->{numSeqs} ? $all->{numSeqs} : 0;
		$self->{n50} = defined $all->{n50} ? $all->{n50} : 0;
		$self->{maxlength} = defined $all->{maxLen} ? $all->{maxLen} : 0;
		$self->{nconts1k} = defined $large->{numSeqs} ? $large->{numSeqs} : 0;
		$self->{totalbp} = defined $all->{numBases} ? $all->{numBases} : 0;
		$self->{totalbp1k} = defined $large->{numBases} ? $large->{numBases} : 0;
		$self->{numNs1k} = defined $large->{numNs} ? $large->{numNs} : 0;
		
		if($self->pstringg =~ m/cov_cutoff/){
			$self->calcAssemblyScore($self->{assmfunc2});
		}
		else {
			$self->calcAssemblyScore($self->{assmfunc});
		}

        return 1;
	}
    return 0;
}

#contigStats
#Original script fa-show.pl by Torsten Seemann (Monash University, Melbourne, Australia)
#Modified by Simon Gladman to suit.
sub contigStats {
	
	my $file = shift;
	my $minsize = shift;
	
	print "In contigStats with $file, $minsize\n" if $interested;
	
	my $numseq=0;
	my $avglen=0;
	my $minlen=1E9;
	my $maxlen=0;
	my @len;
	my $toosmall=0;
	my $nn=0;
	
	my $in = Bio::SeqIO->new(-file => $file, -format => 'Fasta');
	while(my $seq = $in->next_seq()){
		my $L = $seq->length;
		#check > minsize
		if($L < $minsize){
			$toosmall ++;
			next;
		}
		#count Ns
		my $s = $seq->seq;
		my $n = $s =~ s/N/N/gi;
		$n ||= 0;
		$nn += $n;
		#count seqs and other stats
		$numseq ++;
		$avglen += $L;
		$maxlen = $L if $L > $maxlen;
		$minlen = $L if $L < $minlen;
		push @len, $L;
	}
	@len = sort { $a <=> $b } @len;
	my $cum = 0;
	my $n50 = 0;
	for my $i (0 .. $#len){
		$cum += $len[$i];
		if($cum >= $avglen/2) {
			$n50 = $len[$i];
			last;
		}
	}
	
	my %out;
	if($numseq > 0){
		$out{numSeqs} = $numseq;
		$out{numBases} = $avglen;
		$out{numOK} = ($avglen - $nn);
		$out{numNs} = $nn;
		$out{minLen} = $minlen;
		$out{avgLen} = $avglen/$numseq;
		$out{maxLen} = $maxlen;
		$out{n50} = $n50;
		$out{minsize} = $minsize;
		$out{numTooSmall} = $toosmall;
	}
	else {
		$out{$numseq} = 0;
	}
	
	print "Leaving contigstats!\n" if $interested;
	return (\%out);
}


#toString method
sub toString {
    my $self = shift;
    my $tmp = $self->toStringNoV();
    if(defined $self->velvethout){
        $tmp .= "Velveth Output:\n" . $self->velvethout() . "\n";
    }
    if(defined $self->velvetgout){
        $tmp .= "Velvetg Output:\n" . $self->velvetgout() . "\n";
    }
    $tmp .= "**********************************************************\n";
    return $tmp;
}


#toStringNoV method
sub toStringNoV {
    my $self = shift;
    my $tmp = "********************************************************\n";
    if($self->ass_id()){
        $tmp .= "Assembly id: " . $self->ass_id(). "\n";
    }
    if($self->assmscore()){
        $tmp .= "Assembly score: " .$self->assmscore(). "\n";
    }
    if($self->timestamph()){
        $tmp .= "Velveth timestamp: " . $self->timestamph(). "\n";
    }
    if($self->timestampg()){
        $tmp .= "Velvetg timestamp: " . $self->timestampg(). "\n";
    }
    if(defined $self->versionh){
        $tmp .= "Velveth version: " . $self->versionh(). "\n";
    }
    if(defined $self->versiong){
        $tmp .= "Velvetg version: " . $self->versiong(). "\n";
    }
    if(defined $self->readfilename){
        $tmp .= "Readfile(s): " . $self->readfilename(). "\n";
    }
    if(defined $self->pstringh){
        $tmp .= "Velveth parameter string: " . $self->pstringh(). "\n";
    }
    if(defined $self->pstringg){
        $tmp .= "Velvetg parameter string: " . $self->pstringg(). "\n";
    }
    if(defined $self->ass_dir){
        $tmp .= "Assembly directory: " . $self->ass_dir(). "\n";
    }
    if(defined $self->hashval){
        $tmp .= "Velvet hash value: " . $self->hashval(). "\n";
    }
    if(defined $self->rmapfs){
        $tmp .= "Roadmap file size: " . $self->rmapfs(). "\n";
    }
    if(defined $self->sequences){
        $tmp .= "Total number of sequences: " . $self->sequences(). "\n";
    }
    if(defined $self->nconts){
        $tmp .= "Total number of contigs: " . $self->nconts(). "\n";
    }
    if(defined $self->n50){
        $tmp .= "n50: " . $self->n50(). "\n";
    }
    if(defined $self->maxlength){
        $tmp .= "length of longest contig: " . $self->maxlength(). "\n";
    }
    if(defined $self->totalbp){
        $tmp .= "Total bases in contigs: " . $self->totalbp(). "\n";
    }
    if(defined $self->nconts1k){
        $tmp .= "Number of contigs > 1k: " . $self->nconts1k(). "\n";
    }
    if(defined $self->totalbp1k){
        $tmp .= "Total bases in contigs > 1k: " . $self->totalbp1k(). "\n";
    }
    if($self->pstringh =~ /Pair/ && defined $self->pstringg && $self->pstringg =~ /-exp_cov/){
		$tmp .= "Paired Library insert stats:\n";
		my @x = split /\n/, $self->velvetgout;
		foreach(@x){
			chomp;
			if(/Paired-end library \d+ has/){
				s/^\[\d+\.\d+\]\s+//;
				$tmp .= "$_\n";
			}
		}
	}
    $tmp .= "**********************************************************\n";
    return $tmp;
}

sub opt_func_toString {
	my $out = "\nVelvet optimiser assembly optimisation function can be built from the following variables.\n";
	foreach my $key (sort keys %f_opts){
		$out .= "\t$key = " . $f_opts{$key}->{'desc'} . "\n";
	}
	$out .= "Examples are:\n\t'Lbp' = Just the total basepairs in contigs longer than 1kb\n";
	$out .= "\t'n50*Lcon' = The n50 times the number of long contigs.\n";
	$out .= "\t'n50*Lcon/tbp+log(Lbp)' = The n50 times the number of long contigs divided\n\t\tby the total bases in all contigs plus the log of the number of bases\n\t\tin long contigs.\n";
	return $out
}

1;
