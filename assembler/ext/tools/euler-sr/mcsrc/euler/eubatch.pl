#!/usr/bin/perl -w

###########################################################################
# Title:          euler.pl
# Author:         Glenn Tesler
# Description:    Run all components of EULER
# Created:        02/12/2002
# Last modified:  03/19/2002
# Last modified:  03/09/2004 by Haixu Tang for EULER 2.0
#
# Copyright (c) 2001-2004 The Regents of the University of California
# All Rights Reserved
# See file LICENSE for details.
###########################################################################

use Data::Dumper;

sub init_globals {
#    $pgms::bin = "/usr/local/apps/euler";
    $pgms::bin = getbindir();
}

sub print_usage {
    print STDERR <<'ENDusage';
euler.pl   -i fasta_file [optional args]
    -q quality_file     
    -v phrap|direct            default: direct
    -c1 cutoff_phrap_LLR_1     default: 15 (used w/ overlap=phrap)
    -c2 cutoff_phrap_LLR_2     default: 5 (used w/ overlap=phrap)
    -ol minimum_overlap_length default: 20 (used w/ overlap=direct)
    -id minimum_overlap_identity default: 0.95 (used w/ overlap=direct)
    -ws whirl_size             default: 50
    -co minimum_coverage       default: 2
    -realigner		       default: no realigner
    -nohtml  make text output  default: html output
    -script filename           make shell script w/cmds instead of doing them

euler.pl  -i fasta_file -restart [-script filename]
    restart an aborted run of euler (or create shell script to do that)

environment variable EULERBIN = directory with all euler binaries
ENDusage
}

sub getbindir {
    my $bindir = $ENV{'EULERBIN'};
    return $bindir                 if ($bindir);
    return undef                   if (!$0);

    # extract pathname from name of this script

    if ($0 =~ m:(.*)/:) {
	return $1;
    } else {
	return '.';
    }
}

sub main {
    my (@args) = @_;

    my $err = 0;

    init_globals();

    if (!defined $pgms::bin) {
	$err++;
	print STDERR "Environment variable EULERBIN undefined\n";
    }

    my $params = {
	# cmd line
	'fasta_file'      => undef,
	'quality_file'    => undef,
	'cutoff_llr_1'  => 15,
	'cutoff_llr_2'    => 5,
	'min_length'    => 20,
	'min_id'    => 0.95,
	'whirl_size'    => 50,
	'coverage'    => 2,
	'ovlp'            => "direct",
	'nohtml'            => 0,
	'cmethod'	=> "norealigner",
        'report'       => undef,
        'logfile'      => undef,

	# generating script file
	'makescript'      => 0,
	'script_file'     => undef,

	# other
	'bin'             => $pgms::bin,

	# restarting aborted euler
	'curphase'        => 'start',
	'execing'         => 1,        # may not execute all in restart
    };


    ######################################################################
    # non-cmd line parameters
    ######################################################################


    ######################################################################
    # parse command line parameters
    ######################################################################

    while ($#args >= 0 &&
	   ($args[0] =~ /\A-/)) {
	my $sw = shift @args;
	if ($sw eq '-i') {
	    $params->{'fasta_file'} = shift @args;
	} elsif ($sw eq '-q') {
	    $params->{'quality_file'} = shift @args;
	} elsif ($sw eq '-v') {
	    $params->{'ovlp'} = shift @args;
	} elsif ($sw eq '-c1') {
	    $params->{'cutoff_llr_1'} = shift @args;
	} elsif ($sw eq '-c2') {
	    $params->{'cutoff_llr_2'} = shift @args;
	} elsif ($sw eq '-ws') {
	    $params->{'whirl_size'} = shift @args;
	} elsif ($sw eq '-co') {
	    $params->{'coverage'} = shift @args;
	} elsif ($sw eq '-ol') {
	    $params->{'min_length'} = shift @args;
	} elsif ($sw eq '-id') {
	    $params->{'min_id'} = shift @args;
	} elsif ($sw eq '-script') {
	    $params->{'script_file'} = shift @args;
	    $params->{'makescript'} = 1;
	} elsif ($sw eq '-nohtml') {
	    $params->{'nohtml'} = 1;
	} elsif ($sw eq '-realigner') {
	    $params->{'cmethod'} = "realigner";
	} elsif ($sw eq '-restart') {
	    $params = init_restart($params);
	} else {
	    $err++;
	}

	last       if ($err);
    }

    # no extra parameters allowed
    $err++         if ($#args >= 0);

    ######################################################################
    # validate command line parameters
    ######################################################################

    $err ||= validate_params($params);

    if ($err) {
	print_usage();
	die "Invalid parameters";
    }


    ######################################################################
    # run euler
    ######################################################################

    run_euler($params);
}

sub validate_params {
    my ($params) = @_;

    #######################################################################

    if (!defined $params->{'fasta_file'})
    {
	errmsg("Must specify a fasta file with -i option.");
	return 1;
    }

    if (!-e $params->{'fasta_file'})
    {
	errmsg("Fasta file $params->{'fasta_file'} missing.");
	return 1;
    }

    #######################################################################

    if (defined $params->{'quality_file'}
	&& $params->{'quality_file'} eq 'no_file') {
	$params->{'quality_file'} = undef;
    }

    if (defined $params->{'quality_file'}
	&& !-e $params->{'quality_file'})
    {
	errmsg("Quality file $params->{'quality_file'} missing.");
	return 1;
    }

    #######################################################################

    # validate combination of trim and quality file

    my $ovlp = $params->{'ovlp'};

    if ($ovlp ne 'phrap' && $ovlp ne 'direct') {
        $ovlp = $params->{'ovlp'} = "direct";
    }

    return 0;
}

sub errmsg {
    print STDERR @_, "\n";
}


# load restart file
sub init_restart {
    my ($params) = @_;
    my $fasta_file = $params->{'fasta_file'};

    my $restart_file = $params->{'restart_file'} = "${fasta_file}.restart";

    open(RESTART, "<$restart_file") ||
	die "Can't open restart file.";

    my @old_params = <RESTART>;

    (close RESTART) || die "Can't close restart file.";

    # first line is
    #       $VAR1 = {
    $old_params[0] = "{";

    $params = eval join "", @old_params;

    # don't execute until we get to the correct phase
    $params->{'execing'} = 0;

    return $params;
}



# do_cmd($params,$phase,$cmd)
# do_cmd($params,$phase,$cmd,$no_date)
#
# allow euler to be restarted at whatever phase it left off at
# The current phase is a string written in file "${fasta_file}.restart"
# When it is completed, it should be overwritten with a new phase name.
#
# If we are executing, then have shell execute $cmd
# or write it to script file.
# Output of $cmd is redirected to a log file.
# If we are restarting, don't execute until $phase matches the
# phase where we left off.
#
# $no_date: true value means omit date stamp


sub do_cmd {
    my ($params, $phase, $cmd, $no_date) = @_;

    if (!$params->{'execing'}) {
	$params->{'execing'} =
	    $params->{'curphase'} eq $phase;
    }
    return              if (!$params->{'execing'});

    my $report_file = $params->{'report'};
    my $log_file = $params->{'logfile'};
    $cmd .= " >> $log_file";

    if ($params->{'makescript'}) {
	print SCRIPT_FILE "$cmd\n";
    } else {
	newphase($params, $phase);
	date_stamp($params)             if (!$no_date);
	`$cmd`;
#	if ($?) {
#	    `echo Error $?, aborted >> $report_file`;
#	    die;
#	}
    }
}

sub newphase {
    my ($params, $phase) = @_;

    $params->{'curphase'} = $phase;

    if (!$params->{'makescript'}) {
	open(PH_FILE, ">$params->{'restart_file'}") ||
	    die "Can't open restart file.";

	print PH_FILE Dumper($params);
	(close PH_FILE) || die "Can't close restart file.";
    }
}


sub date_stamp {
    my ($params) = @_;

    my $log_file = $params->{'logfile'};

    `echo "" >> $log_file`;
    `date >> $log_file`;
}

    
sub run_euler {
    my ($params) = @_;

    my $EBIN = $params->{'bin'};
#    $EBIN = "echo $EBIN";

    my $fasta_file      = $params->{'fasta_file'};
    my $quality_file    = $params->{'quality_file'};
    my $ovlp            = $params->{'ovlp'};
    my $cutoff_llr_1  = $params->{'cutoff_llr_1'};
    my $cutoff_llr_2  = $params->{'cutoff_llr_2'};
    my $min_length    = $params->{'min_length'};
    my $min_id        = $params->{'min_id'};
    my $whirl_size    = $params->{'whirl_size'};
    my $coverage      = $params->{'coverage'};
    my $cmethod      = $params->{'cmethod'};
    my $nohtml          = $params->{'nohtml'};

    my $EOUT;
    if($nohtml == 0)	{
	    $EOUT = "EULER-report.html";
    } else {
	    $EOUT = "EULER-report.txt";
    }
    my $ELOG = "$fasta_file.log";
    my $report_file     = $params->{'report'} = $EOUT;
    my $log_file     = $params->{'logfile'} = $ELOG;
    my $restart_file    = $params->{'restart_file'}
                        = "${fasta_file}.restart";

    #####################################################################
    # Print report header or scriptfile header
    #####################################################################

    my $describe_quality_file = $quality_file ? $quality_file : "(None)";
    my $resume_msg = "";

    # start script file or log file
    if ($params->{'makescript'}) {
	open(SCRIPT_FILE, ">$params->{'script_file'}") ||
	    die "Can't open script file.";
	print SCRIPT_FILE <<'SCRIPT_HEAD';
#!/bin/sh
SCRIPT_HEAD
    } else {
	# Start report file
	if ($params->{'execing'}) {
	    open(FH_EOUT, ">$EOUT") || die "Can't open report file.";
	} else {
	    open(FH_EOUT, ">>$EOUT") || die "Can't open report file.";
	    $resume_msg = <<ENDmsg;
Restarting at phase:  $params->{'curphase'}
ENDmsg
        }

        if($nohtml == 1)	{
             print FH_EOUT <<"ENDheader";
                            EULER report
EULER V2.0, Copyright (c) 2001-2004, The Regents of the University of California.
All Rights Reserved.
See details in the file LICENSE in the source code.
EULER V2.0 web portal users: see http://nbcr.sdsc.edu/euler

${resume_msg}Fasta file: $fasta_file
Quality file: $describe_quality_file
OverlapMethod: $ovlp<br>
ENDheader
	} else	{
             print FH_EOUT <<"ENDheader";
<html>
<head>
<title>EULER report</title>
<link REL=\"STYLESHEET\" HREF=\"euler.css\" TYPE=\"text/css\">
<meta http-equiv=\"Content-Language\" content=\"en-us\">
<meta http-equiv=\"Content-Type\" content=\"text/html; charset=ISO-8859-1\">
</head>
<body>
<H4>EULER report</H4>
<br>
EULER V2.0, Copyright (c) 2001-2004, The Regents of the University of California.<br>
All Rights Reserved.<br>
See details in the file LICENSE in the source code.<br>
EULER V2.0 web portal users: see http://nbcr.sdsc.edu/euler<br>
<br>
${resume_msg}Fasta file: $fasta_file<br>
Quality file: $describe_quality_file<br>
OverlapMethod: $ovlp<br>
<br>
ENDheader
	}

        if($nohtml == 0)	{
       	 if ($ovlp eq 'phrap') {
	    print FH_EOUT <<"ENDec";
Overlap cutoff 1 (LLR score for phrap): $cutoff_llr_1<br>
Overlap cutoff 2 (LLR score for phrap): $cutoff_llr_2<br>
ENDec
         } elsif ($ovlp eq 'direct') {
	    print FH_EOUT <<"ENDqc";
Overlap cutoff for minimum length: $min_length<br>
Overlap cutoff for minimum identity: $min_id<br>
ENDqc
         }
	 print FH_EOUT <<"ENDab";
Minimum whirl size: $whirl_size<br>
Minimum coverage: $coverage<br>
ENDab
        } else  {
       	 if ($ovlp eq 'phrap') {
	    print FH_EOUT <<"ENDop";
Overlap cutoff 1 (LLR score for phrap): $cutoff_llr_1
Overlap cutoff 2 (LLR score for phrap): $cutoff_llr_2
ENDop
         } elsif ($ovlp eq 'direct') {
	    print FH_EOUT <<"ENDdp";
Overlap cutoff for minimum length: $min_length
Overlap cutoff for minimum identity: $min_id
ENDdp
	 }
	 print FH_EOUT <<"ENDab";
Minimum whirl size: $whirl_size
Minimum coverage: $coverage
ENDab
	}
        close FH_EOUT;
    }

    my $cmdhtml;

    if($nohtml == 1)	{
	$cmdhtml = "";
    } else {
	$cmdhtml = "-H";
    }

    # producing alignment file (overlap detection)
    if ($ovlp eq 'phrap')	{
	do_cmd($params, "crossmatch.manyreads", "$EBIN/phrap/cross_match.manyreads ${fasta_file} >${fasta_file}.cm");
	do_cmd($params, "phrap-ovp.manyreads", "$EBIN/po-prod -s ${fasta_file} -i ${fasta_file}.cm -o phrap.ovp");
	do_cmd($params, "over-align-phrap",
		 "$EBIN/over-align-phrap -i $fasta_file -p phrap.ovp -o $fasta_file.aln -r $cutoff_llr_1 $cmdhtml");
    } else {
        if ( defined $quality_file ) {
	  if(${quality_file} ne "${fasta_file}.qual")	{
		`cp ${quality_file} ${fasta_file}.qual`;
	  }
	  do_cmd($params, "trans_qual", "$EBIN/trans_qual -i $fasta_file -o $fasta_file.reads $cmdhtml ");
	  do_cmd($params, "errcorr-lapper-pair-mem",
		 "$EBIN/errcorr-lapper-pair-mem.pl ${fasta_file}.reads ${fasta_file}.reads.score $min_length $min_id");
	} else {
	  do_cmd($params, "trans_qual", "$EBIN/trans_qual -i $fasta_file -o $fasta_file.reads $cmdhtml");
	  do_cmd($params, "errcorr-lapper-pair-mem",
		 "$EBIN/errcorr-lapper-pair-mem.pl ${fasta_file}.reads ${fasta_file}.reads.score $min_length $min_id");
	}
    }

    # A-Bruijn graph first run
    if ( $ovlp eq "phrap" ) {
       do_cmd($params,"over-repeat-new",
	   "$EBIN/over-repeat-new -i ${fasta_file}.aln -s ${fasta_file} -o ${fasta_file}.gvz -T -300 -e -200 -w $whirl_size -f 0.6 $cmdhtml");
    } else {
       do_cmd($params,"over-repeat-new",
	   "$EBIN/over-repeat-new -i ${fasta_file}.reads.new.aln -s ${fasta_file}.reads.new -o ${fasta_file}.reads.new.gvz -T -300 -e -200 -w $whirl_size -f 0.6 $cmdhtml");
    }

    # add overlaps between end reads
    if ( $ovlp eq "phrap" ) {
       do_cmd($params,"over-repeat-end",
	   "$EBIN/over-repeat-end-po -i ${fasta_file}.aln -s ${fasta_file} -o ${fasta_file}.aln.all -r $cutoff_llr_2 $cmdhtml");
       do_cmd($params,"errcorr_fin_mem",
	   "$EBIN/errcorr_fin_mem -i ${fasta_file}.aln.all -s ${fasta_file} -o ${fasta_file}.fin $cmdhtml");
    } else {
       do_cmd($params,"over-repeat-end",
	   "$EBIN/over-repeat-end -i ${fasta_file}.reads.new.aln -s ${fasta_file}.reads.new -q ${fasta_file}.reads.new.score -o ${fasta_file}.reads.new.aln.all -v 50 -d 0.90 -M 2 $cmdhtml");
       do_cmd($params,"errcorr_fin_mem",
	   "$EBIN/errcorr_fin_mem -i ${fasta_file}.reads.new.aln.all -s ${fasta_file}.reads.new -q ${fasta_file}.reads.new.score -o ${fasta_file}.fin $cmdhtml");
    }

    #A-Bruijn graph for the second run
    do_cmd($params,"over-repeat-new-fin",
	   "$EBIN/over-repeat-new -i ${fasta_file}.fin.aln -s ${fasta_file}.fin -o ${fasta_file}.fin.gvz -T -300 -e -200 -w $whirl_size -f 0.6 $cmdhtml");

    #cleanup temporary files
    do_cmd($params,"over-repeat-new/rm",
	   "rm G_*.gvz", 1);

    # EULER-ET
    do_cmd($params,"euler_et",
	   "$EBIN/euler_et -s ${fasta_file}.fin -o ${fasta_file}.fin.et.gvz -c $coverage -v -E 300 $cmdhtml");

    # EULER-DB
    do_cmd($params,"pairreads",
	   "$EBIN/pairreads -i ${fasta_file}.fin -o ${fasta_file}.fin.pair $cmdhtml");
    do_cmd($params,"euler_db",
	   "$EBIN/euler_db -i ${fasta_file}.fin -m ${fasta_file}.fin.pair -o ${fasta_file}.fin.db.gvz -v -Q 3 $cmdhtml");

    # EULER-Consensus
  if($cmethod eq "realigner")	{
    do_cmd($params,"consensus_et",
	   "$EBIN/consens.pl ${fasta_file}.fin.et.contig.ace ${fasta_file}.fin ${fasta_file}.fin.et.contig.con");
    do_cmd($params,"consensus_db",
	   "$EBIN/consens.pl ${fasta_file}.fin.db.contig.ace ${fasta_file}.fin ${fasta_file}.fin.db.contig.con");
  } else {
    do_cmd($params,"consensus_et",
	   "$EBIN/euler_cons -s ${fasta_file}.fin -g ${fasta_file}.fin.et.graph -i ${fasta_file}.fin.et.intv -e ${fasta_file}.fin.et.edge -c ${fasta_file}.fin.et.contig.con -a ${fasta_file}.fin.et.contig.con.ace");
    do_cmd($params,"consensus_db",
	   "$EBIN/euler_cons -s ${fasta_file}.fin -g ${fasta_file}.fin.db.graph -i ${fasta_file}.fin.db.intv -e ${fasta_file}.fin.db.edge -c ${fasta_file}.fin.db.contig.con -a ${fasta_file}.fin.db.contig.con.ace");
  }
    # OUTPUT multi-alignment of the reads
#    do_cmd($params, "makealn_et",
#       "makealn -s ${fasta_file}.fin -i ${fasta_file}.fin.et.contig.ace -a realigner.out -c ${fasta_file}.fin.et.contig -o ${fasta_file}.fin.et.contig.align");
#    do_cmd($params, "makealn_db",
#       "makealn -s ${fasta_file}.fin -i ${fasta_file}.fin.db.contig.ace -a realigner.out -c ${fasta_file}.fin.db.contig -o ${fasta_file}.fin.db.contig.align");
#
    # Clean temporary files
    # could use perl "unlink", but this way it will generate a script file properly
    if ( $ovlp eq "phrap" ) {
	do_cmd($params,"overlapper/rm",
	       "rm phrap.ovp",
	       1); # omit date stamp
    }

    #cleanup temporary files
  if($cmethod eq "realigner")	{
    do_cmd($params,"consensus/rm",
	   "rm realigner.inp realigner.out", 1);
  }

    # EULER-index
#{
#    my $EBIN = "/users/u4/gptesler/EU/euindex";
    if (defined $quality_file) {
        do_cmd($params,"euindex",
	       "$EBIN/euindex.pl -i ${fasta_file} -q ${quality_file} -o index.html");
    } else {
        do_cmd($params,"euindex",
	       "$EBIN/euindex.pl -i ${fasta_file} -o index.html");
    }
#}

    # print report tailor for html output
    if($nohtml == 0)	{
	open(FH_EOUT, ">>$EOUT") || die "Can't open report file.";
        print FH_EOUT <<"ENDtailor";
</body>
</html>
ENDtailor
	close(FH_EOUT);
    }

    newphase($params,"done");

    if ($params->{'makescript'}) {
	close SCRIPT_FILE || die "Can't close script file";
    } else {
	unlink $restart_file;
    }


}

main(@ARGV);

