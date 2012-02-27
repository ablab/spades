#!/usr/bin/perl -w

###########################################################################
# Title:          euindex.pl
# Author:         Glenn Tesler
# Description:    Generate an index of all the EULER files in a directory
# Created:        02/05/2002
# Last modified:  03/19/2002
# Modified:	  02/16/2004 by Haixu Tang to incorporate the changes
#		  in EULER 2.0
#
# Copyright (c) 2001-2004 The Regents of the University of California
# All Rights Reserved
# See file LICENSE for details.
###########################################################################

#use Data::Dumper;

sub init_globals {
#    $pgms::bin = "/usr/local/apps/euler";
    $pgms::bin = getbindir();

    my $bin = $pgms::bin ? $pgms::bin : ".";
    $pgms::FASTASTAT = "${bin}/seqcount.pl";

    $colors::pagebg = "lemonchiffon";
    $colors::pagetext = "#000000";
    $colors::link = "#0000EE";
    $colors::vlink = "#551A8B";
    $colors::alink = "#FF0000";
    $colors::tablebg = "khaki";
    $colors::tablebg_rh = "bisque";
    $colors::file_missing = "red";
}

sub print_usage {
    print <<'ENDusage';
euindex.pl   -i reads_input_name  [optional args]
   -i reads_input_name: required
   -q quality_input_name: optional, default none
   -d dirname: directory with input files. default = current directory
   -r repname: default = EULER-report.html

   -o output_name: default = index.html

   -m: print names of all other files in the directory missing from index

environment variable EULERBIN = directory with all euler binaries
ENDusage
#  -r repname: default = $reads_input_name.out
#  -or output_report_name: default = reportname.html
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


    my $params = {
	'dir'      => '.',
	'repname'  => undef,
	'xxxx'     => undef,
	'missing'  => 0,
	'qualname' => undef,
	'indexfname' => undef,
#	'repnamehtml' => undef,
    };


    ######################################################################
    # non-cmd line parameters
    ######################################################################

    init_globals();

    if (!defined $pgms::bin) {
	$err++;
	print STDERR "Environment variable EULERBIN undefined\n";
    }

    ######################################################################
    # parse command line parameters
    ######################################################################

    # parse switches
    while ($#args >= 0  &&
	   ($args[0] =~ /\A-/)) {
	my $sw = shift @args;
	if ($sw eq '-d') {
	    $params->{'dir'} = shift @args;
	} elsif ($sw eq '-i') {
	    $params->{'xxxx'} = shift @args;
	} elsif ($sw eq '-r') {
	    $params->{'repname'} = shift @args;
	} elsif ($sw eq '-q') {
	    $params->{'qualname'} = shift @args;
	} elsif ($sw eq '-o') {
	    $params->{'indexfname'} = shift @args;
#	} elsif ($sw eq '-or') {
#	    $params->{'repnamehtml'} = shift @args;
	} elsif ($sw eq '-m') {
	    $params->{'missing'} = 1;
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
	die;
    }

    ######################################################################
    # make index
    ######################################################################

    get_file_names($params);

#    htmlify_report($params);

    print_index($params);
#    print Dumper($params);
}


sub validate_params {
    my ($params) = @_;

    #######################################################################

    if (!defined $params->{'xxxx'})
    {
	errmsg("Must specify a fasta file with -i option.");
	return 1;
    }

    my $xxxx = $params->{'xxxx'};

    #######################################################################

    if (!defined $params->{'repname'})
    {
#	$params->{'repname'} = "$xxxx.out";
	$params->{'repname'} = "EULER-report.html";
    }

    #######################################################################

    if (!defined $params->{'indexfname'}) {
	$params->{'indexfname'} = "index.html";
    }

    #######################################################################

    return 0;
}

sub errmsg {
    print STDERR @_, "\n";
}


###############################################################################
###############################################################################




# Create list of all filenames to index
# Most follow certain naming patterns based on the input file names
# Some are found by examining the directory
sub get_file_names {
    my ($params) = @_;

    my $xxxx = $params->{'xxxx'};
    my $dir = $params->{'dir'};
    my $repname = $params->{'repname'};
    my $qualname = $params->{'qualname'};
#    my $indexfname = $params->{'indexfname'};

    my $fnames = $params->{'fnames'} = {};
    my $desc = $params->{'desc'} = {};

    # read in all file names in directory
    getdir($params);


    ####################################################################
    # inputs
    ####################################################################

    $desc->{'reads_input'} = "Reads-Input";
    $fnames->{'reads_input'} = $xxxx;

    $desc->{'quality_input'} = "Quality-Input";
    $fnames->{'quality_input'} = $qualname;

    $desc->{'name.rul'} = "Mate pairs rules";
    $fnames->{'name.rul'} = "name.rul";

#    $desc->{'step.inp'} = "Error correction rules";
#    $fnames->{'step.inp'} = "step.inp";

    ####################################################################
    # other global
    ####################################################################

    $desc->{'report.txt'} = "Report";
    $fnames->{'report.txt'} = $repname;

    $desc->{'restart_file'} = "Restart file";
    $fnames->{'restart_file'} = "$xxxx.restart";

    ####################################################################
    # outputs: classify reads
    ####################################################################

    $desc->{'reads_used'} = "Reads-Used";
    $fnames->{'reads_used'} = "$xxxx.reads";

    $desc->{'quality_used'} = "Quality-Used";
    $fnames->{'quality_used'} = "$xxxx.reads.qual";

#    $desc->{'reads_fin'} = "Reads-Finishing";
#    $fnames->{'reads_fin'} = "$xxxx.fin";

     $desc->{'reads_fin'} = "Trimmed Reads";
     $fnames->{'reads_fin'} = "$xxxx.fin";

#    $desc->{'quality_fin'} = "Quality-Finishing";
#    $fnames->{'quality_fin'} = "$xxxx.fin.qual";

     $desc->{'quality_fin'} = "Trimmed Quality";
     $fnames->{'quality_fin'} = "$xxxx.fin.qual";

#    $desc->{'reads_long'} = "Reads-Long";
#    $fnames->{'reads_long'} = "$xxxx.long";

#    $desc->{'quality_long'} = "Quality-Long";
#    $fnames->{'quality_long'} = "$xxxx.long.qual";

#    $desc->{'reads_ambig'} = "Reads-Ambiguous";
#    $fnames->{'reads_ambig'} = "$xxxx.ambig";

#    $desc->{'quality_ambig'} = "Quality-Ambiguous";
#    $fnames->{'quality_ambig'} = "$xxxx.ambig.qual";

    ####################################################################
    # outputs: EULER-Trim
    ####################################################################

#    $desc->{'reads_trim'} = "Reads-Trim";
#    $fnames->{'reads_trim'} = "$xxxx.clean";

    ####################################################################
    # outputs: EULER-EC
    ####################################################################

#    $desc->{'reads_ec'} = "Reads-EC";
#    $fnames->{'reads_ec'} = "$xxxx.fix";

    ####################################################################
    # outputs: EULER-ChimDet
    ####################################################################

#    $desc->{'reads_chimgood'} = "Reads-NotChim";
#    $fnames->{'reads_chimgood'} = "$xxxx.fil";

#    $desc->{'reads_chimbad'} = "Reads-Chim";
#    $fnames->{'reads_chimbad'} = "$xxxx.fix.chim";

#    $desc->{'reads_chimunrel'} = "Reads-Unrel";
#    $fnames->{'reads_chimunrel'} = "$xxxx.fix.unrel";

    ####################################################################
    # outputs: FragmentGluer
    ####################################################################

    $desc->{'edge_fg'} = "Edge sequences-FG";
    $fnames->{'edge_fg'} = "$xxxx.fin.edge";

    $desc->{'graph_fg'} = "Graph architecture-FG";
    $fnames->{'graph_fg'} = "$xxxx.fin.graph";

    $desc->{'intv_fg'} = "Fragment intevals-FG";
    $fnames->{'intv_fg'} = "$xxxx.fin.intv";

    $desc->{'comp_fg_gviz'} = "Repeat-Graph-FragmentGluer (GVIZ)";
    $fnames->{'comp_fg_gviz'} = "$xxxx.fin.gvz";

    $desc->{'contigs_fg'} = "Contigs-FG";
    $fnames->{'contigs_fg'} = "$xxxx.fin.contig";

    $desc->{'ace_fg'} = "Contigs-ACE-FG";
    $fnames->{'ace_fg'} = "$xxxx.fin.contig.ace";

    ####################################################################
    # outputs: EULER-ET
    ####################################################################

#    $desc->{'edge_et'} = "Edge sequences-ET";
#    $fnames->{'edge_et'} = "$xxxx.fil.edge";
#
#    $desc->{'graph_et'} = "Graph architecture-ET";
#    $fnames->{'graph_et'} = "$xxxx.fil.graph";
#
#    $desc->{'path_et'} = "Superpaths-ET";
#    $fnames->{'path_et'} = "$xxxx.fil.path";
#
#    $desc->{'comp_et_gviz'} = "Repeat-Graph-ET (GVIZ)";
#    $fnames->{'comp_et_gviz'} =
#	list_num_files($params,"$xxxx.fil_et_comp#.gvz");
#
#
#    $desc->{'contigs_et'} = "Contigs-ET";
#    $fnames->{'contigs_et'} = "$xxxx.fil.contig";
#
#    $desc->{'single_et'} = "Single/double read contigs";
#    $fnames->{'single_et'} = "$xxxx.fil.singleton";
#

    $desc->{'edge_et'} = "Edge sequences-ET";
    $fnames->{'edge_et'} = "$xxxx.fin.et.edge";

    $desc->{'graph_et'} = "Graph architecture-ET";
    $fnames->{'graph_et'} = "$xxxx.fin.et.graph";

    $desc->{'intv_et'} = "Fragment intevals-ET";
    $fnames->{'intv_et'} = "$xxxx.fin.et.intv";

    $desc->{'comp_et_gviz'} = "Repeat-Graph-ET (GVIZ)";
    $fnames->{'comp_et_gviz'} = "$xxxx.fin.et.gvz";

    $desc->{'contigs_et'} = "Contigs-ET";
    $fnames->{'contigs_et'} = "$xxxx.fin.et.contig";

    $desc->{'ace_et'} = "Contigs-ACE-ET";
    $fnames->{'ace_et'} = "$xxxx.fin.et.contig.ace";

    ####################################################################
    # outputs: EULER-DB-preprocess
    ####################################################################

    $desc->{'mate_all'} = "Mate pairs-all";
##    $fnames->{'mate_all'} = "$xxxx.reads.clean.pair";
#    $fnames->{'mate_all'} = "$xxxx.fil.mate";
     $fnames->{'mate_all'} = "$xxxx.fin.pair";


    ####################################################################
    # outputs: EULER-DB
    ####################################################################

#    $desc->{'edge_db'} = "Edge sequences-DB";
#    $fnames->{'edge_db'} = "$xxxx.fil.mate.edge";
#
#    $desc->{'graph_db'} = "Graph architecture-DB";
#    $fnames->{'graph_db'} = "$xxxx.fil.mate.graph";
#
##    $desc->{'path_db'} = "Superpaths-DB";
##    $fnames->{'path_db'} = "$xxxx.fil.mate.path";
#
#    # TODO: Components-DB
#    $desc->{'comp_db_gviz'} = "Components-DB (GVIZ)";
#    $fnames->{'comp_db_gviz'} =
#	list_num_files($params,"$xxxx.fil_db_comp#.gvz");
#
#
#
#    $desc->{'contigs_db'} = "Contigs-DB";
#    $fnames->{'contigs_db'} = "$xxxx.fil.mate.contig";
#
#    $desc->{'mate_db'} = "Mate pairs-DB";
##    $fnames->{'mate_db'} = "$xxxx.mate.pair";
#    $fnames->{'mate_db'} = "$xxxx.fil.et.mate";
#

    $desc->{'edge_db'} = "Edge sequences-DB";
    $fnames->{'edge_db'} = "$xxxx.fin.db.edge";

    $desc->{'graph_db'} = "Graph architecture-DB";
    $fnames->{'graph_db'} = "$xxxx.fin.db.graph";

    $desc->{'intv_db'} = "Fragment intevals-DB";
    $fnames->{'intv_db'} = "$xxxx.fin.db.intv";

    $desc->{'comp_db_gviz'} = "Repeat-Graph-DB (GVIZ)";
    $fnames->{'comp_db_gviz'} = "$xxxx.fin.db.gvz";

    $desc->{'contigs_db'} = "Contigs-DB";
    $fnames->{'contigs_db'} = "$xxxx.fin.db.contig";

    $desc->{'ace_db'} = "Contigs-ACE-DB";
    $fnames->{'ace_db'} = "$xxxx.fin.db.contig.ace";

    $desc->{'comp_sf_gviz'} = "Scaffolding-Graph (GVIZ)";
    $fnames->{'comp_sf_gviz'} = "$xxxx.fin.db.contig.sf.gvz";

    ####################################################################
    # outputs: EULER-Consensus
    ####################################################################

##    $desc->{'contigs_cons'} = "Contigs-Consensus";
##    $fnames->{'contigs_cons'} = "$xxxx.fil.mate.con";
    $desc->{'contigs_et_cons'} = "Contigs-Consensus-ET";
    $fnames->{'contigs_et_cons'} = "$xxxx.fin.et.contig.con";
    $desc->{'contigs_et_align'} = "Contigs-Alignment-ET";
    $fnames->{'contigs_et_align'} = "$xxxx.fin.et.contig.align";

    $desc->{'contigs_db_cons'} = "Contigs-Consensus-DB";
    $fnames->{'contigs_db_cons'} = "$xxxx.fin.db.contig.con";
    $desc->{'contigs_db_align'} = "Contigs-Alignment-DB";
    $fnames->{'contigs_db_align'} = "$xxxx.fin.db.contig.align";

    ####################################################################
    # outputs: EULER-SF (1)   (1st run of EULER-SF)
    ####################################################################

#    $desc->{'scaf1_gviz'} = "Scaffolding 1 (GVIZ)";
#    $fnames->{'scaf1_gviz'} =
##	list_num_files($params,
##		       "$xxxx.fil.mate.contig_sf.gvz.#");
#	list_num_files($params,
#		       "$xxxx.fil.mate.con_sf_comp#.gvz");
#
#
#    ####################################################################
#    # outputs: EULER-Connect
#    ####################################################################
#
#    $desc->{'contigs_connt'} = "Contigs-Connect";
#    $fnames->{'contigs_connt'} = "$xxxx.fil.mate.con.connt";
#
#    $desc->{'chimeric_connt'} = "Contigs-Connect-Chim";
#    $fnames->{'chimeric_connt'} = "$xxxx.fil.mate.con.connt.chim";
#
#    $desc->{'contigs_connt_gviz'} = "Contigs-Connect (GVIZ)";
#    $fnames->{'contigs_connt_gviz'} = "$xxxx.fil.mate.con_connt.gvz";
#
#    ####################################################################
#    # outputs: EULER-SF (2)
#    ####################################################################
#
#    $desc->{'scaf2_gviz'} = "Scaffolding 2 (GVIZ)";
#    $fnames->{'scaf2_gviz'} =
##	list_num_files($params,
#		       "$xxxx.fil.mate.contig.connt_sf.gvz.#");
#
#		       "$xxxx.fil.mate.con.connt_sf_comp#.gvz");
#
##    print Dumper($fnames);
}


# choose_file($params,$name1,$name2,...)
# chooses the first filename that actually exists.
# If none exists, returns $name1.
sub choose_file {
    my ($params,@fnames) = @_;

    my $dir_list_used = $params->{'dir_list_used'};

    foreach my $name (@fnames) {
	return $name         if (defined $dir_list_used->{$name});
    }

    return $fnames[0];
}


# get a list of numbered files
#   list_num_files($params, $pattern)
#      $pattern = 'xxxx#yyyy'
#      All file names in the directory $params->$dir
#      that have a number where the "#" symbol is, will
#      be put into a list
#         [ i1 => name1, i2 => name2, ... ]
#      in numerical order.

sub list_num_files {
    my ($params,$pattern) = @_;

    # read in the directory
    my $allfiles = getdir($params);

    # turn pattern into   \Axxxx(\d+)yyyy\Z
    my $pat2 = $pattern;
    $pat2 =~ s/#/\(\\d+\)/ ;
    $pat2 = "\\A${pat2}\\Z";

    # find the filenames matching the pattern
    my %numfiles = map   { /$pat2/ ? ($1=>$_) : () }    @$allfiles;


    # sort the files in numerical order
    my @keys = sort numeric keys %numfiles;
    my @snumfiles = map { ($_, $numfiles{$_}) } @keys;
    # TODO

    return \@snumfiles;
}

sub numeric { $a <=> $b };


# first call:
#      1. read all filenames in current directory into an array
#      2. make table $used{$filename}=0
#
# subsequent calls: just recall the array of filenames
#
# return reference to the array


sub getdir {
    my ($params) = @_;

    my $dir_list = $params->{'dir_list'};
    return $dir_list      if (defined $dir_list);

    my $dir = $params->{'dir'};
    opendir(CURDIR, $dir)   or   die "Could not open directory $dir";


    my @allfiles = readdir CURDIR;
    closedir CURDIR;

    $dir_list = $params->{'dir_list'} = \@allfiles;

    $params->{'dir_list_used'} = {
	map { ($_=>0) } @allfiles
    };

    return $dir_list;
}

###############################################################################

# generate the brief index
sub brief_table {
    my ($params) = @_;

    my $result = "";

    $result .= "<table border=1 bgcolor=$colors::tablebg>";

    $result .=
	briefrow($params,
		 "Inputs",
		 [['reads_input',optqual($params,'quality_input')],
#		  ['+','stat_fas','reads_input'],
		  ]);

    $result .=
#	briefrow($params,
#		 "Control files",
#		 [['name.rul'],
#		  ['step.inp']
#		  ]);
	briefrow($params,
		 "Control files",
		 [['name.rul']]);

     $result .=
 	briefrow($params,
 		 "Settings",
 		 [['stat_iparams','report.txt']
 		  ]
 		 );

    $result .=
	briefrow($params,
		 "Assemblies",
		 [['@', 'EULER-FragmentGluer:', 'contigs_fg'],
#		  ['+','stat_fas','contigs_fg'],
		  ['-'],
		  ['@', 'EULER-ET:'],
		  ['+', 'contigs_et_cons'],
		  ['+', 'contigs_et_align'],
#		  ['+','stat_fas','contigs_et'],
		  ['-'],
		  ['@','EULER-DB:'],
		  ['+','contigs_db_cons'],
		  ['+','contigs_db_align'],
#		  ['+','+','stat_fas','contigs_db'],
#		  ['+','@','Scaffolding 1:'],
#		  ['+','+','scaf1_gviz'],
#		  ['-'],
#		  ['@','EULER-Connect:'],
#		  ['+','contigs_connt'],
#		  ['+','+','stat_fas','contigs_connt'],
#		  ['+','@','Scaffolding 2:'],
#		  ['+','+','scaf2_gviz']
		  ]
		 );

    $result .=
	briefrow($params,
		 "Report",
		 [['report.txt']]);


    $result .= "</table>";


    return $result;
}


# briefrow($params,
#          "Row name",
#          [[line1], [line2], ...])

sub briefrow {
    my ($params, $rowname, $lines) = @_;

    my $col1 = "<th bgcolor=$colors::tablebg_rh>$rowname</th>";

    my $col2 = fmat_box($params,$lines);
    
    return "<tr>${col1}${col2}</tr>";
}


# full_row($params,
#          "Row name",
#          [[line1], [line2], ...],    # for INPUT
#          [[line1], [line2], ...]     # for OUTPUT
#         )
# same codes as in briefrow

sub full_row {
    my ($params, $rowname, $inputs, $outputs) = @_;

    my $col1 = "<th bgcolor=$colors::tablebg_rh>$rowname</th>";
    my $col2 = fmat_box($params,$inputs);
    my $col3 = fmat_box($params,$outputs);
    
    return "<tr>${col1}${col2}${col3}</tr>";
}




# generate the full index
sub detailed_table {
    my ($params) = @_;

    my $result = "";


    my $thcolor = "sandybrown";
    $result .= <<"TABLE_HEAD";
<table border=1 bgcolor=$colors::tablebg>
<thead>
<th bgcolor=$thcolor>Phase</th>
<th bgcolor=$thcolor>Inputs</th>
<th bgcolor=$thcolor>Outputs</th>
</thead>
<tbody>
TABLE_HEAD

    $result .= full_row($params,
			"Classify reads",
			[['reads_input', optqual($params,'quality_input')],
#			 ['+','stat_fas','reads_input'],
			 ['name.rul']],

#			[['reads_used', optqual($params,'quality_used')],
#			 ['+','stat_fas','reads_used'],

  			 [['reads_fin'],
# 			 [['reads_fin', optqual($params,'quality_fin')],
#			 ['+','stat_fas','reads_fin'],
#
#			 ['reads_long', optqual($params,'quality_long')],
#			 ['+','stat_fas','reads_long'],
#
#			 ['reads_ambig', optqual($params,'quality_ambig')],
#			 ['+','stat_fas','reads_ambig']],
			]
			);

#    $result .= full_row($params,
#			"EULER-Trim",
#			[['reads_used',optqual($params,'quality_used')]],
#			[['reads_trim'],
#			 ['+','stat_fas','reads_trim'],
#			 ]
#			);
#
#    $result .= full_row($params,
#			"EULER-EC",
#			[['reads_trim'],
#			 ['step.inp']],
#			[['reads_ec'],
#			 ['+','stat_fas','reads_ec'],
#			 ]
#			);
#
#    $result .= full_row($params,
#			"EULER-ChimDet",
#			[['reads_ec']],
#			[['reads_chimgood'],
#			 ['+','stat_fas','reads_chimgood'],
#
#			 ['reads_chimbad'],
#			 ['+','stat_fas','reads_chimbad'],
#
#			 ['reads_chimunrel'],
#			 ['+','stat_fas','reads_chimunrel'],
#
#			 ]
#			);
#
     $result .= full_row($params,
 			"EULER-FragmentGluer",
 			[['reads_fin']],
#			[['reads_fin'],['quality_fin']],
			[['@','EULER-FG:'],
			 ['+','edge_fg'],
#			 ['+','+','stat_fas','edge_fg'],
			 ['+','graph_fg'],
			 ['+','+','stat_graph','graph_fg'],
			 ['+','intv_fg'],
			 ['+','@','Repeat-Graph-FragmentGluer (GVIZ):'],
			 ['+','+','comp_fg_gviz'],
			 ['+','contigs_fg'],
#			 ['+','stat_fas','contigs_fg'],
			 ['+','ace_fg'],
#			 ['+','ace_fg'],
			 ['-']],
			);

     $result .= full_row($params,
 			"EULER-ET",
#			[['reads_fin'],['quality_fin']],
			[['reads_fin']],
			[['@','EULER-ET:'],
			 ['+','edge_et'],
#			 ['+','+','stat_fas','edge_et'],
			 ['+','graph_et'],
			 ['+','+','stat_graph','graph_et'],
#			 ['+','path_et'],
#			 ['+','+','stat_path','path_et'],
			 ['+','intv_et'],
			 ['+','@','Repeat-Graph-ET (GVIZ):'],
			 ['+','+','comp_et_gviz'],
			 ['+','contigs_et'],
#			 ['+','stat_fas','contigs_et'],
			 ['+','ace_et'],
			 ['-']],
			);

     $result .= full_row($params,
			 "EULER-DB",
			[['reads_fin'], ['name.rul']],
			[['@','EULER-DB:'],
			 ['+','mate_all'],
#			 ['+','stat_mate','mate_all'],
			 ['+','edge_db'],
#			 ['+','+','stat_fas','edge_db'],
			 ['+','graph_db'],
			 ['+','+','stat_graph','graph_db'],
			 ['+','intv_db'],
			 ['+','@','Repeat-Graph-DB (GVIZ):'],
			 ['+','+','comp_db_gviz'],
			 ['+','contigs_db'],
#			 ['+','stat_fas','contigs_db'],
			 ['+','ace_db'],
			 ['-']]
			);

#    $result .= full_row($params,
#			"EULER-DB (phase 1)",
#			[['reads_fin'],
#			 ['name.rul']],
#			[['mate_all'],
#			 ['+','stat_mate','mate_all']]
#			);
#
#    $result .= full_row($params,
#			"EULER-DB (phase 2)",
#			[['reads_chimgood'],
#			 ['mate_all'],
#			 ['-'],
#			[['@','Graph-DB:'],
#			 ['+','edge_db'],
#			 ['+','+','stat_fas','edge_db'],
#			 ['+','graph_db'],
#			 ['+','+','stat_graph','graph_db'],
##			 ['+','path_db'],
##			 ['+','+','stat_path','path_db'],
#			 ['+','@','Components-DB (GVIZ):'],
#			 ['+','+','comp_db_gviz'],
#			 ['-'],
#			 ['contigs_db'],
#			 ['+','stat_fas','contigs_db'],
#			 ['mate_db'],
#			 ['+','stat_mate','mate_db']]
#			);
#
    $result .= full_row($params,
			"EULER-Consensus",
			[['reads_fin'],
#			 ['contigs_db'],
#			 ['-'],
#			 ['@','Graph-DB:'],
#			 ['+','edge_db'],
#			 ['+','graph_db'],
#			 ['+','path_db']
			 ],
			[['@','Scaffolding-Graph'],
			 ['+','contigs_et_cons'],
			 ['+','contigs_et_align'],
#			 ['+','stat_fas','contigs_et_cons'],
			 ['+','contigs_db_cons'],
			 ['+','contigs_db_align'],
#			 ['+','stat_fas','contigs_db_cons']
			 ]
			);

    $result .= full_row($params,
 			"EULER-SF",
#			"EULER-SF (pass 1)",
			[['reads_fin'],
#			 ['comp_sf_gviz'],
			 ['mate_all']],
			[['@','Scaffolding-Graph'],
			 ['+', 'comp_sf_gviz']]
			);
#
#    $result .= full_row($params,
#			"EULER-Connect",
#			[['reads_used'],
#			 ['contigs_cons']],
#			[['contigs_connt'],
#			 ['+','stat_fas','contigs_connt'],
#			 ['chimeric_connt'],
#			 ['contigs_connt_gviz']]
#			);
#
#    $result .= full_row($params,
#			"EULER-SF (pass 2)",
#			[['reads_chimgood'],
#			 ['contigs_connt'],
#			 ['mate_all']],
#			[['@','Scaffolding 2:'],
#			 ['+', 'scaf2_gviz']]
#			);
#
#

    $result .= <<'TABLE_TAIL';
</tbody>
</table>
TABLE_TAIL


    return $result;
}


# only use quality files if they are present
sub optqual {
    my ($params, $qual_key) = @_;
    my $qualname = $params->{'qualname'};

    return $qualname ? $qual_key : ();
}

sub missing_files {
    my ($params) = @_;

    my $dir_list = $params->{'dir_list'};
    my $dir_used = $params->{'dir_list_used'};
    my $fnames = $params->{'fnames'};


    # omit certain files
    $dir_used->{'.'} = $dir_used->{'..'} = 1;
    $dir_used->{$fnames->{'restart_file'}} = 1;

    my $indexfname = $params->{'indexfname'};
    $indexfname =~ s:.*/::g ;    # get filename portion w/o path
    $dir_used->{$indexfname} = 1;

    

    my @otherfiles =
	grep { !$dir_used->{$_} } @$dir_list;

    return ""         if (!@otherfiles);

    my $result = "<PRE>\n";

    @otherfiles = sort @otherfiles;

    foreach my $fname (@otherfiles) {
	$result .= $fname;
	$result .= "\n";
    }

    $result .= "</PRE>\n";

    return $result;
}


# fmat_box($params, [[line1], [line2], ...])
# each line is:
#     '-':          separator
#     '@', 'str':   title 'str'
#     '+', ...:     indent ..., where ... is one of the codes in this list
#     id, ...:      generate hyperlink for code w/given ID, then do ...
#     'stat_fas', id, ...: generate size stats for FASTA file w/code ID
#     'stat_graph', id, ...:
#                    stats for graph architecture
#     'stat_path', id, ...:
#                    stats for superpath file
#     'stat_mate', id, ...:
#                    stats for matepair file
# TODO:
#     'stat_misc', id, ...:
#                    special statistic parsed from report file


sub fmat_box {
    my ($params, $lines) = @_;

    my $result = "<td valign=top>";

    my $firstline = 1;
    my $linepref = "";
    
    foreach my $rawline (@$lines) {
	if ($firstline) {
	    $firstline = !$firstline;
	    $linepref = "";
	} else {
	    $linepref = "<br>";
	}


	my @line = @$rawline;

	my $firstentry = 1;

	while ($#line >= 0) {
	    my $op = shift @line;

	    if ($op eq '-') {
		$result .= '<hr>';
		$firstline = 1;
		$linepref = "";
	    } elsif ($op eq '+') {
		$result .= $linepref; $linepref = "";
		$result .= "&nbsp;&nbsp;&nbsp;";
	    } elsif ($op eq '@') {
		$result .= $linepref; $linepref = "";
		$result .= "<b>" . (shift @line) . "</b>";
	    } elsif ($op =~ /\Astat_/) {
		my $id = shift @line;

		$result .= $linepref; $linepref = "";
		if ($firstentry) {
		    $firstentry = !$firstentry;
		} else {
		    $result .= "&nbsp;";
		}


		# if already have this kind of comment for file, recall it
		# instead of recomputing it

		my $comment = $params->{'remember'}->{$op}->{$id};

		if (!defined $comment) {
		    if ($op eq 'stat_fas') {
			$comment = stat_fasta($params,$id);
		    } elsif ($op eq 'stat_graph') {
			$comment = stat_graph($params,$id);
		    } elsif ($op eq 'stat_path') {
			$comment = stat_path($params,$id);
		    } elsif ($op eq 'stat_mate') {
			$comment = stat_mate($params,$id);
		    } elsif ($op eq 'stat_iparams') {
			$comment = stat_iparams($params,$id);
		    }

		    $params->{'remember'}->{$op}->{$id} = $comment;
		}
#		else { $comment .= "REMEMBERED"; }

		$result .= $comment;

	    } else {
		$result .= $linepref; $linepref = "";
		if ($firstentry) {
		    $firstentry = !$firstentry;
		} else {
		    $result .= ", ";
		}

		$result .= hlink($params,$op);
	    }
	}
    }

    $result .= "</td>";

    return $result;
}

# hlink($params,$key)
# generate hyperlinks to file(s) given by $key
sub hlink {
    my ($params, $key) = @_;

    my $fnames = $params->{'fnames'};
    my $desc = $params->{'desc'};

    # it's a single filename
    if (!ref($fnames->{$key})) {
	return hlinkfile($params, $desc->{$key}, $fnames->{$key});
    }

    # it's an array of filenames
    my $result = "";
    my $flist_ref = $fnames->{$key};

    my $numnames = scalar(@$flist_ref);

#    print "filename array: $key,  arraysize=$numnames\n";
#    print Dumper($fnames->{$key});


    for (my $i = 0; $i < $numnames; $i += 2) {
	$result .= ", "       if ($i > 0);

	$result .= hlinkfile($params, $flist_ref->[$i], $flist_ref->[$i+1]);
    }

    return $result;
}


# hlinkfile($params, $desc, $filename)
# generate hyperlink of description $desc to file $filename

sub hlinkfile {
    my ($params, $desc, $filename) = @_;

    my $dir = $params->{'dir'};
    my $dir_used = $params->{'dir_list_used'};

    if (!defined $dir_used->{$filename}) {
#	print "desc: $desc\n";
#	print "filename: $filename\n";
	$desc = "<font color=\"$colors::file_missing\">$desc</font>";
    } else {
	$dir_used->{$filename} = 1;
    }

    return "<a href=\"${filename}\">${desc}</a>";
}


###############################################################################
# file statistics

# stat_fasta($params,$key)
# return string with statistics for FASTA file w/id $key
sub stat_fasta {
    my ($params, $key) = @_;

    my $fnames = $params->{'fnames'};
    my $dir = $params->{'dir'};
    my $fname = "$dir/$fnames->{$key}";

    my ($stats,$err) = mysystem("$pgms::FASTASTAT $fname");

    # should be line like
    # 100 seqs, 377783 bp, 29-73320 bp/seq, 3777.83 avg

    if ($err) {
	return "N/A";
    } elsif ($stats =~ m:(\d+) seqs, (\d+) bp, (\d+)-(\d+) bp/seq, ([\d\.]+) avg:) {
	return ($1 == 1) ?
	    "($1 seq, $2 bp, $3-$4 bp/seq, $5 avg)" :
		"($1 seqs, $2 bp, $3-$4 bp/seq, $5 avg)" ;

    } elsif ($stats =~ /\AEmpty/) {
	return "(0 seqs)";
    } else {
	return "N/A";
    }
}


# stat_graph($params,$key)
# return string with statistics for graph architecture file w/id $key
sub stat_graph {
    my ($params, $key) = @_;

    my $fnames = $params->{'fnames'};
    my $dir = $params->{'dir'};
    my $fname = "$dir/$fnames->{$key}";

    my ($stats,$err) = mysystem("head -1 $fname");

    # should be line like
    # Number_of_Vertex 193

    if ($err) {
	return "N/A";
    } elsif ($stats =~ m:Number_of_Vertex (\d+):) {
	return ($1 == 1) ?
	    "($1 vertex)" :
		"($1 vertices)" ;
    } else {
	return "N/A";
    }
}

# stat_path($params,$key)
# return string with statistics for superpath file w/id $key
sub stat_path {
    my ($params, $key) = @_;

    my $fnames = $params->{'fnames'};
    my $dir = $params->{'dir'};
    my $fname = "$dir/$fnames->{$key}";

    my ($stats,$err) = mysystem("wc $fname");

    # should be line like
    #      388    2762    9756 filename

    if ($err) {
	return "N/A";
    } elsif ($stats =~ m:(\d+)\s+(\d+)\s+(\d+):) {
	my $npaths = $1 / 4;

	return ($npaths == 1) ?
	    "($npaths superpath)" :
		"($npaths superpaths)" ;
    } else {
	return "N/A";
    }
}

# stat_mate($params,$key)
# return string with statistics for matepair file w/id $key
sub stat_mate {
    my ($params, $key) = @_;

    my $fnames = $params->{'fnames'};
    my $dir = $params->{'dir'};
    my $fname = "$dir/$fnames->{$key}";

    my ($stats,$err) = mysystem("wc $fname");

    # should be line like
    #      888    3552   19117 filename

    if ($err) {
	return "N/A";
    } elsif ($stats =~ m:(\d+)\s+(\d+)\s+(\d+):) {
	my $npairs = $1;
	return ($npairs == 1) ?
	    "($npairs mate-pair)" :
		"($npairs mate-pairs)" ;

    } else {
	return "N/A";
    }
}


# stat_iparams($params,$key)
# return string with input parameters for report file w/id $key
sub stat_iparams {
    my ($params, $key) = @_;

    my $fnames = $params->{'fnames'};
    my $dir = $params->{'dir'};
    my $fname = "$dir/$fnames->{$key}";

    my ($stats,$err) = mysystem("head -20 $fname");

    # should be lines like
    #                             EULER report
    # EULER V2.0, (c) 2001-2004, The Regents of the University of California.
    # ...
    #
    # Fasta file:           Bac99.fasta.screen
    # Quality file:         Bac99.fasta.qual
    # Euler cutoff:         0.02
    # Quality value cutoff: 15
    # Trim:                 phred

    if ($err) {
	return "N/A";
    } elsif ($stats =~ /Trim:\s+(\S+)/) {
	my $trim = $1 ;
	my $euler_cutoff;
	my $qual_cutoff;

	if ( $stats =~ /Euler cutoff:\s+(\S+)/ ) {
	    $euler_cutoff = $1 ;
	}

	if ( $stats =~ /Quality value cutoff:\s+(\S+)/ ) {
	    $qual_cutoff = $1 ;
	}

	my $comment;
	$comment = <<"ENDip_hdr";
<table border=1 width=100%>
<tr><th align=left>Trim</th>                <td>$trim</td></tr>
ENDip_hdr

        if (defined $euler_cutoff) {
	    $comment .= <<"ENDip_ec";
<tr><th align=left>Euler cutoff</th>        <td>$euler_cutoff</td></tr>
ENDip_ec
        }

	if (defined $qual_cutoff) {
	    $comment .= <<"ENDip_qc";
<tr><th align=left>Quality value cutoff</th><td>$qual_cutoff</td></tr>
ENDip_qc
        }

        $comment .= <<"ENDip_tail";
</table>
ENDip_tail

        return $comment;
    } else {
	return "N/A";
    }
}



# mysystem($cmd)
# returns ($stdout,$stderr)
#   $stdout = output on stdout
#   $stderr = output on stderr
# terminates each line with \n  (normal `cmd` wouldn't do that)

sub mysystem {
    my ($cmd) = @_;

    # From Perl FAQ:
    #  (cmd | sed 's/^/STDOUT:/) 2>&1
    # causes STDOUT and STDERR to be merged,
    # but lines of STDOUT are prefixed by "STDOUT:"

    my $output = `($cmd | sed 's/^/STDOUT:/') 2>&1`;

    my ($stdout,$stderr) = ("","");

    foreach my $line (split "\n", $output) {
	if ($line =~ /^STDOUT:(.*)/) {
	    $stdout .= "$1\n";
	} else {
	    $stderr .= "$line\n";
	}
    }

    return ($stdout,$stderr);
}



###############################################################################

# generate HTML index for files in directory

sub print_index {
    my ($params) = @_;

#    print STDERR Dumper($params);

    my $indexfname = $params->{'indexfname'};

    if (!($indexfname =~ m:/:)) {
	$indexfname = "$params->{'dir'}/$indexfname";
    }

    open INDEX, ">$indexfname"    ||   die "Can't create $indexfname";



    print INDEX <<"HTMLtop_end";
<!DOCTYPE HTML PUBLIC "-//IETF//DTD HTML//EN">
<HTML>
<HEAD>
<TITLE>EULER results index: $params->{'xxxx'}</TITLE>
</HEAD>
<body text="$colors::pagetext" bgcolor=$colors::pagebg link="$colors::link" vlink="$colors::vlink" alink="$colors::alink">
<H1>EULER results index: $params->{'xxxx'}</H1>
<b>EULER V2.0, Copyright &copy; 2001-2004 The Regents of the University of California. All Rights Reserved.
<br>See details in the file LICENSE in the source code.
EULER V2.0 web portal users: see <a href="http://nbcr.sdsc.edu/euler">http://nbcr.sdsc.edu/euler</a></b>
HTMLtop_end


    print INDEX "<h2>EULER user files</h2>";
    print INDEX "<center>", brief_table($params), "</center>";


    print INDEX "<h2>EULER developer files</h2>";
    print INDEX "<center>", detailed_table($params), "</center>";


    if ($params->{'missing'}) {
	my $missing_files = missing_files($params);
	if ($missing_files) {
	    print INDEX "<h2>Other files</h2>\n$missing_files";
	}
    }

    print INDEX <<"HTMLtail_end";
</BODY>
</HTML>
HTMLtail_end


    close(INDEX)    ||   die("Can't close $indexfname");
}


###############################################################################
###############################################################################


# read in the EULER report, output it in HTML

# pick out key results to put in summary table

#sub htmlify_report {
#    my ($params) = @_;
#
#    my $repname = $params->{'repname'};
#    my $repnamehtml = $params->{'repnamehtml'};
#
#    open(RAWREP, "<", $repname);
#}


###############################################################################
###############################################################################

main(@ARGV);
