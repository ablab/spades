#!/usr/bin/perl -w
###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2011) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################
#
# 2011-01   Filipe Ribeiro
#

use strict;

use FindBin;
# Local libraries, from elsewhere in the BroadCRD repository.
use lib "$FindBin::Bin/../";
use lib "$FindBin::Bin/";

use PerlRunTime;    # run_or_die
use PerlUtils;      # ISO_date, abort
use ArachneArgs;    # Command-line ARG=VALUE parser

sub Tag { return ISO_date() . " (QSP): "; }



# CONTROL BEGINS HERE
# Parse command-line options of the form KEY=value.
# This function comes from the ArachneArgs.pm module.

my %args = getCommandArguments
    (HEAD    => { value => undef,
                  help  => "Generates qual stats plots for <HEAD>.qualb."},
     RC      => { value => "0",
                  help  => "0: reads are not RC'ed;  1: reads are RC'ed."},
     OUT_DIR => { value => ".",
                  help  => "Directory where plots are to be stored." },
     GIF     => { value => "0",
                  help  => "Wether to convert .eps to .gif." },);


my $out_dir = $args{OUT_DIR};
abort(Tag() . "Can't find directory '$out_dir'.")
    if (!-d $out_dir);


my @arg_heads = array_from_ref_or_value($args{HEAD});
my @arg_rcs = array_from_ref_or_value($args{RC});

my $n = @arg_heads;
if (@arg_rcs < $n) {
    if (@arg_rcs == 1) {
        push @arg_rcs, map $arg_rcs[0], (1..$n-1);
    }
    else {
        abort(Tag() . "Must specify the same number of RCs as HEADs.");
    }
}


my @gp_orig = <DATA>; # get the gnuplot template from the __DATA__ section at the end 


for (my $i = 0; $i != $n; $i++) {

    my $in_full_head = $arg_heads[$i];
    
    my ($in_dir, $in_head) = (($in_full_head =~ /^(.+)\/([^\/]+)$/) ? 
                              ($1, $2) : 
                              (".", $in_full_head));
    my $rc = $arg_rcs[$i] ? "True" : "False";

    my $out_head = $in_head;
    my $out_full_head = "$out_dir/$out_head";

    my $qualb_fn = "$in_full_head.qualb";
    my $stats_fn = "$out_full_head.stats";
    my $out_fn   = "$out_full_head.out";

    if (!-e $stats_fn) 
    {
        print_action(Tag() . "Generating stats for '$qualb_fn'.");
        my $cmd = ("$FindBin::Bin/QualbStats" . 
                   " HEAD=$in_full_head" .
                   " OUT_HEAD=$out_full_head" .
		   " RC=$rc" .
                   " TEE=$out_fn");
        print "$cmd\n";
        print `$cmd`;
    }
    else 
    {
        print_action(Tag() . "Found precomputed stats for '$qualb_fn'.");
    }

    if (-e "$stats_fn" && executable_exists("gnuplot")) {
        my $gp_fn = "$out_head.gp";

        my @gp = @gp_orig;

        my $orientation = $arg_rcs[$i] ? "outward" : "inward";

        foreach my $line (@gp)
        {
            $line =~ s/__ORIENTATION__/$orientation/g;
            $line =~ s/__HEAD__/$out_head/g;
        }
        file_write_or_die("$out_dir/$gp_fn", @gp);

        print_action(Tag(). "Generating plots.");

        system "cd $out_dir ; chmod --quiet +x $gp_fn ; ./$gp_fn";

        if ($args{GIF} && executable_exists("convert")) {
            foreach my $eps_fn (`ls $out_dir/$out_head.*.eps`)
            {
                chomp $eps_fn;
                my $gif_fn = ($eps_fn =~ /^(.+).eps$/)[0] . ".gif";
                print_subaction(Tag(). "Converting '$eps_fn' to '$gif_fn'\n");
                system "convert -density 120 $eps_fn $gif_fn";
            }
        }

    }



}



sub executable_exists
{
    my ($exe) = @_;
    return (grep -x "$_/$exe", (split ":", $ENV{PATH})) ? 1 : 0;
}






sub file_write_or_die
{
    my ($fn, @l) = @_;
    open(FILE, ">$fn") or die "\n**** can't open '$fn' for writing.\n\n";
    print FILE @l;
    close(FILE);
}






__DATA__
#!/usr/bin/env gnuplot

scale = 1.0


set encoding iso_8859_1

#set terminal postscript {landscape | portrait | eps | default}
#                        {enhanced | noenhanced}
#                        {color | monochrome} {solid | dashed}
#                        {<duplexing>}
#                        {"<fontname>"} {<fontsize>}
#
#set terminal postscript eps enhanced monochrome dashed 'Helvetica' 22*scale
#black = 1; grey = 0
set terminal postscript eps noenhanced color solid 'Helvetica' 22 * scale
black = 7; grey = 9



set border 31 lt black lw 4*scale
#set zeroaxis lt 7 lw 1
set grid  lt grey lw scale

size_h    = 0.80 * scale 
margin_l  = 0.25 * scale
margin_ih = 0.05 * scale
margin_r  = 0.10 * scale

size_v    = 1.00 * scale
margin_t  = 0.20 * scale
margin_iv = 0.20 * scale
margin_b  = 0.20 * scale


set lmargin  0; set rmargin  0
set tmargin  0; set bmargin  0


set tics scale -1.0 * scale



# -----------------------------
# mean q as a function of cycle
# -----------------------------

nh = 1; nv = 1

set output '__HEAD__.mean_q.eps'
set size margin_l + nh*size_h +(nh-1)*margin_ih + margin_r, \
         margin_b + nv*size_v +(nv-1)*margin_iv + margin_t 
set origin 0,0
set multiplot

  set title "'__HEAD__' mean quality score (__ORIENTATION__)"
    
  ih = 0; iv = 0;

  unset log x
  set xlabel 'cycle' font 'Helvetica'
  set format x
  set xtics
  
  unset log y
  set ylabel 'quality score' font 'Helvetica'
  set format y
  set ytics

  set bars small

  set grid xtics ytics lt grey
  set key top right
  set size size_h, size_v
  set origin margin_l + ih*(size_h + margin_ih), margin_b+ iv*(size_v + margin_iv)
  plot [0:] [0:50] \
    '__HEAD__.stats' i 2 u ($1+1):3 t '1st read'  w lp  lt 1 pt 7,\
    ''               i 0 u ($1+1):3     not           w l   lt 1 lw 6,\
    ''               i 3 u ($1+1):3 t '2nd read'  w lp  lt 3 pt 7,\
    ''               i 1 u ($1+1):3     not           w l   lt 3 lw 6

unset multiplot


# --------------------------------
# qualities as a function of cycle
# --------------------------------

nh = 2; nv = 1

set output '__HEAD__.cycle_quals.eps'
set size margin_l + nh*size_h +(nh-1)*margin_ih + margin_r, \
         margin_b + nv*size_v +(nv-1)*margin_iv + margin_t 
set origin 0,0
set multiplot

  set title "'__HEAD__' read 1 (__ORIENTATION__)"
    
  ih = 0; iv = 0;

  unset log x
  set xlabel 'cycle' font 'Helvetica'
  set format x
  set xtics
  
  unset log y
  set ylabel 'quality score' font 'Helvetica'
  set format y
  set ytics

  set grid xtics ytics lt grey
  set key top right
  set size size_h, size_v
  set origin margin_l + ih*(size_h + margin_ih), margin_b+ iv*(size_v + margin_iv)
  plot [0:] [0:50] \
    '__HEAD__.stats' i 2 u ($1+1):2 t 'effective q'    w p  lt 1 pt 7,\
    ''               i 0 u ($1+1):2 not                w l  lt 1 lw 6,\
    ''               i 2 u ($1+1):3 t 'average q'      w p  lt 2 pt 7,\
    ''               i 0 u ($1+1):3 not                w l  lt 2 lw 6,\
    ''               i 2 u ($1+1):5 t 'Q1 Q2 Q3 q'     w p  lt 3 pt 7,\
    ''               i 0 u ($1+1):5 not                w l  lt 3 lw 6,\
    ''               i 2 u ($1+1):4 not                w p  lt 3 pt 7 ps 0.7,\
    ''               i 2 u ($1+1):6 not                w p  lt 3 pt 7 ps 0.7



  set title "'__HEAD__' read 2 (__ORIENTATION__)"

  ih = 1; iv = 0;

  unset log x
  set xlabel 'cycle' font 'Helvetica'
  set format x
  set xtics
  
  unset log y
  unset ylabel
  set format y ''
  set ytics

  set grid xtics ytics lt grey
  set key top right
  set size size_h, size_v
  set origin margin_l + ih*(size_h + margin_ih), margin_b+ iv*(size_v + margin_iv)
  plot [0:] [0:50] \
    '__HEAD__.stats' i 3 u ($1+1):2 t 'effective q'    w p  lt 1 pt 7,\
    ''               i 1 u ($1+1):2 not                w l  lt 1 lw 6,\
    ''               i 3 u ($1+1):3 t 'average q'      w p  lt 2 pt 7,\
    ''               i 1 u ($1+1):3 not                w l  lt 2 lw 6,\
    ''               i 3 u ($1+1):5 t 'Q1 Q2 Q3 q'     w p  lt 3 pt 7,\
    ''               i 1 u ($1+1):5 not                w l  lt 3 lw 6,\
    ''               i 3 u ($1+1):4 not                w p  lt 3 pt 7 ps 0.7,\
    ''               i 3 u ($1+1):6 not                w p  lt 3 pt 7 ps 0.7


unset multiplot


# ------------------------------------------
# frequencies as a function of quality score
# ------------------------------------------

nh = 2; nv = 1

set output '__HEAD__.qual_freqs.eps'
set size margin_l + nh*size_h +(nh-1)*margin_ih + margin_r, \
         margin_b + nv*size_v +(nv-1)*margin_iv + margin_t 
set origin 0,0
set multiplot

  set title "'__HEAD__' read 1 (__ORIENTATION__)"

  ih = 0; iv = 0;

  unset log x
  set xlabel 'quality' font 'Helvetica'
  set format x
  set xtics
  
  unset log y
  set ylabel 'frequency' font 'Helvetica'
  set format y
  set ytics

  set grid xtics ytics lt grey
  set key top left
  set size size_h, size_v
  set origin margin_l + ih*(size_h + margin_ih), margin_b+ iv*(size_v + margin_iv)
  plot [0:40] [0:1] \
    '__HEAD__.stats' i 4 u 1:2 t 'reads' w lp lt 1 lw 4 pt 7,\
    ''               i 4 u 1:3 not       w lp lt 1 lw 4 pt 7,\
    ''               i 4 u 1:4 t 'bases' w lp lt 3 lw 4 pt 7,\
    ''               i 4 u 1:5 not       w lp lt 3 lw 4 pt 7


  set title "'__HEAD__' read 2 (__ORIENTATION__)"

  ih = 1; iv = 0;

  unset log x
  set xlabel 'quality' font 'Helvetica'
  set format x
  set xtics
  
  unset log y
  unset ylabel 
  set format y ''
  set ytics

  set grid xtics ytics lt grey
  set key top left
  set size size_h, size_v
  set origin margin_l + ih*(size_h + margin_ih), margin_b+ iv*(size_v + margin_iv)
  plot [0:40] [0:1] \
    '__HEAD__.stats' i 5 u 1:2 t 'reads' w lp lt 1 lw 4 pt 7,\
    ''               i 5 u 1:3 not       w lp lt 1 lw 4 pt 7,\
    ''               i 5 u 1:4 t 'bases' w lp lt 3 lw 4 pt 7,\
    ''               i 5 u 1:5 not       w lp lt 3 lw 4 pt 7

unset multiplot


#guide for line and point styles:
#  0  ..............  .                    broken line
#  1  --------------  +                    red
#  2  -- -- -- -- --  x                    green
#  3  -  -  -  -  -   *                    blue
#  4  ..............  empty square         magenta
#  5  __.__.__.__.__  full  square         cyan
#  6  _ . _ . _ . _   empty circle         yellow
#  7  - -  - -  - -   full  circle         black
#  8  - - -  - - -    empty up triangle    brown
#  9  - - - -  - - -  full  up triangle    grey
# 10 (1)              empty down triangle
# 11 (2)              full  down triangle
# 12 (3)              empty diamond
# 13 (4)              full  diamond
# 14 (5)              empty pentagon
# 15 (6)              full  pentagon
# 16-31               watches



# ------------
#    done
# ------------
