#!/usr/bin/env perl

$cwd = $ENV{"PWD"};

@dirs = glob("fragment.*");
foreach $dir (@dirs) {
  if ($dir =~ /\.(\d+)/) {
    $fragment = $1;
    print "cd $cwd/$dir/repeat_graphs; ln -sf ../*.fasta .;  ~/projects/mcsrc/assembly/FindResolvingRepeatGraph.pl $fragment.fasta 20 10 > $fragment.repeat_graph_sizes.txt \n";
  }
}
