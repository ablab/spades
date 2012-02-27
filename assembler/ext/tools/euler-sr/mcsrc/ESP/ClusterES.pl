#!/usr/bin/env perl

if ($#ARGV != 1) {
  print "usage: $0 length_file expected_fosmid_size\n";
  print "find clusters of end sequences with technique by Ben Raphael:\n";
  print "if |x_2 - x_1| + |y_2 - y_1| < L, then group es_1 and es_2\n";
  exit(0);
}

$in = shift @ARGV;
$length = shift @ARGV;

open(IN, "$in") or die "cannot open $in\n";

@coords = ();
while(<IN>) {
  @vals = split(/\s+/, $_);
  if ($vals[0] < $vals[2]) {
    $first = $vals[0];
    $second = $vals[2];
  }
  else {
    $first = $vals[2];
    $second = $vals[0];
  }
  push  @coords, [$first, $second];
}

# cluster indices holds the index of the cluster of each point
%clusterIndices = ();

# clusters is an array of arrays that contains each cluster


$cluster = 0;
$np = 0;
$nod = 0;
$nc = 0;
for ($i = 0; $i < $#coords;  $i++ ) {
  $ipairdist = $coords[$i][1] - $coords[$i][0];
  if ($ipairdist > $length) {
#    print "@{$coords[$i]}\n";
    $nod++;
    for ($j = $i+1; $j <= $#coords; $j++) {
      $jpairdist=  $coords[$j][1] - $coords[$j][0];
      if ($jpairdist > $length) {
	# Make sure we're not considering a bes along the diagonal
	$len = abs($coords[$j][0] - $coords[$i][0]) + abs($coords[$j][1] - $coords[$i][1]);
	$nc++;
#	print "$len\n";

	if ($len < $length) {
	$np++;
#	  print "found $i, $j, of length $len\n";
	$ci1 = -1;
	$ci2 = -1;
	  if (exists ($clusterIndices{$i})) {
	    $ci1 = $clusterIndices[$i];
	  }
	  if (exists ($clusterIndices{$j})) {
	    $ci2 = $clusterIndices[$j];
	  }
	  # case 1, neither belongs in a cluster, create a new one.
	  if ($ci1 == -1 and $ci2 == -1) {
#	    print "creating new cluster $i, $j\n";
	    $clusterIndices{$i} = $cluster;
	    $clusterIndices{$j} = $cluster;
#	    push @clusters, [($i, $j)];
	    $cluster++;
	  }
	  # case 2, i is already in a cluster, add j to it
	  elsif ($ci1 != -1 and $ci2 == -1) {
#	    print "add $j to $clusterIndices{$i}\n";
	    $clusterIndices{$j} = $clusterIndices{$i};
#	    push @{$clusters[$clusterIndices{$i}]}, $j;
	  }
	  # case 3, j is already in a cluster
	  elsif ($ci1 == -1 and $ci2 != -1) {
#	    print "add $i to $clusterIndices{$j}\n";
	    $clusterIndices{$i} = $clusterIndices{$j};
#	    push @{$clusters[$clusterIndices{$j}]}, $i;
	  }
	  # case 4, both i and j are in clusters, 
	  # merge the two clusters
	  else {
	    if ($clusterIndices{$i} != $clusterIndices{$j} ) {
#	      print "merging clusters: $clusterIndices{$i} $clusterIndices{$j} $#{$clusters[$i]} $#{$clusters[$j]}\n";
	      $ci = $clusterIndices{$i};
	      $cj = $clusterIndices{$j};
	      # merge the two clusters.  This is an expensive operation, but it shouldn't happen too often
	      foreach $index (keys %clusterIndices) {
		if ($clusterIndices{$index} == $cj) {
		  $clusterIndices{$index} = $ci;
		}
	      }
#	      print "dropping cluster: $cj\n";
	      #	    print "size of $clusterIndices{$i} $#{$clusters[$ci]}\n";
	    }
	  }
	}
      }
    }
  }
}


# now transform the cluster indices into clusters

for ($i = 0; $i < $cluster; $i++ ){
  push @clusters, [()];
}

foreach $index (keys %clusterIndices) {
  push @{$clusters[$clusterIndices{$index}]}, $index;
}

# print the clusters;
$nclust =  $#clusters + 1;

#print "ending with nclust: $nclust should have created $cluster\n";
# find the max clustersize
$maxcs = 0;
for ($i = 0; $i < $nclust; $i++) {
  if ($#{$clusters[$i]} > $maxcs) {
    $maxcs = $#{$clusters[$i]};
  }
}

$nclust = scalar @clusters;
#print "nclust now: $nclust\n";
for ($i = 0; $i < $nclust; $i++) {
  for ($j = 0; $j <= $#{$clusters[$i]}; $j++ ) {
    $c1 = $clusters[$i][$j] + 1;
    print "$c1 ";
  }
  while ($j <= $maxcs) {
    print  "-1 ";
    $j++;
  }
  print "\n";
}

