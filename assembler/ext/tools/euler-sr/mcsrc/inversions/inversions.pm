package inversions;
use common;

sub PrintInversionList {
  my ($invList, $outfile) = @_;
  open(OUT, ">$outfile") or die "cannot open $outfile\n";
  foreach $ref (common::OrderedKeys(\%{$invList})) {
    foreach $qry (common::OrderedKeys(\%{$$invList{$ref}})) {
      print OUT "$ref $qry\n";
      for ($i = 0; $i <= $#{$$invList{$ref}{$qry}}; $i++ ) {
	print OUT "$$invList{$ref}{$qry}[$i][0] $$invList{$ref}{$qry}[$i][0] " .
	  "$$invList{$ref}{$qry}[$i][0] $$invList{$ref}{$qry}[$i][0]\n";
      }
    }
  }
  close OUT;
}

sub RemoveSpecies {
  my ($invList, $specList) = @_;
  my $deleted = 0;
  foreach $ref (common::OrderedKeys(\%{$invList})) {
    $deleted = 0;
    for ($s = 0; $s <= $#$specList; $s++) {
      if ($ref eq $$specList[$s]) {
	delete $$invList{$ref};
	$deleted = 1;
      }
    }
    if ($deletd == 0) {
      foreach $qry (common::OrderedKeys(\%{$$invList{$ref}})) {
	$deleted = 0;
	for ($s = 0; $s <= $#$specList; $s++ ){
	  if ($$specList[$s] == $qry) {
	    delete $$invList{$ref}{$qry};
	    $deleted = 1;
	  }
	}
      }
    }
  }
}  
    

sub ReadInversionList {
  my ($invFile, $invList) = @_;
  open (INV, $invFile) or die "cannot open $invFile\n";
  my $ref = "";
  my $qry = "";
  $numStored = 0;
  while (<INV>) {
    $line = $_;
    if ($line =~ /(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/) {
      if ($ref eq "" or $qry eq "" ) {
	print "Invalid inversion file format, coordinates specified without species \n";
	exit(1);
      }
      my @coords = split(/\s+/, $line);
      my $ninv = StoreLine($invList, $ref, $qry, \@coords);
      ++$numStored;
    }
    elsif ($line =~ /(\S+)\s+(\S+)/) {
      $ref = $1;
      $qry = $2;
    }
  }
  close INV;
}

sub StoreLine {
  my ($invList, $ref, $qry, $coords) = @_;
  if (! exists $$invList{$ref}) {
    %{$$invList{$ref}} = ();
  }
  if (! exists $$invList{$ref}{$qry}) {
    @{${$invList}{$ref}{$qry}} = ();
  }
  push @{${$invList}{$ref}{$qry}}, [$$coords[0], $$coords[1], $$coords[2], $$coords[3], -1, -1];
  $ninv = scalar @{${$invList}{$ref}{$qry}};
  return $ninv;
}

sub GetSpeciesList {
  my ($invList, $species) = @_;
  my %specHash;
  my $ref, $qry;
  foreach $ref (common::OrderedKeys(\%{$invList})) {
    foreach $qry (common::OrderedKeys(\%{$$invList{$ref}})) {
      $specHash{$qry} = 1;
    }
    $specHash{$ref} = 1;
  }
  foreach $spec (common::OrderedKeys(\%specHash)) {
    push @$species, $spec;
  }
}


sub CreateVertices {
  my ($invPairs, $invSpec, $invSites) = @_;
  foreach $spec (@$invSpec) {
    # process as ref
    if (exists ${$invPairs}{$spec}) {
      my $qry;
      foreach $qry (common::OrderedKeys(\%{$$invPairs{$spec}})) {
	my $i;
	my $invSite;
	for ($i = 0; $i <= $#{$$invPairs{$spec}{$qry}}; $i++) {
	  $invSite = StoreInvSite(\@{$$invSites{$spec}},
				  $$invPairs{$spec}{$qry}[$i][0],
				  $$invPairs{$spec}{$qry}[$i][1]);
	  $$invPairs{$spec}{$qry}[$i][4] = $invSite;
	}
      }
    }
    # process as qry
    my $ref;
    foreach $ref (@$invSpec) {
      if (exists $$invPairs{$ref}{$spec}) {
	for ($i = 0; $i <= $#{$$invPairs{$ref}{$spec}}; $i++) {
	  $invSite = StoreInvSite(\@{$$invSites{$spec}},
				  $$invPairs{$ref}{$spec}[$i][2],
				  $$invPairs{$ref}{$spec}[$i][3]);
	  $$invPairs{$ref}{$spec}[$i][5] = $invSite;
	}
      }
    }
  }
}

sub CountVertices {
  my ($vertices, $species)  = @_;
  $nvertex = 0;
  foreach $spec (@$species) {
    $nvertex += scalar @{$$vertices{$spec}};
  }
  return $nvertex;
}

sub StoreInvSite {
  my ($invSites, $start, $end) = @_;
  my $s;
  $s = 0;
#  print " site: $start $end ";
  while ($s <= $#$invSites) {
    my $ovp = common::ComputeOverlap($$invSites[$s][0],
				     $$invSites[$s][1],
				     $start, $end);
    if ($ovp > 0) {
      my $la = $$invSites[$s][1] - $$invSites[$s][0];
      my $lb = $end - $start;
      last;
    }
    $s++;
  }
  if ($s > $#$invSites) {
    push @$invSites, [ ($start, $end, 0) ];
  }
  else {
#    print " overlaps with $$invSites[$s][0] $$invSites[$s][1] $s\n";
    if ($start < $$invSites[$s][0]) {
      $$invSites[$s][0] = $start;
    }
    if ($end > $$invSites[$s][1]) {
      $$invSites[$s][1] = $end;
    }
  }
  return $s;
}


sub InitOrthInvGraph {
  my ($invSites, $specList, $orthInv) = @_;
  foreach $spec (@$specList) {
    # initialize the hash
    @{$orthInv{$spec}} = ();
    for ( $invIndex = 0; $invIndex <= $#{$$invSites{$spec}}; $invIndex++) {
      push @{$orthInv{$spec}}, ();
    }
  }
}

sub CreateEdges {
  my($invPairs, $specList, $invVertex, $invEdge) = @_;
  my $ref, $qry;
  foreach $ref (common::OrderedKeys(\%$invPairs)) {
    foreach $qry (common::OrderedKeys(\%{$$invPairs{$ref}})) {
      for ($i = 0; $i <= $#{$$invPairs{$ref}{$qry}}; $i++ ) {
	$refInvIndex = $$invPairs{$ref}{$qry}[$i][4];
	$qryInvIndex = $$invPairs{$ref}{$qry}[$i][5];
	# Create edges: 
	#   ( {$ref, $refInvIndex}, {$qry}{$qryInvIndex}), and
	#   ( {$qry, $qryInvIndex}, {$ref}{$refInvIndex}), and
	# add the edges
	# each edge has a field of the species it goes to,
	# the index of the inversion for the species,
	# and a marker for whether or not the edge has been traversed
	my $edgePos = FindEdge(\@{$$invEdge{$ref}[$refInvIndex]}, $qry, $qryInvIndex);
	if ( $edgePos == -1) {
	  push @{$$invEdge{$ref}[$refInvIndex]}, [ ($qry, $qryInvIndex, 0)];
#	  print "made edge ref: (($ref, $refInvIndex), ($qry, $qryInvIndex))\n";
	}
#	else {
#	  print "ref edge (($ref, $refInvIndex), ($qry, $qryInvIndex)) exists at $edgePos\n";
#	}
	$edgePos = FindEdge(\@{$$invEdge{$qry}[$qryInvIndex]}, $ref, $refInvIndex);
	if ( $edgePos == -1) {
	  push @{$$invEdge{$qry}[$qryInvIndex]}, [ ($ref, $refInvIndex, 0)];
#	  print "made edge qry: (($qry, $qryInvIndex), ($ref, $refInvIndex))\n"; 
	}
#	else {
#	  print "qry edge:  (($qry, $qryInvIndex), ($ref, $refInvIndex)) exists at $edgePos\n";
#	}
	  
      }
    }
  }
}

sub FindEdge {
  my ($edgeList, $spec, $index) = @_;
  for ($e = 0; $e <= $#{$edgeList}; $e++ ) {
    if ($$edgeList[$e][0] eq $spec and
	$$edgeList[$e][1] == $index) {
      return $e;
    }
  }
  return -1;
}

sub StoreConnectedComponents {
  my ($invVertex, $invEdge, $specList, $components, $vertexList) = @_;
  my $ref, $qry;
  my @component;
  my @vertices;
  foreach $spec (@$specList) {
    for ($inv = 0; $inv < $#{$$invVertex{$spec}}; $inv++) {
      if ($$invVertex{$spec}[$inv][2] == 0) {
	@component = ();
	@vertices = ();
	StoreConnectedComponent($invVertex, $invEdge, $spec, $inv, \@component, \@vertices);
	push @$components, [@component];
	push @$vertexList, [@vertices];
      }
    }
  }
}

sub PrintConnectedComponents {
  my ($invVertex, $invEdge, $specList, $outFile) = @_;
  my $ref, $qry;
  my @component;
  my @vertices;
  $comp = 0;
  foreach $spec (@$specList) {
    for ($inv = 0; $inv < $#{$$invVertex{$spec}}; $inv++) {
      if ($$invVertex{$spec}[$inv][2] == 0) {
	$outName = $outFile . "$comp.locus";
	open (OUT, ">$outName") or die "cannot open $outName\n";
	@component = ();
	@vertices  = ();
#	print "$comp: \n";
	StoreConnectedComponent($invVertex, $invEdge, $spec, $inv,
				\@component, \@vertices);
	my $c;
	UnionBinCoordinates(\@component);
	for ($c = 0; $c <= $#component; $c++ ) {
	  print OUT "$component[$c][0] $component[$c][1] $component[$c][2]\n";
	}
	close OUT;
	$comp++;
      }
    }
  }
}

sub StoreConnectedComponent {
  my($invVertex, $invEdge, $spec, $invIndex, $component, $vertices, $padding) = @_;
  # if this vertex has already been visited, just quit
  if ($$invVertex{$spec}[$invIndex][2] == 1) {
#    print "$padding   skipping ($spec, $invIndex)\n";
    return;
  }

  my $adjlistsize = scalar @{$$invEdge{$spec}[$invIndex]};

  # otherwise, mark this vertex as visited
  $$invVertex{$spec}[$invIndex][2] = 1;
  push @$component, [($spec,
		     $$invVertex{$spec}[$invIndex][0],
		     $$invVertex{$spec}[$invIndex][1]) ];
  push @$vertices, [($spec, $invIndex)];
  my $e;
  $nedge = scalar @{$$invEdge{$spec}[$invIndex]};
#  print "$padding vertex ($spec, $invIndex) has $nedge edges: ";
#  for ($e = 0; $e <= $#{$$invEdge{$spec}[$invIndex]}; $e++ ) {
#    my $qry = $$invEdge{$spec}[$invIndex][$e][0];
#    my $qryInvIndex = $$invEdge{$spec}[$invIndex][$e][1];
#    print " ($qry, $qryInvIndex)";
#  }
#  print "\n";
  
  for ($e = 0; $e <= $#{$$invEdge{$spec}[$invIndex]}; $e++ ) {
    # if this edge has not been traversed, AND the 
    # connected vertex has not been visited, store 
    # the component at that vertex
    if ($$invEdge{$spec}[$invIndex][$e][2] == 0) {
      # mark this edge as traversed
      $$invEdge{$spec}[$invIndex][$e][2] = 1;
      my $qry = $$invEdge{$spec}[$invIndex][$e][0];
      my $qryInvIndex = $$invEdge{$spec}[$invIndex][$e][1];
#      print "$padding   traversing edge ($qry, $qryInvIndex) \n";
      StoreConnectedComponent($invVertex, $invEdge,
			      $qry, $qryInvIndex, $component, $vertices, " $padding");
    }
  }
}

sub UnionBinCoordinates {
  my ($component) = @_;
  my %compStart, %compEnd;
  for ($c = 0; $c <= $#{$component}; $c++ ) {
    if (! exists $compStart{$$component[$c][0]}) {
      $compStart{$$component[$c][0]} = $$component[$c][1];
      $compEnd{$$component[$c][0]}   = $$component[$c][2];
    }
    else {
      if ( $compStart{$$component[$c][0]} > $$component[$c][1]) {
	$compStart{$$component[$c][0]} = $$component[$c][1];
      }
      if ( $compEnd{$$component[$c][0]} < $$component[$c][2]) {
	$compEnd{$$component[$c][0]} = $$component[$c][2];
      }
    }
  }
  @{$component} = ();
  foreach $c (common::OrderedKeys(\%compStart)) {
    push @$component, [ ($c, $compStart{$c}, $compEnd{$c}) ];
  }
}

return 1;
