#!/usr/bin/env perl

if ($#ARGV < 0) {
  print "usage: $0 infile\n";
    print "infile is one tuple per line, with possible extra junk (multiplicity, etc)\n";
}
$infile = shift @ARGV;

open (IN, "$infile") or die "cannot open $infile\n";

%vertices = {};
$edges = {};
while (<IN>) {
  /([actg]+)/i;
  $tuple = $1;
  $tuple =~ tr/A,C,T,G/a,c,t,g/;
  $edges{$tuple} = 1;
  ($l,$r) = GetVertices($tuple);
  $vertices{$l} = 1;
  $vertices{$r} = 1;
}

$numEdges = CountOnes(\%edges);
$numVertices = CountOnes(\%vertices);

while ($numEdges > 0) {
  @k = keys(%edges);
  $i = 0;
  $l = $#k+1;
  $isrepeat = 0;

  while ($i <= $#k && $edges{@k[$i]} == 0) {
    $i++;
  }
  if ($i > $#k) {
    print "error, there are $numEdges edges, but didn't find any $l \n";
    exit(0);
  }
  
  $edge = @k[$i];
  $firstEdge = $edge;
  $edges{$edge} = 0;
  $string = $edge;
  # follow link in left direction.
  ($lv, $rv) = GetVertices($edge);
  #  print "$edge got vertices: $lv, $rv\n";
  @leftEdges = ();
  GetLeftEdges($lv, \%edges, \@leftEdges);
  while (($#leftEdges+1) == 1) {
    # print "loop got left edges @leftEdges , @leftEdges[0]\n";
    $edge = @leftEdges[0];
    $edges[$edge] = 0;
    $string = substr(@leftEdges[0], 0, 1) . $string;
    # print "got seq: $string\n";
    ($lv, $rv) = GetVertices($edge);
    # print "$edge got vertices: $lv, $rv\n";
    GetLeftEdges($lv, \%edges, \@leftEdges);
    $l = length(@leftEdges);
  }
  
  if (($#leftEdges+1) > 1) {
    $isrepeat = 1;
  }

  # follow link in right direction
  $edge = $firstEdge;
  ($lv, $rv) = GetVertices($edge);
  GetRightEdges($rv, \%edges);
  while (($#rightEdges + 1) == 1) {
    $edge = @rightEdges[0];
    $edges[$edge] = 0;
    $string = $string . substr(@rightEdges[0], 0, 1);
    ($lv, $rv) = GetVertices($edge);
    GetRightEdges($rv, \%edges, \@rightEdges);
  }
  
  $ls = length($string);
  
  if (($#rightEdges + 1) > 1) {
    $isrepeat = 1;
  }

  if ($isrepeat == 1) {
    print "$ls $string $i (repeat)\n";
  }
  else {
    print "$ls $string $i\n";
  }
  


  $numEdges = CountOnes(\%edges);
}



# now construct the sequences

sub CountOnes {
  my ($h) = @_;
  $numOnes = 0;
  foreach $ky (keys(%{$h})) {
    if ($$h{$ky} == 1) {
      $numOnes++;
    }
  }
  return $numOnes;
}

sub GetVertices {
  my ($edge) = @_;
  $len = length($edge) -1;
  $edge =~ /[actg]([actg]{$len})/i;  
  $r = $1;
  $edge =~ /([actg]{$len})[actg]/i;
  $l = $1;
  return ($l,$r);
}

sub GetLeftEdges {
  my ($vertex, $edges, $leftEdges) = @_;
  my @nucs = ("a","c","t","g");
  
  @{$leftEdges} = ();
  for $nuc (@nucs) {
    if (exists(${$edges}{"$nuc$vertex"})) {
      push @{$leftEdges}, "$nuc$vertex";
    }
  }
}

sub GetRightEdges {
  my ($vertex, $edges, $rightEdges) = @_;
  @nucs = ["a","c","t","g"];
  @{$rightEdges} = ();
  for $nuc (@nucs) {
    if (exists(${$edges}{"$vertex$nuc"})) {
      push @{$rightEdge}, "$vertex$nuc";
    }
  }
}
