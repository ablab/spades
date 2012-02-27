$numSlots =$nproc;
while (@jobs) {
  while ($numSlots && @jobs) {
    --$numSlots;
    $nj = @jobs;
    unless (fork()) {
      $job = shift @jobs;
      print "child running at slot $numSlots\n";
      print "$job\n";
      exec $job;
    }
    shift @jobs;
  }
  wait;
  ++$numSlots;
}
while ($numSlots < $nproc) {
  wait;
  ++$numSlots;
}
