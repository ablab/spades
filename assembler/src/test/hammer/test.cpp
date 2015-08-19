//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "cute/cute.h"
#include "cute/cute_runner.h"
#include "cute/ide_listener.h"
#include "valid_kmer_generator_test.hpp"

void runSuite() {
  cute::suite s;
  s += ValidKMerGeneratorSuite();
  cute::ide_listener lis;
  cute::makeRunner(lis)(s, "The Suite");
}

int main() {
  runSuite();
  return 0;
}
