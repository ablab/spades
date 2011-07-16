#include "cute/cute.h"
#include "kmer_functions_test.hpp"
#include "ide_listener.h"
#include "cute_runner.h"

void runSuite() {
  cute::suite s;
  s += KMerFunctionsSuite();
  cute::ide_listener lis;
  cute::makeRunner(lis)(s, "The Suite");
}

int main() {
  runSuite();
  return 0;
}
