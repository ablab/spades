#include "cute.h"
#include "ide_listener.h"
#include "cute_runner.h"
#include "seqTest.hpp"
#include "sequenceTest.hpp"
#include "qualityTest.hpp"
#include "nuclTest.hpp"
#include "ireadstreamTest.hpp"
#include "onlineGraphVisualizerTest.hpp"
#include "similarTest.hpp"
#include "cuckooTest.hpp"
#include "single_read_test.hpp"
#include "paired_read_test.hpp"

void runSuite() {
  cute::suite s;
  //TODO add your test here
  s += SeqSuite();
  s += SequenceSuite();
  s += QualitySuite();
  s += NuclSuite();
  s += IReadStreamSuite();
  s += onlineGraphVisualizerSuite();
  s += similarSuite();
  s += CuckooSuite();
  s += SingleReadSuite();
  s += PairedReadSuite();
  cute::ide_listener lis;
  cute::makeRunner(lis)(s, "The Suite");
}

int main() {
  runSuite();
  return 0;
}
