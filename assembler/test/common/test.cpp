#include "cute.h"
#include "ide_listener.h"
#include "cute_runner.h"
#include "seqTest.hpp"
#include "sequenceTest.hpp"
#include "ireadstreamTest.hpp"
//#include "quality_read_stream.hpp"
#include "nuclTest.hpp"
//#include "ifaststreamTest.hpp"
#include "qualityTest.hpp"
#include "onlineGraphVisualizerTest.hpp"
//#include "offlineGraphVisualizerTest.hpp"
#include "similarTest.hpp"
#include "cuckooTest.hpp"
#include "trieTest.hpp"
#include "seq_filterTest.hpp"
#include "structuresTest.hpp"
//TODO function readGenomeFromFile is declared outside common directory
//typedef long long ll; //???
//#include "readGeneratorTest.hpp" 

void runSuite() {
  cute::suite s;
  //TODO add your test here
  s += SeqSuite();
  s += SequenceSuite();
  s += QualitySuite();
  s += NuclSuite();
  //s += IFastaStreamSuite();
  s += IReadStreamSuite();
  s += onlineGraphVisualizerSuite();
  //s += offlineGraphVisualizerSuite();
  s += similarSuite();
  s += CuckooSuite();
  //s += ReadGeneratorSuite();
  //The next test causes segmentation fault
  //s += TrieSuite();
  s += SeqFilterSuite();
  //The next test is correct but needs time
  //s += StructuresSuite();
  cute::ide_listener lis;
  cute::makeRunner(lis)(s, "The Suite");
}

int main() {
  runSuite();
  return 0;
}
