#include "cute.h"
#include "ide_listener.h"
#include "cute_runner.h"
#include "seqTest.hpp"
#include "sequenceTest.hpp"
#include "ireadstreamTest.hpp"
#include "quality_read_stream.hpp"
#include "nuclTest.hpp"
#include "ifaststreamTest.hpp"
#include "qualTest.hpp"
#include "onlineGraphVisualizerTest.hpp"
//#include "offlineGraphVisualizerTest.hpp"
#include "similarTest.hpp"
#include "cuckooTest.hpp"
#include "readGeneratorTest.hpp"

void runSuite() {
	 cute::suite s;
	 //TODO add your test here
	 s += SeqSuite();
	 s += SequenceSuite();
	 s += QualSuite();
	 s += NuclSuite();
	 s += IFastaStreamSuite();
	 s += IReadStreamSuite();
	 s += onlineGraphVisualizerSuite();
//	 s += offlineGraphVisualizerSuite();
	 s += similarSuite();
	 s += CuckooSuite();
	 s += ReadGeneratorSuite();
	 cute::ide_listener lis;
	 cute::makeRunner(lis)(s, "The Suite");
 }

 int main() {
     runSuite();
     return 0;
 }
