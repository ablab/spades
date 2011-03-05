#include "cute.h"
#include "ide_listener.h"
#include "cute_runner.h"
#include "seqTest.hpp"
#include "sequenceTest.hpp"
#include "ireadstreamTest.hpp"
#include "nuclTest.hpp"
#include "ifaststreamTest.hpp"
#include "qualTest.hpp"
#include "onlineGraphVisualizerTest.hpp"
#include "offlineGraphVisualizerTest.hpp"
#include "similarTest.hpp"

void runSuite() {
	 cute::suite s;
	 //TODO add your test here
	 s += SeqSuite();
	 s += SequenceSuite();
	 s += QualSuite();
	 s += NuclSuite();
	 s += IFastaStreamSuite();
	 //s += IReadStreamSuite();
	 s += onlineGraphVisualizerSuite();
	 s += offlineGraphVisualizerSuite();
	 s += similarSuite();
	 cute::ide_listener lis;
	 cute::makeRunner(lis)(s, "The Suite");
 }

 int main() {
     runSuite();
     return 0;
 }
