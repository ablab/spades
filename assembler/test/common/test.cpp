#include "cute.h"
#include "ide_listener.h"
#include "cute_runner.h"
#include "seqTest.hpp"
#include "ireadstreamTest.hpp"
#include "nuclTest.hpp"
#include "ifaststreamTest.hpp"
#include "qualTest.hpp"
#include "onlineGraphVisualizerTest.hpp"
#include "offlineGraphVisualizerTest.hpp"

void runSuite() {
	 cute::suite s;
	 //TODO add your test here
	 s += SeqSuite();
	 s += QualSuite();
	 s += NuclSuite();
	 s += IFastaStreamSuite();
	 s += IReadStreamSuite();
	 s += onlineGraphVisualizerSuite();
	 s += offlineGraphVisualizerSuite();
	 cute::ide_listener lis;
	 cute::makeRunner(lis)(s, "The Suite");
 }

 int main() {
     runSuite();
     return 0;
 }
