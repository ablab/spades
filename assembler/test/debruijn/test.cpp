#include "cute.h"
#include "ide_listener.h"
#include "cute_runner.h"
#include "condensedGraphTest.hpp"

void runSuite() {
	 cute::suite s;
	 //TODO add your test here
	 s += CondensedGraphSuite();
	 cute::ide_listener lis;
	 cute::makeRunner(lis)(s, "De Bruijn Project Test Suite");
 }

 int main() {
     runSuite();
     return 0;
 }
