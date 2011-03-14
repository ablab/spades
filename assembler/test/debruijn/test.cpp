#include "cute.h"
#include "ide_listener.h"
#include "cute_runner.h"
#include "condensedGraphTest.hpp"
#include "debruijnGraphTest.hpp"

void runSuite() {
	 cute::suite s;
	 //TODO add your test here
	 s += CondensedGraphSuite();
	 s += DeBruijnGraphSuite();
	 cute::ide_listener lis;
	 cute::makeRunner(lis)(s, "De Bruijn Project Test Suite");
 }

 int main() {
     runSuite();
     return 0;
 }
