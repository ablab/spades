#include "cute.h"
#include "ide_listener.h"
#include "cute_runner.h"
#include "hashTest.hpp"

void runSuite() {
	 cute::suite s;
	 s += HashSuite();
	 cute::ide_listener lis;
	 cute::makeRunner(lis)(s, "The Suite");
 }

 int main() {
     runSuite();
     return 0;
 }
