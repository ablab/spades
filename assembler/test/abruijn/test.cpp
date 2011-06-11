#include "cute.h"
#include "ide_listener.h"
#include "cute_runner.h"
#include "hashTest.hpp"
#include "spellgenometest.hpp"

#include "common/logging.hpp"
DECL_PROJECT_LOGGER("at")

void runSuite() {
	 cute::suite s;
	 s += HashSuite();
         s += SpellingGenomeSuite();
	 cute::ide_listener lis;
	 cute::makeRunner(lis)(s, "The Suite");
 }

 int main() {
     runSuite();
     return 0;
 }
