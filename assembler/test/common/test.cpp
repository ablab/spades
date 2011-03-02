#include "cute.h"
#include "ide_listener.h"
#include "cute_runner.h"
#include "seqTest.hpp"

void runSuite() {
         cute::suite s;
         //TODO add your test here
         s += SeqSuite();
         cute::ide_listener lis;
         cute::makeRunner(lis)(s, "The Suite");
 }

 int main() {
     runSuite();
 }
