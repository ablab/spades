#include "cute.h"
#include "ide_listener.h"
#include "cute_runner.h"
#include "hashTableTest.hpp"

void runSuite() {
<<<<<<< HEAD
         cute::suite s;
         //TODO add your test here
         s += HashTableSuite();
         cute::ide_listener lis;
         cute::makeRunner(lis)(s, "The Suite");
 }
=======
	cute::suite s;
	//TODO add your test here
	s += HashTableSuite();
	s += CheckStoreVertexSuite();
	s += CheckUniqueWaySuite();
	s += GoUniqueWaySuite();

	cute::ide_listener lis;
	cute::makeRunner(lis)(s, "The Suite");
}
>>>>>>> c9e9a4c50cceaea6b9b70c96835b90cb2ee4055d

 int main() {
     runSuite();
 }
