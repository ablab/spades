#include "cute.h"
#include "ide_listener.h"
#include "cute_runner.h"
#include "hashTableTest.hpp"
#include "graphConstructionTest.hpp"

void runSuite() {
	cute::suite s;
	//TODO add your test here
	s += HashTableSuite();
	s += CheckStoreVertexSuite();
	s += CheckUniqueWaySuite();
	s += GoUniqueWaySuite();

	cute::ide_listener lis;
	cute::makeRunner(lis)(s, "The Suite");
}

int main() {
	runSuite();
}
