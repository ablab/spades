#include <iostream>
#include "toyexamples.hpp"
#include "graph.hpp"
#include "graphBuilder.hpp"
#include "logging.hpp"

LOGGER("a");

int main() {
	INFO("Hello, A Bruijn!");
	GraphBuilder().build();

//	ConstructDeBruijnGraphSimplified ( "TAAACGAAAC", 6, 3 );
	return 0;
}
