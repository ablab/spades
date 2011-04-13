#include <iostream>
#include "toyexamples.hpp"
#include "graph.hpp"
#include "graphBuilder.hpp"
#include "logging.hpp"

LOGGER("a");

int main() {
	INFO("Hello, A Bruijn!");
	GraphBuilder().build();

//	ConstructDeBruijnGraphSimplified ( "ACTGACTGTTGACACTG", 9, 5 );
//	ConstructDeBruijnGraphSimplified ( "ATTGGTACATTGTGGTACGTACTGACT", 11, 5 );
	return 0;
}
