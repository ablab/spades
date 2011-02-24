#include <iostream>
#include "graph.h"
#include "graphBuilder.h"

int main() {
	std::cout << "Hello, A Bruijn!";
	CGraph graph = GraphBuilder().build();
}
