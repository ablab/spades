#include <iostream>
#include "graph.hpp"
#include "graphBuilder.hpp"

int main() {
	std::cout << "Hello, A Bruijn!";
	CGraph graph = GraphBuilder().build();
}
