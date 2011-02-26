#include <iostream>
#include "toyexamples.hpp"
#include "graph.hpp"
#include "graphBuilder.hpp"

int main() {
  //std::cout << "Hello, A Bruijn!";
  //CGraph graph = GraphBuilder().build();

  ConstructDeBruijnGraph ( "ACTGTACGTACCTGT", 7, 4 );
}
