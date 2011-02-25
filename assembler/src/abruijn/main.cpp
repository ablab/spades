#include <iostream>
#include "graph.h"
//#include "graphBuilder.h"
#include "toyexamples.h"

int main() {
  std::cout << "Hello, A Bruijn!";
  //CGraph graph = GraphBuilder().build();

  ConstructDeBruijnGraph ( "ACTGTACGTAC", 4 );

}
