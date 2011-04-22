#include "toyexamples.hpp"
#include "graphBuilder.hpp"
#include "logging.hpp"

using namespace abruijn;
using namespace std;

LOGGER("a");

int main(int argc, char* argv[]) {
	INFO("Hello, A Bruijn!");

	for (int i = 1; i < argc; ++i) {
		string s(argv[i]);
	}

	GraphBuilder().build();
//	ConstructDeBruijnGraphSimplified ( "ATGCATTGCACTGCA", 6, 3 );
	return 0;
}
