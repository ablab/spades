#ifndef _toy_examples_h_
#define _toy_examples_h_

#include <string>
using namespace std;

//typename std::string;

/// Given a genomic string and k, prints 
/// the corresponding de Bruijn graph to a given file.
/// The implementation is highly INEFFICIENT.
extern void ConstructDeBruijnGraph ( string genome, unsigned read_size, unsigned k );

extern void ConstructDeBruijnGraphSimplified ( string genome, unsigned read_size, unsigned k );

extern void ABruijnGraphWithGraphVisualizer ( string genome, unsigned read_size, unsigned k );

#endif

