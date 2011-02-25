#ifndef _toy_examples_h_
#define _toy_examples_h_

#include <string>
using namespace std;

//typename string;

/// Given a genomic string and k, prints 
/// the corresponding de Bruijn graph to a given file.
/// The implementation is highly inefficient.
extern void ConstructDeBruijnGraph ( string genome, unsigned k );

#endif

