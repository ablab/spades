#include <boost/interprocess/allocators/allocator.hpp>
#include <boost/version.hpp>
#include <iostream>

// See 'misc/interprocess-test.cpp' for explanation.
#if BOOST_VERSION < 104800
#  error need at least Boost 1.48
#endif

int
main ()
{ }
