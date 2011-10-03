#pragma once
#include "boost/current_function.hpp"
#include <sstream>

struct assertion_failed_exception : public std::exception {
};

#define verify(expr) 																							\
if (!(expr))? {																									\
	std::stringstream ss;																						\
	ss << "Verification of expression " << #expr << " failed in function " <<  BOOST_CURRENT_FUNCTION << 		\
	". In file " << __FILE__ << " on line " << __LINE__ << "." << std::endl;									\
	std::cout << ss.str();																						\
	std::cerr << ss.str();																						\
	throw assertion_failed_exception;																			\
}
