//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once
#include "boost/current_function.hpp"
#include <sstream>
#include <iostream>
#include "stacktrace.hpp"

struct assertion_failed_exception : public std::exception {
};

#define VERIFY(expr) 																							\
if (!(expr)) {																									\
	std::stringstream ss;																						\
	ss << "Verification of expression '" << #expr << "' failed in function '" <<  BOOST_CURRENT_FUNCTION << 	\
	"'. In file '" << __FILE__ << " on line " << __LINE__ << "'." << std::endl;									\
	std::cout << ss.str();																						\
	std::cerr << ss.str();																						\
	print_stacktrace();																						    \
	throw((assertion_failed_exception()));																		\
}

#define VERIFY_MSG(expr, msg) 																					\
if (!(expr)) {																									\
	std::stringstream ss;																						\
	ss << "Verification of expression '" << #expr << "' failed in function '" <<  BOOST_CURRENT_FUNCTION << 	\
	"'. In file '" << __FILE__ << "' on line " << __LINE__ << ". Message '" << msg << "'." << std::endl;		\
	std::cout << ss.str();																						\
	std::cerr << ss.str();																						\
	throw ((assertion_failed_exception()));																		\
}
