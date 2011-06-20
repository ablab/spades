/*
 * logging.hpp
 *
 *  Created on: 02.03.2011
 *      Author: Mikhail Dvorkin
 */

#ifndef LOGGING_HPP_
#define LOGGING_HPP_

#include <iostream>
#include "log4cxx/logger.h"

log4cxx::LoggerPtr __scope_logger();

#define DECL_PROJECT_LOGGER(category)						                \
log4cxx::LoggerPtr __scope_logger()								            \
{															                \
	static log4cxx::LoggerPtr logger = log4cxx::Logger::getLogger(category);\
	return logger;											                \
}

#define DECL_LOGGER(category)                                         \
static log4cxx::LoggerPtr __scope_logger()					          \
{																      \
	static log4cxx::LoggerPtr logger = log4cxx::Logger::getLogger(    \
	        ::__scope_logger()->getName() + "." +	category);        \
	return logger;													  \
}

#define DEBUG(message)                      LOG4CXX_DEBUG(__scope_logger(), message)
#define INFO(message)                       LOG4CXX_INFO (__scope_logger(), message)
#define VERBOSE(n, message)                 if (n % 10000 == 0 && n > 0) INFO(n << message)
#define TRACE(message)                      LOG4CXX_TRACE(__scope_logger(), message)
#define WARN(message)                       LOG4CXX_WARN (__scope_logger(), message)
#define ERROR(message)                      LOG4CXX_ERROR(__scope_logger(), message)
#define LOG_ASSERT(condition, message)      LOG4CXX_ASSERT(__scope_logger(), condition, message)
#define FATAL_ASSERT(condition, message)    LOG4CXX_ASSERT(__scope_logger(), condition, message) std::cerr << "ASSERTION FAILED: " << message; exit(1)
#define FATAL(message)                      LOG4CXX_FATAL(__scope_logger(), message) std::cerr << "FATAL ERROR: " << message; exit(1)

#endif /* LOGGING_HPP_ */
