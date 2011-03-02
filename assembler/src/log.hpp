/*
 * log.hpp
 *
 *  Created on: 02.03.2011
 *      Author: vyahhi
 */

#ifndef LOG_HPP_
#define LOG_HPP_

template <size_t N, size_t base = 2>
struct log {
	const static size_t value = 1 + log<N/base, base>::value;
};

template <size_t base>
struct log<1, base> {
	const static size_t value = 0;
};

template <size_t base>
struct log<0, base> {
	const static size_t value = 0;
};

#endif /* LOG_HPP_ */
