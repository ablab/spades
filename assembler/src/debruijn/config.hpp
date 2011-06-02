/*
 * config.hpp
 *
 *  Created on: 03.05.2011
 *      Author: vyahhi
 */

#ifndef CONFIG_HPP_
#define CONFIG_HPP_

#include "libs/ConfigFile/ConfigFile.h"

/*
 * Most of run-time configurations is located in config.inp file (CONFIG_FILE below).
 */
#define CONFIG_FILENAME "./src/debruijn/config.inp"
#define K 35
//#define R 100
//#define I 220

ConfigFile CONFIG(CONFIG_FILENAME);

#endif /* CONFIG_HPP_ */
