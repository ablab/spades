/*
 * simple_tools.cpp
 *
 *  Created on: 20.06.2011
 *      Author: vyahhi
 */

#include "common/simple_tools.hpp"
#include "common/logging.hpp"
#include <fstream>

bool fileExists(std::string filename) {
	return std::ifstream(filename);
}

void checkFileExistenceFATAL(std::string filename) {
	if (!fileExists(filename)) {
		FATAL("File " << filename << " doesn't exists or can't be read!\n");
	}
}
