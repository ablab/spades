#pragma once

#include "standard_vis.hpp"


bool CheckFileExists(const string& file) {
	if (!fs::is_regular_file(file)) {
        ERROR("The file " << file << " does not exist.");
        return false;
    }
    return true;
}
