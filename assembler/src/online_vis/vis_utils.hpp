#pragma once

#include "standard_vis.hpp"


bool CheckFileExists(const string& file) {
	if (!fs::exists(file)) {
        ERROR("The file " << file << " does not exist.");
        return false;
    }
    return true;
}
