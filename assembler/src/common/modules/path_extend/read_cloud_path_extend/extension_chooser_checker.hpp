#pragma once

#include "modules/path_extend/extension_chooser.hpp"

namespace path_extend {
    class TenXExtensionChecker {
        TenXExtensionChooser& chooser_;
        string genome_file_;

        TenXExtensionChecker(const TenXExtensionChooser& chooser, const string& genome_file) :
                chooser_(chooser), genome_file_(genome_file) {}
    };

}

