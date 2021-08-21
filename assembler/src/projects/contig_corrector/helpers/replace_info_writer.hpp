#pragma once

#include "replace_info_fwd.hpp"

#include <string>
#include <memory>
#include <fstream>

namespace helpers {

class ReplaceInfoWriter {
    std::ofstream out;
    ReplaceInfoWriter(std::string file);
public:
    void Write(std::string const & name, RangeType type, unsigned long long from, unsigned long long len);

    static std::unique_ptr<ReplaceInfoWriter> stream;
    static void SetStream(std::string file);
};


} // namespace helpers
