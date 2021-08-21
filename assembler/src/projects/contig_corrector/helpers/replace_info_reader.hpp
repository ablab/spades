#pragma once

#include "replace_info_fwd.hpp"
#include "aligner_output_reader.hpp"

#include <string>
#include <fstream>

namespace helpers {

template<ReplaceInfoColumns el>
struct type_getter<ReplaceInfoColumns, el> {
    using type = std::conditional_t<
                        el == ReplaceInfoColumns::from ||
                        el == ReplaceInfoColumns::len,
                        unsigned long long,
                        std::conditional_t<
                            el == ReplaceInfoColumns::seq_type,
                            RangeType,
                            std::string>>;
};

template<>
inline RangeType CastTo<RangeType>(std::string && s) {
    if (s == "origin")
        return RangeType::origin;
    if (s == "edge")
        return RangeType::edge;
    if (s == "path")
        return RangeType::path;
    throw std::invalid_argument("Unknown type value '" + s + "' in replace info; suppoted values: 'origin'/'edge'/'path'");
}

template<ReplaceInfoColumns ... columns>
Records<ReplaceInfoColumns, columns ...> ReadReplaceInfoDump(std::string const & file, FilterType<ReplaceInfoColumns, columns ...> filter) {
    std::ifstream inp(file);
    if (!inp.is_open())
        throw "Cannot open '" + file + "' file";
    Records<ReplaceInfoColumns, columns ...> records;
    RecordPusher<ReplaceInfoColumns, columns ...> pusher(records, filter);
    std::string line;
    GetNextNonemptyLine(inp, line); // skipping header
    while (GetNextNonemptyLine(inp, line))
        pusher.Push(line);
    return records;
}

} // namespace helpers
