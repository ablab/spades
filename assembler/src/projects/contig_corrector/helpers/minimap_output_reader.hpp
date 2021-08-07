#pragma once
#include "aligner_output_reader.hpp"

namespace helpers {

enum class MColumns : size_t {
    Q_name      = 0,
    Q_size      = 1,
    Q_start     = 2,
    Q_end       = 3,
    strand      = 4,
    T_name      = 5,
    T_size      = 6,
    T_start     = 7,
    T_end       = 8,
    match       = 9,
    bases       = 10,
    quality     = 11,
    TOTAL_COLUMNS_SIZE = 12
};

template<MColumns el>
struct type_getter<MColumns, el> {
    using type = std::conditional_t<
                        el == MColumns::Q_name ||
                        el == MColumns::T_name,
                        std::string,
                 std::conditional_t<
                        el == MColumns::strand,
                        char,
                 long long>>;
};

template<MColumns ... columns>
Records<MColumns, columns ...> ReadMinimapOutput(std::istream & inp, FilterType<MColumns, columns ...> const & filter) {
    Records<MColumns, columns ...> records;
    RecordPusher<MColumns, columns ...> pusher(records, filter);
    std::string line;
    size_t total_lines = 0;
    size_t accepted_lines = 0;
    while (GetNextNonemptyLine(inp, line)) {
        accepted_lines += pusher.Push(line);
        ++total_lines;
    }
    INFO("Total line read: " << total_lines);
    INFO("Accepted lines: " << accepted_lines);
    return records;
}

} // namespace helpers
