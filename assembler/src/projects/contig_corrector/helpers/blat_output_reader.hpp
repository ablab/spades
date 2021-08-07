#pragma once
#include "aligner_output_reader.hpp"

namespace helpers {

enum class BColumns : size_t {
    match       = 0,
    mismatch    = 1,
    rep_match   = 2,
    Ns          = 3,
    Q_gap_count = 4,
    Q_gap_bases = 5,
    T_gap_count = 6,
    T_gap_bases = 7,
    strand      = 8,
    Q_name      = 9,
    Q_size      = 10,
    Q_start     = 11,
    Q_end       = 12,
    T_name      = 13,
    T_size      = 14,
    T_start     = 15,
    T_end       = 16,
    block_count = 17,
    blockSizes  = 18,
    qStarts     = 19,
    tStarts     = 20,
    TOTAL_COLUMNS_SIZE = 21
};

template<BColumns el>
struct type_getter<BColumns, el> {
    using type = std::conditional_t<
                        el == BColumns::Q_name ||
                        el == BColumns::T_name ||
                        el == BColumns::blockSizes ||
                        el == BColumns::qStarts ||
                        el == BColumns::tStarts,
                        std::string,
                 std::conditional_t<
                        el == BColumns::strand,
                        char,
                 long long>>;
};

inline void SkipHeader(std::istream & inp) {
    std::string line;
    if (!GetNextNonemptyLine(inp, line))
        throw std::string("Empty file");
    if (line != "psLayout version 3")
        throw std::string("Sorry, unsupported blat output version");
    GetNextNonemptyLine(inp, line); // column's_names_upper_part
    GetNextNonemptyLine(inp, line); // column's_names_lower_part
    GetNextNonemptyLine(inp, line); // -------------------------
    if (line != std::string(159, '-'))
        throw std::string("blat output file is corrupted");
}

template<BColumns ... columns>
Records<BColumns, columns ...> ReadBlatOutput(std::istream & inp, FilterType<BColumns, columns ...> const & filter) {
    Records<BColumns, columns ...> records;
    RecordPusher<BColumns, columns ...> pusher(records, filter);
    std::string line;
    SkipHeader(inp);
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
