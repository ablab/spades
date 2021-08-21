#pragma once

namespace helpers {

enum class ReplaceInfoColumns{
    contig_name        = 0,
    seq_type           = 1,
    from               = 2,
    len                = 3,
    TOTAL_COLUMNS_SIZE = 4
};

enum class RangeType {
    origin = 0,
    edge   = 1,
    path   = 2
};

} // namespace helpers
