#pragma once

namespace error_analyzer {

enum class RangeType {
    origin = 0,
    edge   = 1,
    path   = 2
};

enum class RangeEndsType {
    origin_head = 3,
    origin_tail = 4
};

} // namespace error_analyzer
