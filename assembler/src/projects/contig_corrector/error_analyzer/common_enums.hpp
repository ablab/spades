#pragma once
#include "helpers/replace_info_fwd.hpp"

namespace error_analyzer {

using RangeType = helpers::RangeType;

enum class RangeEndsType {
    origin_head = 3,
    origin_tail = 4
};

enum class BoundStatType {
    on_bound = 3
};

enum class ErrorType {
    mismatch,
    insertion,
    deletion
};

} // namespace error_analyzer
