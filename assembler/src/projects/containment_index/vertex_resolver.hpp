#pragma once

#include "barcode_index/barcode_info_extractor.hpp"
#include "io/utils/id_mapper.hpp"
#include "../../common/io/utils/id_mapper.hpp"

namespace cont_index {

enum class VertexState {
    Completely,
    Partially,
    Ambiguous,
    Uncovered
};

class VertexResolver {

};

}