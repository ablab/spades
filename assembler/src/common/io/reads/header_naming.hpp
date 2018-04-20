//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <string>
#include <cstdlib>

namespace io {

inline std::string MakeContigId(size_t number, const std::string& prefix = "NODE") {
    return prefix.empty() ? std::to_string(number) : (prefix + "_" + std::to_string(number));
}

inline std::string MakeContigId(size_t number, size_t length, const std::string& prefix = "NODE") {
    return MakeContigId(number, prefix) + "_length_" + std::to_string(length);
}

inline std::string MakeContigId(size_t number, size_t length, double coverage, const std::string& prefix = "NODE") {
    return MakeContigId(number, length, prefix) + "_cov_" + std::to_string(coverage);
}

inline std::string MakeContigId(size_t number, size_t length, double coverage, size_t id, const std::string& prefix = "NODE") {
    return MakeContigId(number, length, coverage, prefix) + "_ID_" +  std::to_string(id);
}

inline std::string MakeRNAContigId(size_t number, size_t length, double coverage, size_t gene_id, size_t isoform_id, const std::string& prefix = "NODE") {
    return MakeContigId(number, length, coverage, prefix) + "_g" + std::to_string(gene_id)  + "_i" + std::to_string(isoform_id);
}

inline std::string AddComponentId(const std::string& s, size_t component_id) {
    return s + "_component_" + std::to_string(component_id);
}

}
