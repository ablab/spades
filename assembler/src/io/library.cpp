//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "io/library.hpp"

#include <yaml-cpp/yaml.h>

#include <string>
#include <iostream>

using namespace io;

namespace YAML {
template<>
struct convert<LibraryOrientation> {
  static Node encode(const LibraryOrientation &rhs) {
    switch (rhs) {
      case LibraryOrientation::FR:
        return Node("fr");
      case LibraryOrientation::RF:
        return Node("rf");
      case LibraryOrientation::FF:
        return Node("ff");
      case LibraryOrientation::RR:
        return Node("rr");
      case LibraryOrientation::Undefined:
      default:
        return Node("undefined");
    }
  }

  static bool decode(const Node& node, LibraryOrientation& rhs) {
    std::string orientation = node.as<std::string>("");

    if (orientation == "fr")
      rhs = LibraryOrientation::FR;
    else if (orientation == "rf")
      rhs = LibraryOrientation::RF;
    else if (orientation == "ff")
      rhs = LibraryOrientation::FF;
    else if (orientation == "rr")
      rhs = LibraryOrientation::RR;
    else
      rhs = LibraryOrientation::Undefined;

    return true;
  }
};

template<>
struct convert<LibraryType> {
  static Node encode(const LibraryType &rhs) {
    switch (rhs) {
      case LibraryType::PairedEnd:
        return Node("paired-end");
      case LibraryType::SingleReads:
        return Node("single");
      case LibraryType::MatePairs:
        return Node("mate-pairs");
      case LibraryType::HQMatePairs:
        return Node("hq-mate-pairs");
      case LibraryType::PacBioReads:
        return Node("pacbio");
      case LibraryType::SangerReads:
        return Node("sanger");
      case LibraryType::NanoporeReads:
        return Node("nanopore");
      case LibraryType::TrustedContigs:
        return Node("trusted-contigs");
      case LibraryType::UntrustedContigs:
        return Node("untrusted-contigs");
      case LibraryType::PathExtendContigs:
        return Node("path-extend-contigs");
      default:
        return Node();
    }
  }

  static bool decode(const Node& node, LibraryType& rhs) {
    std::string type = node.as<std::string>();

    if (type == "paired-end")
      rhs = LibraryType::PairedEnd;
    else if (type == "mate-pairs")
      rhs = LibraryType::MatePairs;
    else if (type == "hq-mate-pairs")
      rhs = LibraryType::HQMatePairs;
    else if (type == "pacbio")
      rhs = LibraryType::PacBioReads;
    else if (type == "single")
      rhs = LibraryType::SingleReads;
    else if (type == "sanger")
      rhs = LibraryType::SangerReads;
    else if (type == "nanopore")
      rhs = LibraryType::NanoporeReads;
    else if (type == "trusted-contigs")
      rhs = LibraryType::TrustedContigs;
    else if (type == "untrusted-contigs")
      rhs = LibraryType::UntrustedContigs;
    else if (type == "path-extend-contigs")
      rhs = LibraryType::PathExtendContigs;
    else
      return false;
    return true;
  }

};

Node convert<SequencingLibraryBase>::encode(const io::SequencingLibraryBase& rhs) {
    Node node;

    node["orientation"] = rhs.orientation();
    node["type"] = rhs.type();

    for (const auto& read_pair : rhs.paired_reads()) {
      node["left reads"].push_back(read_pair.first);
      node["right reads"].push_back(read_pair.second);
    }
    for (const auto& reads : rhs.single_reads())
      node["single reads"].push_back(reads);

    return node;
}

bool convert<SequencingLibraryBase>::decode(const Node& node, SequencingLibraryBase& rhs) {
    rhs.load(node);
    return true;
}

Node convert<io::SequencingLibrary<> >::encode(const io::SequencingLibrary<>& rhs) {
  return convert<io::SequencingLibraryBase>::encode(rhs);
}

bool convert<io::SequencingLibrary<> >::decode(const Node& node, io::SequencingLibrary<>& rhs) {
  rhs.load(node);
  return true;
}

} // namespace YAML

void SequencingLibraryBase::load(const YAML::Node &node) {
  orientation_ = node["orientation"].as<io::LibraryOrientation>(LibraryOrientation::Undefined);
  type_ = node["type"].as<LibraryType>();

  switch (type_) {
    case LibraryType::PairedEnd:
    case LibraryType::MatePairs:
    case LibraryType::HQMatePairs:
      left_paired_reads_ = node["left reads"].as<std::vector<std::string> >();
      right_paired_reads_ = node["right reads"].as<std::vector<std::string> >();

      if (left_paired_reads_.size() != right_paired_reads_.size())
        throw("Left and right reads lists should have equal length");

      if (orientation_ == LibraryOrientation::Undefined)
        throw("Orientation for paired reads should be specified");

      // FALLTHROUGH in case of single reads present!
      if (!node["single reads"])
        break;
    case LibraryType::SingleReads:
    case LibraryType::PacBioReads:
    case LibraryType::SangerReads:
    case LibraryType::NanoporeReads:
    case LibraryType::TrustedContigs:
    case LibraryType::UntrustedContigs:
    case LibraryType::PathExtendContigs:
      single_reads_ = node["single reads"].as<std::vector<std::string> >();
      break;
    default:
      // Impossible
      std::cerr << node << std::endl;
      throw("Unsupported library type");
  }
}
