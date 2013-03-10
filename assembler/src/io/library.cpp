#include "io/library.hpp"

#include <yaml-cpp/yaml.h>

#include <string>
#include <iostream>

using namespace io;

namespace YAML {
template<>
struct convert<DataSet> {
  static bool decode(const Node& node, DataSet& rhs) {
    if (!node.IsSequence())
      return false;

    rhs.load(node);
    return true;
  }
};

template<>
struct convert<SequencingLibrary> {
  static bool decode(const Node& node, SequencingLibrary& rhs) {
    rhs.load(node);
    return true;
  }
};

template<>
struct convert<LibraryOrientation> {
  static bool decode(const Node& node, LibraryOrientation& rhs) {
    std::string orientation = node.as<std::string>();

    if (orientation == "fr")
      rhs = LibraryOrientation::FR;
    else if (orientation == "rf")
      rhs = LibraryOrientation::RF;
    else if (orientation == "ff")
      rhs = LibraryOrientation::FF;
    else if (orientation == "rr")
      rhs = LibraryOrientation::RR;
    else
      return false;

    return true;
  }
};

template<>
struct convert<LibraryType> {
  static bool decode(const Node& node, LibraryType& rhs) {
    std::string type = node.as<std::string>();

    if (type == "paired-end")
      rhs = LibraryType::PairedEnd;
    else if (type == "single")
      rhs = LibraryType::SingleReads;
    else
      return false;

    return true;
  }

};
}

void DataSet::load(const std::string &filename) {
  YAML::Node config = YAML::LoadFile(filename);

  *this = config.as<DataSet>();
}

void SequencingLibrary::load(const YAML::Node &node) {
  orientation_ = node["orientation"].as<io::LibraryOrientation>();
  type_ = node["type"].as<LibraryType>();

  switch (type_) {
    case LibraryType::PairedEnd:
      left_paired_reads_ = node["left reads"].as<std::vector<std::string> >();
      right_paired_reads_ = node["right reads"].as<std::vector<std::string> >();

      if (left_paired_reads_.size() != right_paired_reads_.size())
        throw("Left and right reads lists should have equal length");

      // FALLTHROUGH in case of single reads present!
      if (!node["single reads"])
        break;
    case LibraryType::SingleReads:
      single_reads_ = node["single reads"].as<std::vector<std::string> >();
      break;
    default:
      // Impossible
      std::cerr << node << std::endl;
      throw("Unsupported library type");
  }
}

void DataSet::load(const YAML::Node &node) {
  for (YAML::const_iterator it = node.begin(); it != node.end(); ++it)
    libraries_.push_back(it->as<SequencingLibrary>());
}
