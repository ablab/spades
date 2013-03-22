#include "io/library.hpp"

#include <yaml-cpp/yaml.h>

#include <string>
#include <iostream>
#include <fstream>

using namespace io;

namespace YAML {
template<>
struct convert<DataSet> {
  static Node encode(const DataSet& rhs) {
    Node node;
    for (auto it = rhs.library_begin(), et = rhs.library_end(); it != et; ++it)
      node.push_back(*it);

    return node;
  }

  static bool decode(const Node& node, DataSet& rhs) {
    if (!node.IsSequence())
      return false;

    rhs.load(node);
    return true;
  }
};

template<>
struct convert<SequencingLibrary> {
  static Node encode(const SequencingLibrary& rhs) {
    Node node;

    node["orientation"] = rhs.orientation();
    node["type"] = rhs.type();
    if (rhs.insert_size())
      node["insert size"] = rhs.insert_size();

    for (auto it = rhs.paired_begin(), et = rhs.paired_end(); et != it; ++it) {
      node["left reads"].push_back(it->first);
      node["right reads"].push_back(it->second);
    }
    for (auto it = rhs.single_begin(), et = rhs.single_end(); et != it; ++it)
      node["single"].push_back(*it);

    return node;
  }

  static bool decode(const Node& node, SequencingLibrary& rhs) {
    rhs.load(node);
    return true;
  }
};

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
        return Node("undefined");
    }
  }

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
      case LibraryType::LongSingleReads:
        return Node("long-single");
    }
  }

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

void DataSet::save(const std::string &filename) const {
  std::ofstream ofs(filename.c_str());
  ofs << YAML::Node(*this);
}

void SequencingLibrary::load(const YAML::Node &node) {
  orientation_ = (node["orientation"] ?
                  node["orientation"].as<io::LibraryOrientation>() : LibraryOrientation::Undefined);
  type_ = node["type"].as<LibraryType>();

  switch (type_) {
    case LibraryType::PairedEnd:
      left_paired_reads_ = node["left reads"].as<std::vector<std::string> >();
      right_paired_reads_ = node["right reads"].as<std::vector<std::string> >();

      if (left_paired_reads_.size() != right_paired_reads_.size())
        throw("Left and right reads lists should have equal length");

      if (orientation_ == LibraryOrientation::Undefined)
        throw("Orientation for paired reads should be specified");

      if (node["insert size"])
        insert_size_ = node["insert size"].as<unsigned>();

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
