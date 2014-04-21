#ifndef __IO_LIBRARY_HPP__
#define __IO_LIBRARY_HPP__

#include "adt/chained_iterator.hpp"

#include <boost/iterator/iterator_facade.hpp>
#include <yaml-cpp/yaml.h>

#include <string>
#include <vector>
#include <utility>
#include <iostream>
#include <fstream>
#include <functional>
#include <algorithm>
#include "path_helper.hpp"

namespace io {

class DataSetReader;

enum class LibraryType {
  SingleReads,
  PairedEnd,
  MatePairs,
  HQMatePairs,
  PacBioReads,
  SangerReads,
  TrustedContigs,
  UntrustedContigs,
};

enum class LibraryOrientation {
  FR,
  FF,
  RF,
  RR,
  Undefined
};

struct update_relative_filename : public std::binary_function<std::string, std::string, std::string> {
  std::string operator() (const std::string &filename, const std::string &input_dir) const {
    if (filename[0] == '/')
      return filename;
    return input_dir + filename;
  }
};

class SequencingLibraryBase {
 public:
  class paired_reads_iterator :
      public boost::iterator_facade<paired_reads_iterator,
                                    std::pair<std::string, std::string>,
                                    boost::forward_traversal_tag,
                                    std::pair<std::string, std::string> > {

    typedef std::vector<std::string>::const_iterator inner_iterator;

   public:
    paired_reads_iterator(inner_iterator left, inner_iterator right)
        : left_(left), right_(right){}

   private:
    friend class boost::iterator_core_access;

    void increment() { ++left_; ++right_; }
    bool equal(const paired_reads_iterator &other) const {
      return this->left_ == other.left_ && this->right_ == other.right_;
    }
    std::pair<std::string, std::string> dereference() const {
      return std::make_pair(*left_, *right_);
    }

    inner_iterator left_;
    inner_iterator right_;
  };

  typedef chained_iterator<std::vector<std::string>::const_iterator> single_reads_iterator;

  SequencingLibraryBase()
      : type_(LibraryType::PairedEnd), orientation_(LibraryOrientation::FR) {}

  void load(const YAML::Node &node);

  LibraryType type() const { return type_; }
  void set_type(LibraryType type) { type_ = type; }
  LibraryOrientation orientation() const { return orientation_; }
  void set_orientation(LibraryOrientation orientation) { orientation_ = orientation; }

  void clear() {
    left_paired_reads_.clear();
    right_paired_reads_.clear();
    single_reads_.clear();
  }

  void update_relative_reads_filenames(const std::string &input_dir) {
      std::transform(left_paired_reads_.begin(), left_paired_reads_.end(), left_paired_reads_.begin(), std::bind2nd(update_relative_filename(), input_dir));
      std::transform(right_paired_reads_.begin(), right_paired_reads_.end(), right_paired_reads_.begin(), std::bind2nd(update_relative_filename(), input_dir));
      std::transform(single_reads_.begin(), single_reads_.end(), single_reads_.begin(), std::bind2nd(update_relative_filename(), input_dir));
  }

  void push_back_single(const std::string &reads) {
    single_reads_.push_back(reads);
  }

  void push_back_paired(const std::string &left, const std::string &right) {
    left_paired_reads_.push_back(left);
    right_paired_reads_.push_back(right);
  }

  paired_reads_iterator paired_begin() const {
    return paired_reads_iterator(left_paired_reads_.begin(), right_paired_reads_.begin());
  }
  paired_reads_iterator paired_end() const {
    return paired_reads_iterator(left_paired_reads_.end(), right_paired_reads_.end());
  }

  single_reads_iterator reads_begin() const {
    // NOTE: We have a contract with single_end here. Single reads always go last!
    single_reads_iterator res(left_paired_reads_.begin(), left_paired_reads_.end());
    res.join(right_paired_reads_.begin(), right_paired_reads_.end());
    res.join(single_reads_.begin(), single_reads_.end());

    return res;
  }
  single_reads_iterator reads_end() const {
    // NOTE: Do not forget about the contract with single_begin here!
    return single_reads_iterator(single_reads_.end(), single_reads_.end());
  }

  single_reads_iterator single_begin() const {
    return single_reads_iterator(single_reads_.begin(), single_reads_.end());
  }
  single_reads_iterator single_end() const {
    // NOTE: Do not forget about the contract with single_begin here!
    return single_reads_iterator(single_reads_.end(), single_reads_.end());
  }

  bool is_graph_contructable() const {
    return (type_ == io::LibraryType::PairedEnd ||
            type_ == io::LibraryType::SingleReads ||
            type_ == io::LibraryType::HQMatePairs);
  }

  bool is_paired() const {
    return (type_ == io::LibraryType::PairedEnd ||
            type_ == io::LibraryType::MatePairs||
            type_ == io::LibraryType::HQMatePairs);
  }


  bool is_repeat_resolvable() const {
    return (type_ == io::LibraryType::PairedEnd ||
            type_ == io::LibraryType::HQMatePairs ||
            type_ == io::LibraryType::MatePairs ||
            type_ == io::LibraryType::PacBioReads ||
            type_ == io::LibraryType::SangerReads ||
            type_ == io::LibraryType::TrustedContigs ||
            type_ == io::LibraryType::UntrustedContigs);
  }

  bool is_pacbio_alignable() const {
    return (type_ == io::LibraryType::PacBioReads ||
            type_ == io::LibraryType::SangerReads ||
            type_ == io::LibraryType::TrustedContigs ||
            type_ == io::LibraryType::UntrustedContigs);
  }

 private:
  LibraryType type_;
  LibraryOrientation orientation_;

  std::vector<std::string> left_paired_reads_;
  std::vector<std::string> right_paired_reads_;
  std::vector<std::string> single_reads_;
};

struct NoData {};

template<class Data = NoData>
class SequencingLibrary: public SequencingLibraryBase {
 public:
  const Data& data() const {
    return data_;
  }
  Data& data() {
    return data_;
  }

 private:
  Data data_;
};

// Just convenient wrapper to "unwrap" the iterators over libraries.
template<class Data = NoData>
class DataSet {
  typedef SequencingLibrary<Data> Library;
  typedef std::vector<Library> LibraryStorage;

 public:
  typedef typename LibraryStorage::iterator iterator;
  typedef typename LibraryStorage::const_iterator const_iterator;
  typedef chained_iterator<typename Library::single_reads_iterator> single_reads_iterator;
  typedef chained_iterator<typename Library::paired_reads_iterator> paired_reads_iterator;

  DataSet() {}
  explicit DataSet(const std::string &path) { load(path); }
  DataSet(const YAML::Node &node) { load(node); }

  void load(const std::string &filename) {
    YAML::Node config = YAML::LoadFile(filename);

    std::string input_dir = path::parent_path(filename);
    if (input_dir[input_dir.length() - 1] != '/')
        input_dir += '/';

    load(config);

    for (auto it = libraries_.begin(); it != libraries_.end(); ++it)
      it->update_relative_reads_filenames(input_dir);
  }

  void save(const std::string &filename) const {
    std::ofstream ofs(filename.c_str());
    ofs << YAML::Node(*this);
  }

  void load(const YAML::Node &node) {
    clear();
    for (YAML::const_iterator it = node.begin(); it != node.end(); ++it)
      libraries_.push_back(it->as<Library>());
  }

  void clear() { libraries_.clear(); }
  void push_back(const Library &lib) {
    libraries_.push_back(lib);
  }
  Library& operator[](size_t n) { return libraries_[n]; }
  const Library& operator[](size_t n) const { return libraries_[n]; }
  size_t lib_count() const { return libraries_.size(); }
  iterator library_begin() { return libraries_.begin(); }
  const_iterator library_begin() const { return libraries_.cbegin(); }
  const_iterator library_cbegin() const { return libraries_.cbegin(); }
  iterator library_end() { return libraries_.end(); }
  const_iterator library_end() const { return libraries_.cend(); }
  const_iterator library_cend() const { return libraries_.cend(); }

  single_reads_iterator reads_begin() const {
    auto it = libraries_.begin();
    single_reads_iterator res(it->reads_begin(), it->reads_end());
    ++it;
    for (auto end = libraries_.end(); it != end; ++it)
      res.join(it->reads_begin(), it->reads_end());

    return res;
  }

  single_reads_iterator reads_end() const {
    return single_reads_iterator(libraries_.back().reads_end(), libraries_.back().reads_end());
  }

  single_reads_iterator single_begin() const {
    auto it = libraries_.begin();
    single_reads_iterator res(it->single_begin(), it->single_end());
    ++it;
    for (auto end = libraries_.end(); it != end; ++it)
      res.join(it->single_begin(), it->single_end());

    return res;
  }

  single_reads_iterator single_end() const {
    return single_reads_iterator(libraries_.back().single_end(), libraries_.back().single_end());
  }

  paired_reads_iterator paired_begin() const {
    auto it = libraries_.begin();
    paired_reads_iterator res(it->paired_begin(), it->paired_end());
    ++it;
    for (auto end = libraries_.end(); it != end; ++it)
      res.join(it->paired_begin(), it->paired_end());

    return res;
  }

  paired_reads_iterator paired_end() const {
    return paired_reads_iterator(libraries_.back().paired_end(), libraries_.back().paired_end());
  }

 private:
  LibraryStorage libraries_;
};

}

namespace YAML {
template<>
struct convert<io::SequencingLibraryBase > {
  static Node encode(const io::SequencingLibraryBase& rhs);
  static bool decode(const Node& node, io::SequencingLibraryBase& rhs);
};

template<>
struct convert<io::SequencingLibrary<> > {
  static Node encode(const io::SequencingLibrary<>& rhs);
  static bool decode(const Node& node, io::SequencingLibrary<>& rhs);
};

template<class Data>
struct convert<io::DataSet<Data> > {
  static Node encode(const io::DataSet<Data>& rhs) {
    Node node;

    for (auto it = rhs.library_begin(), et = rhs.library_end(); it != et; ++it)
      node.push_back(*it);

    return node;
  }
};
}


#endif // __IO_LIBRARY_HPP__
