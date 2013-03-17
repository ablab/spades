#ifndef __IO_LIBRARY_HPP__
#define __IO_LIBRARY_HPP__

#include "adt/chained_iterator.hpp"

#include <boost/iterator/iterator_facade.hpp>

#include <string>
#include <vector>
#include <utility>
#include <iostream>

namespace YAML {
class Node;
};

namespace io {

class DataSetReader;

enum class LibraryType {
  SingleReads,
  PairedEnd,
  MatePairs,
  LongSingleReads
};

enum class LibraryOrientation {
  FR,
  FF,
  RF,
  RR,
  Undefined
};

class SequencingLibrary {
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

  SequencingLibrary()
      : type_(LibraryType::PairedEnd), orientation_(LibraryOrientation::FR) {}

  void load(const YAML::Node &node);

  LibraryType type() const { return type_; }
  LibraryOrientation orientation() const { return orientation_; }

  void clear() {
    left_paired_reads_.clear();
    right_paired_reads_.clear();
    single_reads_.clear();
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


 private:
  LibraryType type_;
  LibraryOrientation orientation_;

  std::vector<std::string> left_paired_reads_;
  std::vector<std::string> right_paired_reads_;
  std::vector<std::string> single_reads_;

  friend class DataSet;
};

// Just convenient wrapper to "unwrap" the iterators over datasets.
class DataSet {
  typedef std::vector<SequencingLibrary> LibraryStorage;

 public:
  typedef LibraryStorage::iterator iterator;
  typedef LibraryStorage::const_iterator const_iterator;
  typedef chained_iterator<SequencingLibrary::single_reads_iterator> single_reads_iterator;
  typedef chained_iterator<SequencingLibrary::paired_reads_iterator> paired_reads_iterator;

  DataSet() {}
  explicit DataSet(const std::string &path) { load(path); }
  DataSet(const YAML::Node &node) { load(node); }

  void load(const std::string &);
  void load(const YAML::Node &node);
  void print() const;

  void clear() { libraries_.clear(); }
  void push_back(const SequencingLibrary &lib) {
    libraries_.push_back(lib);
  }
  SequencingLibrary& operator[](size_t n) { return libraries_[n]; }
  const SequencingLibrary& operator[](size_t n) const { return libraries_[n]; }
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


};

#endif // __IO_LIBRARY_HPP__
