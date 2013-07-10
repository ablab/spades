#ifndef __ADT_CHAINED_ITERATOR_HPP__
#define __ADT_CHAINED_ITERATOR_HPP__

#include <boost/iterator/iterator_facade.hpp>

#include <iterator>
#include <vector>

template<class It>
class chained_iterator :
      public boost::iterator_facade<chained_iterator<It>,
                                    typename std::iterator_traits<It>::value_type,
                                    boost::forward_traversal_tag,
                                    typename std::iterator_traits<It>::value_type> {
 public:
  chained_iterator(It begin, It end) :
      section_(0), current_(begin) {
    join(begin, end);
  }

  void join(It begin, It end) {
    begins_.push_back(begin);
    ends_.push_back(end);
    skip_empty();
  }

 private:
  friend class boost::iterator_core_access;

  void skip_empty() {
    while ((section_ + 1) < begins_.size() &&
           current_ == ends_[section_])
      current_ = begins_[++section_];
  }

  void increment() {
    skip_empty();
    ++current_;
    skip_empty();
  }

  bool equal(const chained_iterator &other) const {
    return current_ == other.current_;
  }

  typename std::iterator_traits<It>::value_type dereference() const {
    return *current_;
  }

  size_t section_;
  It current_;
  std::vector<It> begins_;
  std::vector<It> ends_;
};


#endif
