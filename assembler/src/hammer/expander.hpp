#ifndef __HAMMER_EXPANDER_HPP__
#define __HAMMER_EXPANDER_HPP__

class KMerData;
class Read;

#include <cstring>

class Expander {
  KMerData &data_;
  size_t changed_;
  
 public:
  Expander(KMerData &data) 
      : data_(data), changed_(0) {}

  size_t changed() const { return changed_; }

  bool operator()(const Read &r);
};

#endif
