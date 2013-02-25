#pragma once

#include "io/ireader.hpp"
#include "io/single_read.hpp"

namespace io {

template <class Read>
class SequenceReader : public IReader<Read> {
 public:
  explicit SequenceReader(const Sequence &sequence, const std::string &name = "")
      : sequence_(sequence),
        name_(name),
        opened_(true),
        eof_(false) {
  }

  virtual ~SequenceReader() {
  }

  virtual bool is_open() {
    return opened_;
  }

  virtual bool eof() {
    return eof_;
  }

  virtual void close() {
    opened_ = false;
  }

  void reset() {
    eof_ = false;
    opened_ = true;
  }

  ReadStat get_stat() const {
        return ReadStat();
  }

  SequenceReader<Read> &operator>>(Read &read);

 private:
  Sequence sequence_;
  std::string name_;
  bool opened_;
  bool eof_;
};

template <>
SequenceReader<SingleRead> &SequenceReader<SingleRead>::operator>>(SingleRead &read) {
  if (!eof_) {
    read = SingleRead(name_, sequence_.str());
    eof_ = true;
  }
  return *this;
}

template <>
SequenceReader<SingleReadSeq> &SequenceReader<SingleReadSeq>::operator>>(SingleReadSeq &read) {
  if (!eof_) {
    read = SingleReadSeq(sequence_);
    eof_ = true;
  }
  return *this;
}

}
