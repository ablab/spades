#ifndef HAMMER_VALIDKMERGENERATOR_HPP_
#define HAMMER_VALIDKMERGENERATOR_HPP_
template<uint32_t kK>
class ValidKMerGenerator {
 public:
  explicit ValidKMerGenerator(const Read &read, int bad_quality_threshold = 2) :
      bad_quality_threshold_(bad_quality_threshold),
      pos_(-1),
      end_(-1),
      has_more_(true),
      first(true),
      kmer_(),
      seq_(read.getSequenceString()),
      qual_(read.getQualityString()){
    TrimBadQuality();
    Next();
  }
  bool HasMore() const {
    return has_more_;
  }
  const Seq<kK>& kmer() {
    return kmer_;
  }
  void Next();
 private:
  void TrimBadQuality();
  uint32_t bad_quality_threshold_;
  uint32_t pos_;
  uint32_t end_;
  bool has_more_;
  bool first;
  Seq<kK> kmer_;
  const std::string &seq_;
  const std::string &qual_;
};

template<uint32_t kK>
void ValidKMerGenerator<kK>::TrimBadQuality() {
  pos_ = 0;
  for (; pos_ < qual_.size(); ++pos_) {
    if ((uint32_t)qual_[pos_] > bad_quality_threshold_)
      break;
  }
  end_ = qual_.size();
  for (; end_ > pos_; --end_) {
    if ((uint32_t)qual_[end_ - 1] > bad_quality_threshold_)
      break;
  }
}

template<uint32_t kK>
void ValidKMerGenerator<kK>::Next() {
  if (pos_ + kK > end_) {
    has_more_ = false;
  } else if (first || !is_nucl(seq_[pos_ + kK - 1])) {
    // in this case we have to look for new k-mer
    uint32_t start_hypothesis = pos_;
    uint32_t i = pos_;
    for (; i < seq_.size(); ++i) {
      if (i == kK + start_hypothesis) {
        break;
      }
      if (!is_nucl(seq_[i])) {
        start_hypothesis = i + 1;
      }
    }
    if (i == kK + start_hypothesis) {
      kmer_ = Seq<kK>(seq_.data() + start_hypothesis, false);
      pos_ = start_hypothesis + 1;
    } else {
      has_more_ = false;
    }
  } else {
    // good case we can just cyclic shift our answe
    kmer_ = kmer_ << seq_[pos_ + kK - 1];
    ++pos_;
  }
  first = false;
}
#endif //  HAMMER_VALIDKMERGENERATOR_HPP__
