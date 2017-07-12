//
// Created by Vasily Ershov on 19.03.16.
//

#ifndef PROJECT_READ_CORRECTOR_INFO_H
#define PROJECT_READ_CORRECTOR_INFO_H

#include <functional>
#include <queue>
#include "hkmer.hpp"

namespace hammer {
namespace correction {

namespace numeric = boost::numeric::ublas;
using HRun = HomopolymerRun;

template <class Moveable>
inline Moveable pop_queue(std::priority_queue<Moveable>& queue) {
  Moveable result(std::move(const_cast<Moveable&&>(queue.top())));
  queue.pop();
  return result;
}

struct IonEvent {
  
  IonEvent(const char nucl = 0, const char observed_size = 0,
            const char fixed_size = 0, const bool to_good_correction = false)
      : nucl_(nucl),
        overserved_size_(observed_size),
        fixed_size_(fixed_size),
        is_to_good_correction_(to_good_correction) {}

  IonEvent(const IonEvent& other) = default;

  char nucl_;
  char overserved_size_;
  char fixed_size_;
  bool is_to_good_correction_;

  inline HRun FixedHRun() const {
    return HRun((uint8_t)nucl_, (uint8_t)fixed_size_);
  }

  inline HRun ObservedHRun() const {
    return HRun((uint8_t)nucl_, (uint8_t)overserved_size_);
  }
};

class CorrectedRead {
private:
  std::vector<HRun> runs_;
  std::shared_ptr<CorrectedRead> previous_;

 public:
  CorrectedRead() : previous_(nullptr) {}

  CorrectedRead(std::shared_ptr<CorrectedRead> previous)
      : previous_(previous) {}

  CorrectedRead(std::vector<HRun>&& runs,
                 std::shared_ptr<CorrectedRead> previous)
      : runs_(std::move(runs)), previous_(previous) {}

  inline void Add(const HRun hrun) { runs_.push_back(hrun); }

  size_t Size() const {
    size_t size = previous_ != nullptr ? previous_->Size() : 0;
    for (auto hrun : runs_) {
      size += hrun.len;
    }
    return size;
  }

  inline void Fill(std::string& result) const {
    if (previous_ != nullptr) {
      previous_->Fill(result);
    }

    for (auto hrun : runs_) {
      result += hrun.str();
    }
  }

  inline std::string ToString() const {
    std::string result;
    result.reserve(Size() + 10);
    Fill(result);
    return result;
  }
};

template <class PenaltyState>
class CorrectionState {
  template <class>
  friend class MoveToNextDivergence;
  template <class>
  friend class StateBuilder;

 private:
  PenaltyState penalty_state;
  HKMer kmer_;
  std::shared_ptr<CorrectedRead> current_read_ = std::shared_ptr<CorrectedRead>(nullptr);
  int16_t cursor_ = 0;
  int16_t corrections_ = 0;

 public:
  const HKMer& GetHKMer() const { return kmer_; }

  inline double Penalty() const { return penalty_state.Penalty(); }

  inline size_t TotalCorrections() const { return (size_t)corrections_; }

  const CorrectedRead* Read() const { return current_read_.get(); }

  unsigned Position() const { return (unsigned)cursor_; }
};

class CorrectionContext {
 private:
  std::vector<char> read_;
  std::vector<uint8_t> hrun_sizes_;
  const KMerData& data_;
  bool reversed_;

  inline void FillHRunSizes(const std::vector<char>& read,
                            std::vector<uint8_t>& hrun_sizes) const {
    size_t offset = 0;
    hrun_sizes.resize(read.size());

    while (offset < read.size()) {
      size_t cursor = offset;
      while (cursor < read.size() && read[cursor] == read[offset]) {
        ++cursor;
      };
      uint8_t sz = (uint8_t)(cursor - offset);
      while (sz > 0) {
        hrun_sizes[offset++] = sz;
        --sz;
      }
    }
  }

 public:
  CorrectionContext(const KMerData& data, const std::string& read,
                     bool reverse)
      : data_(data)
      , reversed_(reverse) {
    read_.resize(read.size());
    for (size_t i = 0; i < read.size(); ++i) {
      read_[i] = dignucl(read[i]);
    }

    FillHRunSizes(read_, hrun_sizes_);
  }

  inline const std::vector<char>& GetRead() const { return read_; }

  inline size_t GetOriginalOffset(const size_t offset) const {
    if (reversed_) {
      return read_.size() - offset;
    }
    return offset;
  }

  inline bool IsReversed() const { return reversed_; }

  inline HRun GetHRun(size_t offset) const {
    return HRun((uint8_t)read_[offset], (uint8_t)hrun_sizes_[offset]);
  }

  inline KMerStat const* TryGetKMerStats(const HKMer& kmer) const {
    auto idx = data_.checking_seq_idx(kmer);
    return idx == -1ULL ? nullptr : &data_[kmer];
  }

  inline bool Skip(const HKMer& kmer) const {
    auto stat = TryGetKMerStats(kmer);
    return stat != nullptr ? stat->skip() : false;
  }
};

//
template <class PenaltyCalcer>
class StateBuilder {
  using State = CorrectionState<typename PenaltyCalcer::PenaltyState>;
  const State& previous_;
  const PenaltyCalcer& penalty_calcer_;
  const CorrectionContext& context_;
  State next_;

 public:
  StateBuilder(const State& previous,
               const PenaltyCalcer& penalty_calcer,
                const CorrectionContext& context)
      : previous_(previous),
        penalty_calcer_(penalty_calcer),
        context_(context),
        next_() {
    next_.current_read_.reset(new CorrectedRead(previous_.current_read_));
    next_.kmer_ = previous_.kmer_;
    next_.penalty_state = previous_.penalty_state;
    next_.cursor_ = previous_.cursor_;
    next_.corrections_ = previous_.corrections_;
  }

  inline void AddEvent(const IonEvent& event) {
    if (event.fixed_size_ != 0) {
      const HRun run = event.FixedHRun();
      next_.kmer_ <<= run;
      next_.current_read_->Add(run);
    }

    next_.cursor_ = (int16_t)(next_.cursor_ + event.overserved_size_);
    penalty_calcer_.Update(next_.penalty_state, event,
                         context_.TryGetKMerStats(next_.kmer_));

    if (event.fixed_size_ != event.overserved_size_) {
      next_.corrections_++;
    }
  }

  inline State Build() { return next_; }

  static State Initial(const CorrectionContext& context,
                       const PenaltyCalcer& penalty,
                       unsigned skip) {
    State state;
    state.penalty_state = PenaltyCalcer::CreateState(
        context.IsReversed(), (unsigned)context.GetRead().size());
    state.current_read_.reset(new CorrectedRead());
    size_t offset = 0;
    size_t minSkip = 0;

    for (unsigned i = 0; i < hammer::K; ++i) {
      minSkip += context.GetHRun(minSkip).len;
      if (minSkip >= context.GetRead().size()) {
        break;
      }
    }

    if (minSkip > skip) {
      skip = (unsigned)minSkip;
    }
    state.cursor_ = (int16_t)skip;

    while (offset < skip) {
      HRun run = context.GetHRun(offset);
      state.kmer_ <<= run;
      state.current_read_->Add(context.GetHRun(offset));
      penalty.UpdateInitial(state.penalty_state,
                            IonEvent(run.nucl, run.len, run.len, true),
                            context.TryGetKMerStats(state.kmer_));
      offset += run.len;
    }
    return state;
  }
};

template <class PenaltyCalcer>
class MoveToNextDivergence {
  using State = CorrectionState<typename PenaltyCalcer::PenaltyState>;
  std::vector<IonEvent> Proceeded;
  State& state_;
  const CorrectionContext& context_;
  const PenaltyCalcer& calcer_;
  unsigned cursor_;

 public:
  MoveToNextDivergence(State& state,
                       const CorrectionContext& context,
                       const PenaltyCalcer& calcer)
      : state_(state),
        context_(context),
        calcer_(calcer),
        cursor_((unsigned)state.cursor_) {}

  inline bool FindNextDivergence() {
    const auto& context = context_;
    const size_t readSize = context.GetRead().size();
    HKMer currentHKMer = state_.kmer_;

    while (cursor_ < readSize) {
      const HRun hrun = context.GetHRun(cursor_);
      currentHKMer <<= hrun;

      if (calcer_.Skip(currentHKMer)) {
        Proceeded.push_back({hrun.Nucl(), hrun.Len(), hrun.Len(), true});
        cursor_ += hrun.len;
      } else {
        break;
      }
    }
    return cursor_ != (unsigned)state_.cursor_;
  }

  // we'll use it only while we move in branch…
  inline void Move() {
    for (unsigned i = 0; i < Proceeded.size(); ++i) {
      state_.current_read_->Add(Proceeded[i].FixedHRun());
      state_.kmer_ <<= Proceeded[i].FixedHRun();
      calcer_.Update(state_.penalty_state, Proceeded[i],
                    context_.TryGetKMerStats(state_.kmer_));
    }
    state_.cursor_ = (int16_t)cursor_;
  }
};

template <class PenaltyCalcer>
class SkipMayBeBadHRun {
private:
  using TState = CorrectionState<typename PenaltyCalcer::PenaltyState>;
  const TState& previous_;
  const CorrectionContext& context_;
  const PenaltyCalcer& calcer_;

 public:
  SkipMayBeBadHRun(const TState& previous,
                   const CorrectionContext& context,
                    const PenaltyCalcer& calcer)
      : previous_(previous)
        , context_(context)
        , calcer_(calcer) {}

  inline TState State() {
    StateBuilder<PenaltyCalcer> nextBuilder(previous_, calcer_, context_);
    const auto hrun = context_.GetHRun(previous_.Position());
    nextBuilder.AddEvent(IonEvent(hrun.nucl, hrun.len, hrun.len, false));
    return nextBuilder.Build();
  }
};

class HRunSizeSearcher {
 private:
  HKMer hkmer_;
  const uint8_t observed_nucl_;
  const char observed_size_;
  const std::function<bool(const hammer::HKMer&)>& is_good_func;

 public:
  HRunSizeSearcher(const HKMer& prev,
                    HRun run,
                    std::function<bool(const hammer::HKMer&)>& good)
      : hkmer_(prev),
        observed_nucl_(run.nucl),
        observed_size_(run.len),
        is_good_func(good) {
    assert(hkmer_[K - 1].nucl != run.nucl);
    hkmer_ <<= run;
  }

  inline IonEvent WithoutCorrection() {
    hkmer_[K - 1].len = observed_size_ & 0x3F;
    return IonEvent(observed_nucl_, observed_size_, observed_size_, is_good_func(hkmer_));
  }

  inline std::vector<IonEvent> TryFindInsertions(char max_error_size = 3,
                                                  const bool greedy = true) {
    std::vector<IonEvent> results;
    results.reserve(max_error_size);

    const char nucl = hkmer_[K - 1].nucl;
    for (char i = 1; i <= max_error_size; ++i) {
      hkmer_[K - 1].len = (observed_size_ + i) & 0x3F;
      if (is_good_func(hkmer_)) {
        results.push_back(
            IonEvent(nucl, observed_size_, (uint8_t)(observed_size_ + i), true));
        if (greedy) {
          break;
        }
      }
    }
    return results;
  }

  inline std::vector<IonEvent> TryFindAllDeletions(const char max_error_size = 3,
                                                    const bool greedy = true) {
    std::vector<IonEvent> results;
    results.reserve(max_error_size);

    const char nucl = hkmer_[K - 1].nucl;

    const char start = (const char)std::max(1, observed_size_ - max_error_size);

    for (char i = (char)(observed_size_ - 1); i >= start; --i) {
      hkmer_[K - 1].len = i & 0x3F;
      if (is_good_func(hkmer_)) {
        results.push_back(IonEvent(nucl, observed_size_, i, true));
        if (greedy) {
          break;
        }
      }
    }
    return results;
  }

  inline IonEvent TryFindInsertion(char max_error_size = 3) {
    const char nucl = hkmer_[K - 1].nucl;
    bool found = false;
    for (char i = 1; i <= max_error_size; ++i) {
      hkmer_[K - 1].len = (observed_size_ + i) & 0x3F;
      if (is_good_func(hkmer_)) {
        found = true;
        break;
      }
    }
    return IonEvent(nucl, observed_size_,
                     (const char)(found ? hkmer_[K - 1].len : observed_size_ + 1),
                     found);
  }

  inline IonEvent TryFindDeletion(const char max_error_size = 3) {
    const char nucl = hkmer_[K - 1].nucl;
    bool found = false;

    const char start = (const char)std::max(1, observed_size_ - max_error_size);

    for (char i = (char)(observed_size_ - 1); i >= start; --i) {
      hkmer_[K - 1].len = i & 0x3F;
      if (is_good_func(hkmer_)) {
        found = true;
        break;
      }
    }
    return IonEvent(nucl, observed_size_,
                     (const char)(found ? hkmer_[K - 1].len : observed_size_ - 1),
                     found);
  }

  inline std::vector<IonEvent> Find(const char max_error_size = 3) {
    std::vector<IonEvent> events;

    IonEvent without = WithoutCorrection();
    if (without.is_to_good_correction_) {
      events.push_back(without);
      return events;
    }

    IonEvent insertion = TryFindInsertion(max_error_size);
    if (insertion.is_to_good_correction_) {
      events.push_back(insertion);
    }

    IonEvent deletion = TryFindDeletion(max_error_size);
    if (deletion.is_to_good_correction_) {
      events.push_back(deletion);
    }

    return events;
  }
};

template <class PenaltyCalcer>
class CorrectLastHRun {
  using TState = CorrectionState<typename PenaltyCalcer::PenaltyState>;
  const TState& previous_;
  const CorrectionContext& context_;
  const PenaltyCalcer& calcer_;
  std::function<bool(const hammer::HKMer&)> is_good_function_;

  const unsigned kMaxFulldel = cfg::get().max_full_del;
  const unsigned kMaxInDel = cfg::get().max_indel;
  const unsigned kMaxFromZeroInsertion = cfg::get().max_from_zero_insertion;
  const unsigned kMaxSecondIndel = cfg::get().max_second_indel;

 private:
  inline bool AddAnotherNuclInsertions(const HRun run,
                                       const TState& previous,
                                       std::priority_queue<TState>& corrections) {
    bool found = false;
    const auto& kmer = previous.GetHKMer();

    for (uint8_t c = 0; c < 4; ++c) {
      if (c == run.nucl || c == kmer[K - 1].nucl) {
        continue;
      }

      HKMer another_nucl_insertion = kmer;
      another_nucl_insertion <<= HRun(c, 1);

      for (unsigned i = 1; i <= kMaxFromZeroInsertion; ++i) {
        another_nucl_insertion[K - 1].len = i & 0x3F;
        if (is_good_function_(another_nucl_insertion)) {
          HRunSizeSearcher rest_searcher(another_nucl_insertion, run, is_good_function_);
          auto events = rest_searcher.Find((const char)kMaxSecondIndel);
          for (auto& event : events) {
            if (event.is_to_good_correction_) {
              StateBuilder<PenaltyCalcer> builder(previous, calcer_, context_);
              builder.AddEvent(IonEvent(c, 0, (const char)i, true));  // new insertion
              builder.AddEvent(event);
              corrections.emplace(builder.Build());
              found = true;
            }
          }
          break;
        }
      }
    }
    return found;
  }

 public:
  CorrectLastHRun(const TState& previous,
                  const CorrectionContext& context,
                   const PenaltyCalcer& calcer)
      : previous_(previous),
        context_(context),
        calcer_(calcer),
        is_good_function_(calcer_.Good()) {}

  inline void AddOnlySimpleCorrections(std::priority_queue<TState>& corrections,
                                       unsigned indel_size = 1) {
    const unsigned cursor = previous_.Position();
    const HRun run = context_.GetHRun(cursor);

    if (!is_good_function_(previous_.GetHKMer())) {
      return;
    }

    {
      HRunSizeSearcher searcher(previous_.GetHKMer(), run, is_good_function_);
      auto insertion = searcher.TryFindInsertion((char)indel_size);
      if (insertion.is_to_good_correction_) {
        StateBuilder<PenaltyCalcer> builder(previous_, calcer_, context_);
        builder.AddEvent(insertion);
        corrections.emplace(builder.Build());
      }

      auto deletion = searcher.TryFindDeletion((const char)indel_size);
      if (deletion.is_to_good_correction_) {
        StateBuilder<PenaltyCalcer> builder(previous_, calcer_, context_);
        builder.AddEvent(deletion);
        corrections.emplace(builder.Build());
      }
    }
    //
    if (run.len == 1 && (cursor + 1 < context_.GetRead().size())) {
      auto nextRun = context_.GetHRun(cursor + 1);
      {
        for (char c = 0; c < 4; ++c) {
          if (c == run.nucl || c == nextRun.nucl) {
            continue;
          }

          HKMer kmer = previous_.GetHKMer();
          kmer <<= HRun((uint8_t)c, 1);

          if (is_good_function_(kmer)) {
            StateBuilder<PenaltyCalcer> builder(previous_, calcer_, context_);
            builder.AddEvent(IonEvent(run.nucl, run.len, 0, true));
            builder.AddEvent(IonEvent(c, 0, 1, true));
            corrections.emplace(builder.Build());
          }
        }
      }
    } else if (run.len > 2) {
      for (char c = 0; c < 4; ++c) {
        if (c == run.nucl) {
          continue;
        }

        HKMer kmer = previous_.GetHKMer();
        kmer <<= HRun(run.nucl, (uint8_t)(run.len - 1));
        kmer <<= HRun(c, 1);
        kmer <<= HRun(run.nucl, 1);

        const unsigned maxLen = (unsigned)(run.len - 2);
        for (unsigned i = 0; i < maxLen; ++i) {
          kmer[K - 3].len = (i + 1) & 0x3F;
          kmer[K - 1].len = (maxLen - i) & 0x3F;
          if (is_good_function_(kmer)) {
            StateBuilder<PenaltyCalcer> builder(previous_, calcer_, context_);
            builder.AddEvent(IonEvent(run.nucl, (char)(i + 2), (char)(i + 1), true));
            builder.AddEvent(IonEvent(c, (char)0, (char)1, true));
            builder.AddEvent(
                IonEvent(run.nucl, (char)(maxLen - i), (char)(maxLen - i), true));
            corrections.emplace(builder.Build());
          }
        }
      }
    }
  }

  inline bool AddPossibleCorrections(std::priority_queue<TState>& corrections) {
    const unsigned cursor = previous_.Position();
    const HRun run = context_.GetHRun(cursor);
    bool found = false;

    if (is_good_function_(previous_.GetHKMer())) {
      HRunSizeSearcher searcher(previous_.GetHKMer(), run, is_good_function_);
      {
        auto insertions = searcher.TryFindInsertions((char)kMaxInDel);
        for (const auto& insertion : insertions) {
          if (insertion.is_to_good_correction_) {
            StateBuilder<PenaltyCalcer> builder(previous_, calcer_, context_);
            builder.AddEvent(insertion);
            corrections.emplace(builder.Build());
            found = true;
          }
        }
      }

      {
        auto deletions = searcher.TryFindAllDeletions((const char)std::max((int)run.len, 1));
        if (deletions.size()) {
          for (const auto& deletion : deletions) {
            const uint8_t restSize =
                (uint8_t)(deletion.overserved_size_ - deletion.fixed_size_);
            if (restSize <= kMaxInDel) {
              StateBuilder<PenaltyCalcer> builder(previous_, calcer_, context_);
              builder.AddEvent(deletion);
              corrections.emplace(builder.Build());
            }

            // Try insertion after part of hrun. Errors of type aaaaa -> aaa g
            // aa
            if (restSize > 1) {
              StateBuilder<PenaltyCalcer> indel_builder(previous_,
                                                        calcer_,
                                                        context_);
              const IonEvent partDel = IonEvent(
                  deletion.nucl_, deletion.fixed_size_, deletion.fixed_size_, true);
              indel_builder.AddEvent(partDel);
              const TState state = indel_builder.Build();
              found |= AddAnotherNuclInsertions(HRun(deletion.nucl_, restSize),
                                                state, corrections);
            }
          }
          found = true;
        }
      }

      if (!found) {
        found |= AddAnotherNuclInsertions(run, previous_, corrections);

        int read_size = (int)context_.GetRead().size();
        const int next_cursor = cursor + run.len;

        if (next_cursor >= read_size) {
          return found;
        }
        const HRun next_run = context_.GetHRun((size_t)next_cursor);

        // try full deletion of hrun.
        if (run.len <= kMaxFulldel) {
          if (next_run.nucl != previous_.GetHKMer()[K - 1].nucl) {
            StateBuilder<PenaltyCalcer> builder(previous_, calcer_, context_);
            builder.AddEvent(IonEvent(run.nucl, run.len, 0, true));  // full deletion
            corrections.emplace(builder.Build());
            found = true;
          } else {
            {
              StateBuilder<PenaltyCalcer> builder(previous_, calcer_, context_);
              builder.AddEvent(IonEvent(run.nucl, run.len, 0, true));  // full deletion
              builder.AddEvent(IonEvent(next_run.nucl, next_run.len, 0, true));  // full deletion
              corrections.emplace(builder.Build());
            }
            {
              StateBuilder<PenaltyCalcer> builder(previous_, calcer_, context_);
              builder.AddEvent(IonEvent(run.nucl, run.len, 0, true));  // full deletion
              auto state = builder.Build();
              found |= AddAnotherNuclInsertions(next_run, state, corrections);
            }
          }
        }
      }
    } else {
      {
        HKMer test = previous_.GetHKMer();
        HRun fixed = run;
        fixed.len = (fixed.len + 1) & 0x3F;
        test <<= fixed;
        size_t local_cursor = cursor + run.len;

        for (unsigned i = 0; i < (K - 1); ++i) {
          if (local_cursor >= context_.GetRead().size()) {
            break;
          }
          const HRun cursorRun = context_.GetHRun(local_cursor);
          test <<= cursorRun;
          local_cursor += cursorRun.len;

          if (is_good_function_(test)) {
            found = true;
            StateBuilder<PenaltyCalcer> builder(previous_, calcer_, context_);
            builder.AddEvent(IonEvent(run.nucl, run.len, (char) (run.len + 1), false));
            corrections.emplace(builder.Build());
            break;
          }
        }
      }

      if (run.len > 1) {
        HKMer test = previous_.GetHKMer();
        HRun fixed = run;
        fixed.len = (fixed.len - 1) & 0x3F;
        test <<= fixed;

        size_t local_cursor = cursor + run.len;

        for (unsigned i = 0; i < (K - 1); ++i) {
          if (local_cursor >= context_.GetRead().size()) {
            break;
          }
          const HRun cursorRun = context_.GetHRun(local_cursor);
          test <<= cursorRun;
          local_cursor += cursorRun.len;

          if (is_good_function_(test)) {
            found = true;
            StateBuilder<PenaltyCalcer> builder(previous_, calcer_, context_);
            builder.AddEvent(IonEvent(run.nucl, run.len, (uint8_t)(run.len - 1), false));
            corrections.emplace(builder.Build());
            break;
          }
        }
      }
    }
    return found;
  }
};
}  // namespace correction
}  // namespace hammer

namespace std {
using namespace hammer::correction;

template <class PenaltyState>
struct less<CorrectionState<PenaltyState> > {
  bool operator()(const CorrectionState<PenaltyState>& left,
                  const CorrectionState<PenaltyState>& right) const {
    return left.Penalty() < right.Penalty() ||
           (left.Penalty() == right.Penalty() &&
            left.Position() < right.Position());
  }
};

}  // namespace std

#endif  // PROJECT_READ_CORRECTOR_INFO_H
