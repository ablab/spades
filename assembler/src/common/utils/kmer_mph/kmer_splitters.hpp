//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "kmer_splitter.hpp"
#include "io/reads/io_helper.hpp"
#include "adt/iterator_range.hpp"

namespace utils {

using RtSeqKMerSplitter = kmers::KMerSortingSplitter<RtSeq>;

template<class KmerFilter>
class DeBruijnKMerSplitter : public RtSeqKMerSplitter {
 private:
  KmerFilter kmer_filter_;
 protected:
  size_t read_buffer_size_;
 protected:
  bool FillBufferFromSequence(const Sequence &seq,
                              unsigned thread_id) {
      if (seq.size() < this->K_)
        return false;

      RtSeq kmer = seq.start<RtSeq>(this->K_) >> 'A';
      bool stop = false;
      for (size_t j = this->K_ - 1; j < seq.size(); ++j) {
        kmer <<= seq[j];
        if (!kmer_filter_.filter(kmer))
          continue;

        stop |= this->push_back_internal(kmer, thread_id);
      }

      return stop;
  }

    bool FillBufferFromSequence(const RtSeq &seq,
                                unsigned thread_id) {
      if (seq.size() < this->K_)
        return false;

      RtSeq kmer = seq.start(this->K_) >> 'A';
      bool stop = false;
      for (size_t j = this->K_ - 1; j < seq.size(); ++j) {
        kmer <<= seq[j];
        if (!kmer_filter_.filter(kmer))
          continue;

        stop |= this->push_back_internal(kmer, thread_id);
      }

      return stop;
  }

 public:
  DeBruijnKMerSplitter(fs::TmpDir work_dir,
                       unsigned K, KmerFilter kmer_filter, size_t read_buffer_size = 0)
      : RtSeqKMerSplitter(work_dir, K), kmer_filter_(kmer_filter), read_buffer_size_(read_buffer_size) {
  }
 protected:
  DECL_LOGGER("DeBruijnKMerSplitter");
};

template<class Read, class KmerFilter>
class DeBruijnReadKMerSplitter : public DeBruijnKMerSplitter<KmerFilter> {
  io::ReadStreamList<Read>& streams_;

  template<class ReadStream>
  size_t
  FillBufferFromStream(ReadStream& stream, unsigned thread_id);

 public:
  using typename DeBruijnKMerSplitter<KmerFilter>::RawKMers;
  DeBruijnReadKMerSplitter(fs::TmpDir work_dir,
                           unsigned K,
                           io::ReadStreamList<Read>& streams,
                           size_t read_buffer_size = 0,
                           KmerFilter filter = KmerFilter())
      : DeBruijnKMerSplitter<KmerFilter>(work_dir, K, filter, read_buffer_size),
      streams_(streams) {}

  RawKMers Split(size_t num_files, unsigned nthreads) override;
};

template<class Read, class KmerFilter> template<class ReadStream>
size_t
DeBruijnReadKMerSplitter< Read, KmerFilter>::FillBufferFromStream(ReadStream &stream,
                                                                  unsigned thread_id) {
  typename ReadStream::ReadT r;
  size_t reads = 0;

  while (!stream.eof()) {
    stream >> r;
    reads += 1;

    if (this->FillBufferFromSequence(r.sequence(), thread_id))
      break;
  }

  return reads;
}

template<class Read, class KmerFilter>
typename DeBruijnReadKMerSplitter<Read, KmerFilter>::RawKMers
DeBruijnReadKMerSplitter<Read, KmerFilter>::Split(size_t num_files, unsigned nthreads) {
  auto out = this->PrepareBuffers(num_files, nthreads, this->read_buffer_size_);

  size_t counter = 0, n = 15;
  streams_.reset();
  while (!streams_.eof()) {
#   pragma omp parallel for num_threads(nthreads) reduction(+ : counter)
    for (unsigned i = 0; i < (unsigned)streams_.size(); ++i) {
      counter += FillBufferFromStream(streams_[i], omp_get_thread_num());
    }

    this->DumpBuffers(out);

    if (counter >> n) {
      INFO("Processed " << counter << " reads");
      n += 1;
    }
  }

  this->ClearBuffers();
  INFO("Used " << counter << " reads");
  return out;
}

template<class KmerFilter, class KMerIterator>
class DeBruijnKMerKMerSplitter : public DeBruijnKMerSplitter<KmerFilter> {
  using kmer_range = adt::iterator_range<KMerIterator>;
  unsigned K_source_;
  std::vector<kmer_range> kmers_;
  bool add_rc_;

  size_t FillBufferFromKMers(kmer_range &range, size_t thread_id);

 public:
  using typename DeBruijnKMerSplitter<KmerFilter>::RawKMers;

  DeBruijnKMerKMerSplitter(fs::TmpDir work_dir,
                           unsigned K_target, unsigned K_source, bool add_rc, size_t read_buffer_size = 0)
      : DeBruijnKMerSplitter<KmerFilter>(work_dir, K_target, KmerFilter(), read_buffer_size),
        K_source_(K_source), add_rc_(add_rc) {}

  void AddKMers(kmer_range range) {
    kmers_.emplace_back(std::move(range));
  }

  RawKMers Split(size_t num_files, unsigned nthreads) override;
};

template<class KmerFilter, class KMerIterator>
inline size_t DeBruijnKMerKMerSplitter<KmerFilter, KMerIterator>::FillBufferFromKMers(kmer_range &range,
                                                                                      size_t thread_id) {
  size_t seqs = 0;
  for (auto &it = range.begin(); it != range.end(); ++it) {
    RtSeq nucls(K_source_, it->first); // FIXME: temporary
    seqs += 1;

    bool stop = this->FillBufferFromSequence(nucls, unsigned(thread_id));
    if (add_rc_)
      stop |= this->FillBufferFromSequence(!nucls, unsigned(thread_id));

    if (stop)
      break;
  }

  return seqs;
}

template<class KmerFilter, class KMerIterator>
typename DeBruijnKMerKMerSplitter<KmerFilter,  KMerIterator>::RawKMers
DeBruijnKMerKMerSplitter<KmerFilter, KMerIterator>::Split(size_t num_files, unsigned nthreads) {
  auto out = this->PrepareBuffers(num_files, nthreads, this->read_buffer_size_);

  size_t counter = 0, n = 10;

  while (!std::all_of(kmers_.begin(), kmers_.end(),
                      [](const kmer_range &r) { return r.begin() == r.end(); })) {
#   pragma omp parallel for num_threads(nthreads) reduction(+ : counter)
    for (size_t i = 0; i < kmers_.size(); ++i)
        counter += FillBufferFromKMers(kmers_[i], omp_get_thread_num());

    this->DumpBuffers(out);

    if (counter >> n) {
      INFO("Processed " << counter << " kmers");
      n += 1;
    }
  }

  INFO("Used " << counter << " kmers.");

  this->ClearBuffers();

  return out;
}


}
