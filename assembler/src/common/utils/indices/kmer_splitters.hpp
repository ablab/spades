//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "io/reads/io_helper.hpp"
#include "storing_traits.hpp"

#include "utils/file_limit.hpp"
#include "utils/mph_index/kmer_index_builder.hpp"

namespace debruijn_graph {

template<class StoringType>
struct StoringTypeFilter {
};

template<>
struct StoringTypeFilter<SimpleStoring> {
    template<class Kmer>
    bool filter(const Kmer &/*kmer*/) const {
        return true;
    }
};

template<>
struct StoringTypeFilter<InvertableStoring> {
    template<class Kmer>
    bool filter(const Kmer &kmer) const {
        return kmer.IsMinimal();
    }
};

using RtSeqKMerSplitter = ::KMerSortingSplitter<RtSeq>;

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

 public:
  DeBruijnKMerSplitter(const std::string &work_dir,
                       unsigned K, KmerFilter kmer_filter, size_t read_buffer_size = 0, uint32_t seed = 0)
      : RtSeqKMerSplitter(work_dir, K, seed), kmer_filter_(kmer_filter), read_buffer_size_(read_buffer_size) {
  }
 protected:
  DECL_LOGGER("DeBruijnKMerSplitter");
};

struct ReadStatistics {
  size_t reads_;
  size_t max_read_length_;
  size_t bases_;
};

template<class Read, class KmerFilter>
class DeBruijnReadKMerSplitter : public DeBruijnKMerSplitter<KmerFilter> {
  io::ReadStreamList<Read> &streams_;
  io::SingleStream *contigs_;

  template<class ReadStream>
  ReadStatistics
  FillBufferFromStream(ReadStream& stream, unsigned thread_id);

  ReadStatistics rs_;

 public:
  DeBruijnReadKMerSplitter(const std::string &work_dir,
                           unsigned K, uint32_t seed,
                           io::ReadStreamList<Read>& streams,
                           io::SingleStream* contigs_stream = 0,
                           size_t read_buffer_size = 0)
      : DeBruijnKMerSplitter<KmerFilter>(work_dir, K, KmerFilter(), read_buffer_size, seed),
      streams_(streams), contigs_(contigs_stream), rs_({0 ,0 ,0}) {}

  path::files_t Split(size_t num_files) override;

  size_t read_length() const { return rs_.max_read_length_; }
  ReadStatistics stats() const { return rs_; }
};

template<class Read, class KmerFilter> template<class ReadStream>
ReadStatistics
DeBruijnReadKMerSplitter<Read, KmerFilter>::FillBufferFromStream(ReadStream &stream,
                                                                 unsigned thread_id) {
  typename ReadStream::ReadT r;
  size_t reads = 0, rl = 0, bases = 0;

  while (!stream.eof()) {
    stream >> r;
    rl = std::max(rl, r.size());
    reads += 1;
    bases += r.size();

    if (this->FillBufferFromSequence(r.sequence(), thread_id))
      break;
  }
  return { reads, rl, bases };
}

template<class Read, class KmerFilter>
path::files_t DeBruijnReadKMerSplitter<Read, KmerFilter>::Split(size_t num_files) {
  unsigned nthreads = (unsigned) streams_.size();

  INFO("Splitting kmer instances into " << num_files << " buckets. This might take a while.");
  path::files_t out = this->PrepareBuffers(num_files, nthreads, this->read_buffer_size_);

  size_t counter = 0, rl = 0, bases = 0, n = 15;
  streams_.reset();
  while (!streams_.eof()) {
#   pragma omp parallel for num_threads(nthreads) reduction(+ : counter) reduction(+ : bases) shared(rl)
    for (unsigned i = 0; i < nthreads; ++i) {
      ReadStatistics stats = FillBufferFromStream(streams_[i], i);
      counter += stats.reads_;
      bases += stats.bases_;

      // There is no max reduction in C/C++ OpenMP... Only in FORTRAN :(
#     pragma omp flush(rl)
      if (stats.max_read_length_ > rl)
#     pragma omp critical
      {
        rl = std::max(rl, stats.max_read_length_);
      }
    }

    this->DumpBuffers(out);

    if (counter >> n) {
      INFO("Processed " << counter << " reads");
      n += 1;
    }
  }

  if (contigs_) {
    INFO("Adding contigs from previous K");
    unsigned cnt = 0;
    contigs_->reset();
    while (!contigs_->eof()) {
      FillBufferFromStream(*contigs_, cnt);
      this->DumpBuffers(out);
      if (++cnt >= nthreads)
        cnt = 0;
    }
  }

  this->ClearBuffers();

  INFO("Used " << counter << " reads. Maximum read length " << rl);
  INFO("Average read length " << double(bases) / double(counter));
  rs_ = { counter, rl, bases };

  return out;
}

template<class Graph, class KmerFilter>
class DeBruijnGraphKMerSplitter : public DeBruijnKMerSplitter<KmerFilter> {
  typedef typename Graph::ConstEdgeIt EdgeIt;
  typedef typename Graph::EdgeId EdgeId;

  const Graph &g_;

  size_t FillBufferFromEdges(EdgeIt &edge, unsigned thread_id);

 public:
  DeBruijnGraphKMerSplitter(const std::string &work_dir,
                            unsigned K, const Graph &g, size_t read_buffer_size = 0)
      : DeBruijnKMerSplitter<KmerFilter>(work_dir, K, KmerFilter(), read_buffer_size), g_(g) {}

  path::files_t Split(size_t num_files) override;
};

template<class Graph, class KmerFilter>
size_t
DeBruijnGraphKMerSplitter<Graph, KmerFilter>::FillBufferFromEdges(EdgeIt &edge,
                                                                  unsigned thread_id) {
  size_t seqs = 0;
  for (; !edge.IsEnd(); ++edge) {
    const Sequence &nucls = g_.EdgeNucls(*edge);

    seqs += 1;
    if (this->FillBufferFromSequence(nucls, thread_id))
      break;
  }

  return seqs;
}

template<class Graph, class KmerFilter>
path::files_t DeBruijnGraphKMerSplitter<Graph, KmerFilter>::Split(size_t num_files) {
  INFO("Splitting kmer instances into " << num_files << " buckets. This might take a while.");

  path::files_t out = this->PrepareBuffers(num_files, 1, this->read_buffer_size_);

  size_t counter = 0, n = 10;
  for (auto it = g_.ConstEdgeBegin(); !it.IsEnd(); ) {
    counter += FillBufferFromEdges(it, 0);

    this->DumpBuffers(out);

    if (counter >> n) {
      INFO("Processed " << counter << " edges");
      n += 1;
    }
  }

  INFO("Used " << counter << " sequences.");

  this->ClearBuffers();
  
  return out;
}


template<class KmerFilter>
class DeBruijnKMerKMerSplitter : public DeBruijnKMerSplitter<KmerFilter> {
  typedef MMappedFileRecordArrayIterator<RtSeq::DataType> kmer_iterator;

  unsigned K_source_;
  std::vector<std::string> kmers_;
  bool add_rc_;

  size_t FillBufferFromKMers(kmer_iterator &kmer,
                             unsigned thread_id);

 public:
  DeBruijnKMerKMerSplitter(const std::string &work_dir,
                           unsigned K_target, unsigned K_source, bool add_rc, size_t read_buffer_size = 0)
      : DeBruijnKMerSplitter<KmerFilter>(work_dir, K_target, KmerFilter(), read_buffer_size),
        K_source_(K_source), add_rc_(add_rc) {}

  void AddKMers(const std::string &file) {
    kmers_.push_back(file);
  }

  path::files_t Split(size_t num_files) override;
};

template<class KmerFilter>
inline size_t DeBruijnKMerKMerSplitter<KmerFilter>::FillBufferFromKMers(kmer_iterator &kmer,
                                                                        unsigned thread_id) {
  size_t seqs = 0;
  for (; kmer.good(); ++kmer) {
    Sequence nucls(RtSeq(K_source_, *kmer));
    seqs += 1;

    bool stop = this->FillBufferFromSequence(nucls, thread_id);
    if (add_rc_)
      stop |= this->FillBufferFromSequence(!nucls, thread_id);

    if (stop)
      break;
  }

  return seqs;
}

template<class KmerFilter>
path::files_t DeBruijnKMerKMerSplitter<KmerFilter>::Split(size_t num_files) {
  unsigned nthreads = (unsigned) kmers_.size();

  INFO("Splitting kmer instances into " << num_files << " buckets. This might take a while.");

  path::files_t out = this->PrepareBuffers(num_files, nthreads, this->read_buffer_size_);

  size_t counter = 0, n = 10;
  std::vector<kmer_iterator> its;
  its.reserve(nthreads);
  for (auto it = kmers_.begin(), et = kmers_.end(); it != et; ++it)
    its.emplace_back(*it, RtSeq::GetDataSize(K_source_));

  while (std::any_of(its.begin(), its.end(),
                     [](const kmer_iterator &it) { return it.good(); })) {
#   pragma omp parallel for num_threads(nthreads) reduction(+ : counter)
    for (unsigned i = 0; i < nthreads; ++i)
      counter += FillBufferFromKMers(its[i], i);

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
