//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "io/reads/io_helper.hpp"
#include "utils/filesystem/file_limit.hpp"

namespace utils {

template<class Seq>
class KMerSplitter {
public:
    typedef typename Seq::hash hash_function;

    KMerSplitter(const std::string &work_dir, unsigned K, uint32_t seed = 0)
            : work_dir_(work_dir), K_(K), seed_(seed) {}

    virtual ~KMerSplitter() {}

    virtual fs::files_t Split(size_t num_files, unsigned nthreads) = 0;

    size_t kmer_size() const {
        return Seq::GetDataSize(K_) * sizeof(typename Seq::DataType);
    }

    unsigned K() const { return K_; }

protected:
    const std::string &work_dir_;
    hash_function hash_;
    unsigned K_;
    uint32_t seed_;

    DECL_LOGGER("K-mer Splitting");
};

template<class Seq>
class KMerSortingSplitter : public KMerSplitter<Seq> {
public:
    KMerSortingSplitter(const std::string &work_dir, unsigned K, uint32_t seed = 0)
            : KMerSplitter<Seq>(work_dir, K, seed), cell_size_(0), num_files_(0) {}

protected:
    using SeqKMerVector = adt::KMerVector<Seq>;
    using KMerBuffer = std::vector<SeqKMerVector>;

    std::vector<KMerBuffer> kmer_buffers_;
    size_t cell_size_;
    size_t num_files_;

    fs::files_t PrepareBuffers(size_t num_files, unsigned nthreads, size_t reads_buffer_size) {
        num_files_ = num_files;

        // Determine the set of output files
        fs::files_t out;
        for (unsigned i = 0; i < num_files_; ++i)
            out.push_back(this->GetRawKMersFname(i));

        size_t file_limit = num_files_ + 2*nthreads;
        size_t res = limit_file(file_limit);
        if (res < file_limit) {
            WARN("Failed to setup necessary limit for number of open files. The process might crash later on.");
            WARN("Do 'ulimit -n " << file_limit << "' in the console to overcome the limit");
        }

        if (reads_buffer_size == 0) {
            reads_buffer_size = 536870912ull;
            size_t mem_limit =  (size_t)((double)(utils::get_free_memory()) / (nthreads * 3));
            INFO("Memory available for splitting buffers: " << (double)mem_limit / 1024.0 / 1024.0 / 1024.0 << " Gb");
            reads_buffer_size = std::min(reads_buffer_size, mem_limit);
        }
        cell_size_ = reads_buffer_size / (num_files_ * this->kmer_size());
        // Set sane minimum cell size
        if (cell_size_ < 16384)
            cell_size_ = 16384;

        INFO("Using cell size of " << cell_size_);
        kmer_buffers_.resize(nthreads);
        for (unsigned i = 0; i < nthreads; ++i) {
            KMerBuffer &entry = kmer_buffers_[i];
            entry.resize(num_files_, adt::KMerVector<Seq>(this->K_, (size_t) (1.1 * (double) cell_size_)));
        }

        return out;
    }

    bool push_back_internal(const Seq &seq, unsigned thread_id) {
        KMerBuffer &entry = kmer_buffers_[thread_id];

        size_t idx = this->GetFileNumForSeq(seq, (unsigned)num_files_);
        entry[idx].push_back(seq);
        return entry[idx].size() > cell_size_;
    }

    void DumpBuffers(const fs::files_t &ostreams) {
        VERIFY(ostreams.size() == num_files_ && kmer_buffers_[0].size() == num_files_);

#   pragma omp parallel for
        for (unsigned k = 0; k < num_files_; ++k) {
            // Below k is thread id!

            size_t sz = 0;
            for (size_t i = 0; i < kmer_buffers_.size(); ++i)
                sz += kmer_buffers_[i][k].size();

            adt::KMerVector<Seq> SortBuffer(this->K_, sz);
            for (auto & entry : kmer_buffers_) {
                const auto &buffer = entry[k];
                for (size_t j = 0; j < buffer.size(); ++j)
                    SortBuffer.push_back(buffer[j]);
            }
            libcxx::sort(SortBuffer.begin(), SortBuffer.end(), typename adt::KMerVector<Seq>::less2_fast());
            auto it = std::unique(SortBuffer.begin(), SortBuffer.end(), typename adt::KMerVector<Seq>::equal_to());

#     pragma omp critical
            {
                size_t cnt =  it - SortBuffer.begin();

                // Write k-mers
                FILE *f = fopen(ostreams[k].c_str(), "ab");
                VERIFY_MSG(f, "Cannot open temporary file to write");
                fwrite(SortBuffer.data(), SortBuffer.el_data_size(), cnt, f);
                fclose(f);

                // Write index
                f = fopen((ostreams[k] + ".idx").c_str(), "ab");
                VERIFY_MSG(f, "Cannot open temporary file to write");
                fwrite(&cnt, sizeof(cnt), 1, f);
                fclose(f);
            }
        }

        for (auto & entry : kmer_buffers_)
            for (auto & eentry : entry)
                eentry.clear();
    }

    void ClearBuffers() {
        for (auto & entry : kmer_buffers_)
            for (auto & eentry : entry) {
                eentry.clear();
                eentry.shrink_to_fit();
            }
    }

    std::string GetRawKMersFname(unsigned suffix) const {
        return fs::append_path(this->work_dir_, "kmers.raw." + std::to_string(suffix));
    }

    unsigned GetFileNumForSeq(const Seq &s, unsigned total) const {
        return (unsigned)(this->hash_(s, this->seed_) % total);
    }

};

using RtSeqKMerSplitter = KMerSortingSplitter<RtSeq>;

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

  fs::files_t Split(size_t num_files, unsigned nthreads) override;

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
fs::files_t DeBruijnReadKMerSplitter<Read, KmerFilter>::Split(size_t num_files, unsigned nthreads) {
  fs::files_t out = this->PrepareBuffers(num_files, nthreads, this->read_buffer_size_);

  size_t counter = 0, rl = 0, bases = 0, n = 15;
  streams_.reset();
  while (!streams_.eof()) {
#   pragma omp parallel for num_threads(nthreads) reduction(+ : counter) reduction(+ : bases) shared(rl)
    for (unsigned i = 0; i < (unsigned)streams_.size(); ++i) {
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

    fs::files_t Split(size_t num_files, unsigned nthreads) override;
};

template<class KmerFilter>
inline size_t DeBruijnKMerKMerSplitter<KmerFilter>::FillBufferFromKMers(kmer_iterator &kmer,
                                                                        unsigned thread_id) {
  size_t seqs = 0;
  for (; kmer.good(); ++kmer) {
    RtSeq nucls(K_source_, *kmer);
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
fs::files_t DeBruijnKMerKMerSplitter<KmerFilter>::Split(size_t num_files, unsigned nthreads) {
  unsigned nit = (unsigned) kmers_.size();

  fs::files_t out = this->PrepareBuffers(num_files, nthreads, this->read_buffer_size_);

  size_t counter = 0, n = 10;
  std::vector<kmer_iterator> its;
  its.reserve(nit);
  for (auto it = kmers_.begin(), et = kmers_.end(); it != et; ++it)
    its.emplace_back(*it, RtSeq::GetDataSize(K_source_));

  while (std::any_of(its.begin(), its.end(),
                     [](const kmer_iterator &it) { return it.good(); })) {
#   pragma omp parallel for num_threads(nthreads) reduction(+ : counter)
    for (unsigned i = 0; i < nit; ++i)
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
