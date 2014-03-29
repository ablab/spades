#pragma once
/*
 * kmer_splitters.hpp
 *
 *  Created on: May 24, 2013
 *      Author: anton
 */

#include "io/io_helper.hpp"
#include "storing_traits.hpp"

#include "file_limit.hpp"

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

// used for temporary reads storage during parallel reading
static const size_t READS_BUFFER_SIZE = 536870912; // 512 MB in bytes

typedef ::KMerSplitter<runtime_k::RtSeq> RtSeqKMerSplitter;

typedef KMerVector<runtime_k::RtSeq> RtSeqKMerVector;
typedef std::vector<RtSeqKMerVector> KMerBuffer;

template<class KmerFilter>
class DeBruijnKMerSplitter : public RtSeqKMerSplitter {
 private:
  bool skip_not_minimal_;
  KmerFilter kmer_filter_;
 protected:
  size_t read_buffer_size_;
 protected:
  size_t FillBufferFromSequence(const Sequence &seq,
                                KMerBuffer &buffer, unsigned num_files) const {
      size_t kmers = 0;

      if (seq.size() < this->K_)
        return kmers;

      runtime_k::RtSeq kmer = seq.start<runtime_k::RtSeq>(this->K_) >> 'A';
      for (size_t j = this->K_ - 1; j < seq.size(); ++j) {
        kmer <<= seq[j];
        if(kmer_filter_.filter(kmer)) {
            buffer[this->GetFileNumForSeq(kmer, num_files)].push_back(kmer);
            kmers++;
        }
      }
      return kmers;
  }


  void DumpBuffers(size_t num_files, size_t nthreads,
                   std::vector<KMerBuffer> &buffers,
                   const path::files_t &ostreams) const{
    # pragma omp parallel for
      for (unsigned k = 0; k < num_files; ++k) {
        size_t sz = 0;
        for (size_t i = 0; i < nthreads; ++i)
          sz += buffers[i][k].size();

        KMerVector<runtime_k::RtSeq> SortBuffer(this->K_, sz);
        for (size_t i = 0; i < nthreads; ++i) {
          KMerBuffer &entry = buffers[i];
          for (size_t j = 0; j < entry[k].size(); ++j)
            SortBuffer.push_back(entry[k][j]);
        }
        libcxx::sort(SortBuffer.begin(), SortBuffer.end(), KMerVector<runtime_k::RtSeq>::less2_fast());
        auto it = std::unique(SortBuffer.begin(), SortBuffer.end(), KMerVector<runtime_k::RtSeq>::equal_to());

    #   pragma omp critical
        {
          FILE *f = fopen(ostreams[k].c_str(), "ab");
          VERIFY_MSG(f, "Cannot open temporary file to write");
          fwrite(SortBuffer.data(), SortBuffer.el_data_size(), it - SortBuffer.begin(), f);
          fclose(f);
        }
      }

      for (unsigned i = 0; i < nthreads; ++i) {
        for (unsigned j = 0; j < num_files; ++j) {
          buffers[i][j].clear();
        }
      }
  }

 public:
  DeBruijnKMerSplitter(const std::string &work_dir,
                       unsigned K, KmerFilter kmer_filter, size_t read_buffer_size = 0, uint32_t seed = 0)
      : RtSeqKMerSplitter(work_dir, K, seed), kmer_filter_(kmer_filter), read_buffer_size_(read_buffer_size) {
  }
 protected:
  DECL_LOGGER("DeBruijnKMerSplitter");
};

template<class Read, class KmerFilter>
class DeBruijnReadKMerSplitter : public DeBruijnKMerSplitter<KmerFilter> {
  io::ReadStreamList<Read> &streams_;
  io::SingleStream *contigs_;

  template<class ReadStream>
  std::pair<size_t, size_t>
  FillBufferFromStream(ReadStream& stream,
                       KMerBuffer &tmp_entries,
                       unsigned num_files, size_t cell_size) const;

  size_t rl_;

 public:
  DeBruijnReadKMerSplitter(const std::string &work_dir,
                           unsigned K, uint32_t seed,
                           io::ReadStreamList<Read>& streams,
                           io::SingleStream* contigs_stream = 0,
                           size_t read_buffer_size = 0)
      : DeBruijnKMerSplitter<KmerFilter>(work_dir, K, KmerFilter(), read_buffer_size, seed),
        streams_(streams), contigs_(contigs_stream), rl_(0) {
  }

  virtual path::files_t Split(size_t num_files);

  size_t read_length() const { return rl_; }
};

template<class Read, class KmerFilter> template<class ReadStream>
std::pair<size_t, size_t>
DeBruijnReadKMerSplitter<Read, KmerFilter>::FillBufferFromStream(ReadStream &stream,
                                                     KMerBuffer &buffer,
                                                     unsigned num_files, size_t cell_size) const {
  typename ReadStream::ReadT r;
  size_t reads = 0, kmers = 0, rl = 0;

  while (!stream.eof() && kmers < num_files * cell_size) {
    stream >> r;
    rl = std::max(rl, r.size());
    reads += 1;

    kmers += this->FillBufferFromSequence(r.sequence(), buffer, num_files);
  }
  return std::make_pair(reads, rl);
}

template<class Read, class KmerFilter>
path::files_t DeBruijnReadKMerSplitter<Read, KmerFilter>::Split(size_t num_files) {
  unsigned nthreads = (unsigned) streams_.size();

  INFO("Splitting kmer instances into " << num_files << " buckets. This might take a while.");

  // Determine the set of output files
  path::files_t out;
  for (unsigned i = 0; i < num_files; ++i)
    out.push_back(this->GetRawKMersFname(i));

  size_t file_limit = num_files + 2*nthreads;
  size_t res = limit_file(file_limit);
  if (res < file_limit) {
    WARN("Failed to setup necessary limit for number of open files. The process might crash later on.");
    WARN("Do 'ulimit -n " << file_limit << "' in the console to overcome the limit");
  }

  size_t reads_buffer_size = DeBruijnKMerSplitter<KmerFilter>::read_buffer_size_;
  if (reads_buffer_size == 0) {
    reads_buffer_size = READS_BUFFER_SIZE;
    size_t mem_limit =  (size_t)((double)(get_free_memory()) / (nthreads * 3));
    INFO("Memory available for splitting buffers: " << (double)mem_limit / 1024.0 / 1024.0 / 1024.0 << " Gb");
    reads_buffer_size = std::min(reads_buffer_size, mem_limit);
  }
  size_t cell_size = reads_buffer_size /
                     (num_files * runtime_k::RtSeq::GetDataSize(this->K_) * sizeof(runtime_k::RtSeq::DataType));

  // Set sane minimum cell size
  if (cell_size < 16384)
    cell_size = 16384;
  INFO("Using cell size of " << cell_size);

  std::vector<KMerBuffer> tmp_entries(nthreads);
  for (unsigned i = 0; i < nthreads; ++i) {
    KMerBuffer &entry = tmp_entries[i];
    entry.resize(num_files, RtSeqKMerVector(this->K_, (size_t) (1.1 * (double) cell_size)));
  }

  size_t counter = 0, rl = 0, n = 15;
  streams_.reset();
  while (!streams_.eof()) {
#   pragma omp parallel for num_threads(nthreads) reduction(+ : counter) shared(rl)
    for (size_t i = 0; i < nthreads; ++i) {
      std::pair<size_t, size_t> stats = FillBufferFromStream(streams_[i], tmp_entries[i], (unsigned) num_files, cell_size);
      counter += stats.first;

      // There is no max reduction in C/C++ OpenMP... Only in FORTRAN :(
#     pragma omp flush(rl)
      if (stats.second > rl)
#     pragma omp critical
      {
        rl = std::max(rl, stats.second);
      }
    }

    this->DumpBuffers(num_files, nthreads, tmp_entries, out);

    if (counter >> n) {
      INFO("Processed " << counter << " reads");
      n += 1;
    }
  }

  if (contigs_) {
    INFO("Adding contigs from previous K");
    size_t cnt = 0;
    contigs_->reset();
    while (!contigs_->eof()) {
      FillBufferFromStream(*contigs_, tmp_entries[cnt], (unsigned) num_files, cell_size);
      this->DumpBuffers(num_files, nthreads, tmp_entries, out);
      if (++cnt >= nthreads)
        cnt = 0;
    }
  }

  INFO("Used " << counter << " reads. Maximum read length " << rl);
  rl_ = rl;

  return out;
}

template<class Graph, class KmerFilter>
class DeBruijnGraphKMerSplitter : public DeBruijnKMerSplitter<KmerFilter> {
  typedef typename Graph::ConstEdgeIt EdgeIt;
  typedef typename Graph::EdgeId EdgeId;

  const Graph &g_;

  size_t FillBufferFromEdges(EdgeIt &edge,
                             KMerBuffer &tmp_entries,
                             unsigned num_files, size_t cell_size) const;

 public:
  DeBruijnGraphKMerSplitter(const std::string &work_dir,
                            unsigned K, const Graph &g, size_t read_buffer_size = 0)
      : DeBruijnKMerSplitter<KmerFilter>(work_dir, K, KmerFilter(), read_buffer_size), g_(g) {}

  virtual path::files_t Split(size_t num_files);
};

template<class Graph, class KmerFilter>
size_t
DeBruijnGraphKMerSplitter<Graph, KmerFilter>::FillBufferFromEdges(EdgeIt &edge,
                                                      KMerBuffer &buffer,
                                                      unsigned num_files, size_t cell_size) const {
  size_t seqs = 0;
  for (size_t kmers = 0; !edge.IsEnd() && kmers < num_files * cell_size; ++edge) {
    const Sequence &nucls = g_.EdgeNucls(*edge);

    kmers += this->FillBufferFromSequence(nucls, buffer, num_files);
    seqs += 1;
  }

  return seqs;
}

template<class Graph, class KmerFilter>
path::files_t DeBruijnGraphKMerSplitter<Graph, KmerFilter>::Split(size_t num_files) {
  INFO("Splitting kmer instances into " << num_files << " buckets. This might take a while.");

  // Determine the set of output files
  path::files_t out;
  for (unsigned i = 0; i < num_files; ++i)
    out.push_back(this->GetRawKMersFname(i));

  size_t file_limit = num_files + 2*16;
  size_t res = limit_file(file_limit);
  if (res < file_limit) {
    WARN("Failed to setup necessary limit for number of open files. The process might crash later on.");
    WARN("Do 'ulimit -n " << file_limit << "' in the console to overcome the limit");
  }

  size_t reads_buffer_size = DeBruijnKMerSplitter<KmerFilter>::read_buffer_size_;
  if (reads_buffer_size == 0) {
    reads_buffer_size = READS_BUFFER_SIZE;
    size_t mem_limit =  (size_t)((double)(get_free_memory()) / (3));
    INFO("Memory available for splitting buffers: " << (double)mem_limit / 1024.0 / 1024.0 / 1024.0 << " Gb");
    reads_buffer_size = std::min(reads_buffer_size, mem_limit);
  }
  size_t cell_size = reads_buffer_size /
                     (num_files * runtime_k::RtSeq::GetDataSize(this->K_) * sizeof(runtime_k::RtSeq::DataType));
  INFO("Using cell size of " << cell_size);

  std::vector<KMerBuffer> tmp_entries(1);
  KMerBuffer &entry = tmp_entries[0];
  entry.resize(num_files, RtSeqKMerVector(this->K_, (size_t) (1.1 * (double) cell_size)));

  size_t counter = 0, n = 10;
  for (auto it = g_.ConstEdgeBegin(); !it.IsEnd(); ) {
    counter += FillBufferFromEdges(it, tmp_entries[0], (unsigned) num_files, cell_size);

    this->DumpBuffers(num_files, 1, tmp_entries, out);

    if (counter >> n) {
      INFO("Processed " << counter << " edges");
      n += 1;
    }
  }

  INFO("Used " << counter << " sequences.");

  return out;
}


template<class KmerFilter>
class DeBruijnKMerKMerSplitter : public DeBruijnKMerSplitter<KmerFilter> {
  typedef MMappedFileRecordArrayIterator<runtime_k::RtSeq::DataType> kmer_iterator;

  unsigned K_source_;
  std::vector<std::string> kmers_;
  bool add_rc_;

  size_t FillBufferFromKMers(kmer_iterator &kmer,
                             KMerBuffer &tmp_entries,
                             unsigned num_files, size_t cell_size) const;

 public:
  DeBruijnKMerKMerSplitter(const std::string &work_dir,
                           unsigned K_target, unsigned K_source, bool add_rc, size_t read_buffer_size = 0)
      : DeBruijnKMerSplitter<KmerFilter>(work_dir, K_target, KmerFilter(), read_buffer_size), K_source_(K_source), add_rc_(add_rc) {}

  void AddKMers(const std::string &file) {
    kmers_.push_back(file);
  }

  virtual path::files_t Split(size_t num_files);
};

template<class KmerFilter>
inline size_t DeBruijnKMerKMerSplitter<KmerFilter>::FillBufferFromKMers(kmer_iterator &kmer,
                                                     KMerBuffer &buffer,
                                                     unsigned num_files, size_t cell_size) const {
  size_t seqs = 0;
  for (size_t kmers = 0; kmer.good() && kmers < num_files * cell_size; ++kmer) {
    Sequence nucls(runtime_k::RtSeq(K_source_, *kmer));
    kmers += this->FillBufferFromSequence(nucls, buffer, num_files);
    if(add_rc_)
      kmers += this->FillBufferFromSequence(!nucls, buffer, num_files);
    seqs += 1;
  }

  return seqs;
}

template<class KmerFilter>
inline path::files_t DeBruijnKMerKMerSplitter<KmerFilter>::Split(size_t num_files) {
  unsigned nthreads = (unsigned) kmers_.size();
  INFO("Splitting kmer instances into " << num_files << " buckets. This might take a while.");

  // Determine the set of output files
  path::files_t out;
  for (unsigned i = 0; i < num_files; ++i)
    out.push_back(this->GetRawKMersFname(i));

  size_t file_limit = num_files + 2*nthreads;
  size_t res = limit_file(file_limit);
  if (res < file_limit) {
    WARN("Failed to setup necessary limit for number of open files. The process might crash later on.");
    WARN("Do 'ulimit -n " << file_limit << "' in the console to overcome the limit");
  }

  size_t reads_buffer_size = DeBruijnKMerSplitter<KmerFilter>::read_buffer_size_;
  if (reads_buffer_size == 0) {
    reads_buffer_size = READS_BUFFER_SIZE;
    size_t mem_limit =  (size_t)((double)(get_free_memory()) / (nthreads * 3));
    INFO("Memory available for splitting buffers: " << (double)mem_limit / 1024.0 / 1024.0 / 1024.0 << " Gb");
    reads_buffer_size = std::min(reads_buffer_size, mem_limit);
  }
  size_t cell_size = reads_buffer_size /
                     (num_files * runtime_k::RtSeq::GetDataSize(this->K_) * sizeof(runtime_k::RtSeq::DataType));
  // Set sane minimum cell size
  if (cell_size < 16384)
    cell_size = 16384;
  INFO("Using cell size of " << cell_size);

  std::vector<KMerBuffer> tmp_entries(nthreads);
  for (unsigned i = 0; i < nthreads; ++i) {
    KMerBuffer &entry = tmp_entries[i];
    entry.resize(num_files, RtSeqKMerVector(this->K_, (size_t) (1.1 * (double) cell_size)));
  }

  size_t counter = 0, n = 10;
  std::vector<kmer_iterator> its;
  its.reserve(nthreads);
  for (auto it = kmers_.begin(), et = kmers_.end(); it != et; ++it)
    its.emplace_back(*it, runtime_k::RtSeq::GetDataSize(K_source_));

  bool anygood = false;
  do {
#   pragma omp parallel for num_threads(nthreads) reduction(+ : counter)
    for (size_t i = 0; i < nthreads; ++i)
      counter += FillBufferFromKMers(its[i], tmp_entries[i], (unsigned) num_files, cell_size);

    this->DumpBuffers(num_files, nthreads, tmp_entries, out);

    if (counter >> n) {
      INFO("Processed " << counter << " kmers");
      n += 1;
    }

    anygood = false;
    for (auto it = its.begin(), et = its.end(); it != et; ++it)
      anygood |= it->good();
  } while (anygood);

  INFO("Used " << counter << " kmers.");

  return out;
}


}
