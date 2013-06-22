#pragma once

#include "debruijn_edge_index.hpp"

namespace debruijn_graph {

template <class Seq>
class OldDeBruijnEdgeIndexBuilder;

//template<class Graph, class Seq = runtime_k::RtSeq, class traits = kmer_index_traits<Seq>>
//class KmerStoringDeBruijnEdgeIndex : public DeBruijnKMerIndex<EdgeInfo<typename Graph::EdgeId>, traits> {
//    typedef typename Graph::EdgeId IdType;
//    typedef DeBruijnKMerIndex<EdgeInfo<IdType>, traits> base;
//    const Graph &graph_;
//  public:
//    typedef typename traits::SeqType KMer;
//    typedef KMerIndex<traits>        KMerIndexT;
//
//    KmerStoringDeBruijnEdgeIndex(unsigned K, const Graph &graph, const std::string &workdir)
//            : base(K, workdir), graph_(graph) {}
//
//    ~KmerStoringDeBruijnEdgeIndex() {}
//
//    bool contains(const KMer& kmer) const {
//        typename base::KMerIdx idx = base::seq_idx(kmer);
//        return contains(idx, kmer);
//    }
//
//    KMer kmer(typename base::KMerIdx idx) const {
//        VERIFY(idx < this->size());
//        const typename base::KMerIndexValueType &entry = base::operator[](idx);
//        VERIFY(entry.offset_ != -1);
//        return KMer(this->K_, graph_.EdgeNucls(entry.edgeId_), entry.offset_);
//    }
//
//    bool contains(typename base::KMerIdx idx, const KMer &k) const {
//        // Sanity check
//        if (idx == base::InvalidKMerIdx || idx >= this->size())
//            return false;
//
//        const typename base::KMerIndexValueType &entry = base::operator[](idx);
//
//        if (entry.offset_ == -1)
//            return false;
//
//        return k == KMer(this->K_, graph_.EdgeNucls(entry.edgeId_), entry.offset_);
//    }
//
//    //    bool contains(KMerIdx idx, const KMer &k,
//    //                  bool check_push_back = true) const {
//    //               // Sanity check
//    //               if (idx == InvalidKMerIdx || idx >= size())
//    //                       return false;
//    //
//    //               if (idx < data_.size()) {
//    //                       auto it = kmers->begin() + idx;
//    //                       return (typename traits::raw_equal_to()(k, *it));
//    //               }
//    //
//    //               if (check_push_back) {
//    //                       auto it = push_back_index_.right.find(idx - data_.size());
//    //                       return (it != push_back_index_.right.end() && it->second == k);
//    //               }
//    //
//    //               return false;
//    //       }
//
//  public:
//    /**
//     * Number of edges coming into param edge's end
//     */
//    unsigned RivalEdgeCount(const KMer &kmer) const {
//        KMer kmer2 = kmer << 'A';
//        unsigned res = 0;
//        for (char c = 0; c < 4; ++c)
//            if (RandomAccessContains(kmer2 >> c))
//                res += 1;
//
//        return res;
//    }
//  private:
//    bool RandomAccessContains(const KMer &kmer) const {
//        VERIFY(false);
//        return false;
//    }
//
//    KMer RandomAccesskmer(typename base::KMerIdx idx) const {
//        VERIFY(false);
//        return KMer();
//    }
//
//  public:
//    unsigned RivalEdgeCount(typename base::KMerIdx idx) const {
//        KMer kmer2 = RandomAccesskmer(idx) << 'A';
//        unsigned res = 0;
//        for (char c = 0; c < 4; ++c)
//            if (RandomAccessContains(kmer2 >> c))
//                res += 1;
//
//        return res;
//    }
//
//    /**
//     * Number of edges going out of the param edge's end
//     */
//    unsigned NextEdgeCount(const KMer &kmer) const {
//        unsigned res = 0;
//        for (char c = 0; c < 4; ++c)
//            if (RandomAccessContains(kmer << c))
//                res += 1;
//
//        return res;
//    }
//
//    unsigned NextEdgeCount(typename base::KMerIdx idx) const {
//        KMer kmer = RandomAccesskmer(idx);
//
//        unsigned res = 0;
//        for (char c = 0; c < 4; ++c)
//            if (RandomAccessContains(kmer << c))
//                res += 1;
//
//        return res;
//    }
//
//    KMer NextEdge(const KMer &kmer) const { // returns any next edge
//        for (char c = 0; c < 4; ++c) {
//            KMer s = kmer << c;
//            typename base::KMerIdx idx = base::seq_idx(s);
//            if (RandomAccessContains(idx))
//                return RandomAccesskmer(idx);
//        }
//
//        VERIFY_MSG(false, "Couldn't find requested edge!");
//        return KMer(base::K());
//        // no next edges (we should request one here).
//    }
//
//    KMer NextEdge(typename base::KMerIdx idx) const { // returns any next edge
//        KMer kmer = RandomAccesskmer(idx);
//
//        for (char c = 0; c < 4; ++c) {
//            KMer s = kmer << c;
//            if (RandomAccessContains(s))
//                return s;
//        }
//
//        VERIFY_MSG(false, "Couldn't find requested edge!");
//        return KMer(base::K());
//        // no next edges (we should request one here).
//    }
//
//    std::pair<IdType, size_t> get(const KMer &kmer) const {
//        typename base::KMerIdx idx = base::seq_idx(kmer);
//        VERIFY(this->contains(idx, kmer));
//
//        const EdgeInfo<IdType> &entry = base::operator[](idx);
//        return std::make_pair(entry.edgeId_, (size_t)entry.offset_);
//    }
//
//    std::pair<IdType, size_t> get(typename base::KMerIdx idx) const {
//        const EdgeInfo<IdType> &entry = base::operator[](idx);
//        return std::make_pair(entry.edgeId_, (size_t)entry.offset_);
//    }
//
//    bool DeleteIfEqual(const KMer &kmer, IdType id) {
//        typename base::KMerIdx idx = base::seq_idx(kmer);
//        if (!contains(idx, kmer))
//            return false;
//
//        EdgeInfo<IdType> &entry = base::operator[](idx);
//        if (entry.edgeId_ == id) {
//            entry.offset_ = -1;
//            return true;
//        }
//        return false;
//    }
//
//    void RenewKMers(const Sequence &nucls, IdType id) {
//        VERIFY(nucls.size() >= base::K());
//        KMer kmer(base::K(), nucls);
//
//        PutInIndex(kmer, id, 0);
//        for (size_t i = base::K(), n = nucls.size(); i < n; ++i) {
//            kmer <<= nucls[i];
//            PutInIndex(kmer, id, i - base::K() + 1);
//        }
//    }
//
//    void DeleteKMers(const Sequence &nucls, IdType id) {
//        VERIFY(nucls.size() >= base::K());
//        KMer kmer(base::K(), nucls);
//        DeleteIfEqual(kmer, id);
//        for (size_t i = base::K(), n = nucls.size(); i < n; ++i) {
//            kmer <<= nucls[i];
//            DeleteIfEqual(kmer, id);
//        }
//    }
//
//    //todo WTF refiew
//    template<class Writer>
//    void BinWrite(Writer &writer) const {
//        base::index_.serialize(writer);
//        size_t sz = base::data_.size();
//        writer.write((char*)&sz, sizeof(sz));
//        for (size_t i = 0; i < sz; ++i)
//            writer.write((char*)&(base::data_[i].count_), sizeof(base::data_[0].count_));
//        traits::raw_serialize(writer, base::kmers);
//    }
//
//    //todo WTF refiew
//    template<class Reader>
//    void BinRead(Reader &reader, const std::string &FileName) {
//        base::clear();
//        base::index_.deserialize(reader);
//        size_t sz = 0;
//        reader.read((char*)&sz, sizeof(sz));
//        base::data_.resize(sz);
//        for (size_t i = 0; i < sz; ++i)
//            reader.read((char*)&(base::data_[i].count_), sizeof(base::data_[0].count_));
//        base::kmers = traits::raw_deserialize(reader, FileName);
//    }
//
//    void PutInIndex(const KMer &kmer, IdType id, int offset) {
//        size_t idx = base::seq_idx(kmer);
//        if (contains(idx, kmer)) {
//            EdgeInfo<IdType> &entry = base::operator[](idx);
//            entry.edgeId_ = id;
//            entry.offset_ = offset;
//        }
//    }
//};

template<class ValueType, class Seq = runtime_k::RtSeq,
    class traits = kmer_index_traits<Seq> >
class DeBruijnKMerIndex {
 public:
  typedef Seq                      KMer;
  typedef KMerIndex<KMer, traits>  KMerIndexT;

 protected:
  typedef ValueType                       KMerIndexValueType;
  typedef std::vector<KMerIndexValueType> KMerIndexStorageType;
  typedef boost::bimap<KMer, size_t>      KMerPushBackIndexType;

  unsigned K_;
  std::string workdir_;
  KMerIndexT index_;
  KMerIndexStorageType data_;
  typename traits::RawKMerStorage *kmers;
  KMerPushBackIndexType push_back_index_;
  KMerIndexStorageType push_back_buffer_;


 public:
  typedef typename KMerIndexStorageType::iterator value_iterator;
  typedef typename KMerIndexStorageType::const_iterator const_value_iterator;
  typedef typename traits::RawKMerStorage::iterator kmer_iterator;
  typedef typename traits::RawKMerStorage::const_iterator const_kmer_iterator;
  typedef size_t KMerIdx;
  static const size_t InvalidKMerIdx = SIZE_MAX;

  DeBruijnKMerIndex(unsigned K, const std::string &workdir)
      : K_(K), index_(K), kmers(NULL) {
    workdir_ = path::make_temp_dir(workdir, "kmeridx");
  }
  ~DeBruijnKMerIndex() {
    delete kmers;
    path::remove_dir(workdir_);
  }

  void clear() {
    index_.clear();
    data_.clear();
    KMerIndexStorageType().swap(data_);
    push_back_index_.clear();
    push_back_buffer_.clear();
    delete kmers;
    kmers = NULL;
  }

  unsigned K() const { return K_; }

  const KMerIndexValueType &operator[](KMerIdx idx) const {
    if (idx < data_.size())
      return data_[idx];

    return push_back_buffer_[idx - data_.size()];
  }
  KMerIndexValueType &operator[](const KMer &s) {
    return operator[](index_.seq_idx(s));
  }
  const KMerIndexValueType &operator[](const KMer &s) const {
    return operator[](index_.seq_idx(s));
  }
  KMerIndexValueType &operator[](KMerIdx idx) {
    if (idx < data_.size())
      return data_[idx];

    return push_back_buffer_[idx - data_.size()];
  }
  KMerIdx seq_idx(const KMer &s) const {
    KMerIdx idx = index_.seq_idx(s);

    // First, check whether we're insert index itself.
    if (contains(idx, s, /* check push back */ false))
      return idx;

    // Maybe we're inside push_back buffer then?
    auto it = push_back_index_.left.find(s);
    if (it != push_back_index_.left.end())
      return data_.size() + it->second;

    return InvalidKMerIdx;
  }

  size_t size() const { return data_.size() + push_back_buffer_.size(); }

  value_iterator value_begin() {
    return data_.begin();
  }
  const_value_iterator value_begin() const {
    return data_.begin();
  }
  const_value_iterator value_cbegin() const {
    return data_.cbegin();
  }
  value_iterator value_end() {
    return data_.end();
  }
  const_value_iterator value_end() const {
    return data_.end();
  }
  const_value_iterator value_cend() const {
    return data_.cend();
  }

  kmer_iterator kmer_begin() {
    return kmers->begin();
  }
  const_kmer_iterator kmer_begin() const {
    return kmers->cbegin();
  }
  kmer_iterator kmer_end() {
    return kmers->end();
  }
  const_kmer_iterator kmer_end() const {
    return kmers->cend();
  }

  KMerIdx kmer_idx_begin() const {
    return 0;
  }

  KMerIdx kmer_idx_end() const {
    return data_.size();
  }

  bool contains(const KMer &k) const {
    KMerIdx idx = seq_idx(k);

    return idx != InvalidKMerIdx;
  }

  bool contains(KMerIdx idx) const {
    return idx < size();
  }

  size_t insert(const KMer &s, const KMerIndexValueType &value) {
    size_t idx = push_back_buffer_.size();
    push_back_index_.insert(typename KMerPushBackIndexType::value_type(s, idx));
    push_back_buffer_.push_back(value);

    return idx;
  }

  KMer kmer(KMerIdx idx) const {
    VERIFY(contains(idx));

    if (idx < data_.size()) {
      auto it = kmers->begin() + idx;
      return (typename traits::raw_create()(K_, *it));
    }

    idx -= data_.size();
    return push_back_index_.right.find(idx)->second;
  }

  template<class Writer>
  void BinWrite(Writer &writer) const {
    index_.serialize(writer);
    size_t sz = data_.size();
    writer.write((char*)&sz, sizeof(sz));
    writer.write((char*)&data_[0], sz * sizeof(data_[0]));
    sz = push_back_buffer_.size();
    writer.write((char*)&sz, sizeof(sz));
    writer.write((char*)&push_back_buffer_[0], sz * sizeof(push_back_buffer_[0]));
    for (auto it = push_back_index_.left.begin(), e = push_back_index_.left.end(); it != e; ++it) {
      size_t idx = it->second;
      KMer::BinWrite(writer, it->first);
      writer.write((char*)&idx, sizeof(idx));
      sz -= 0;
    }
    VERIFY(sz == 0);
    traits::raw_serialize(writer, kmers);
  }

  template<class Reader>
  void BinRead(Reader &reader, const std::string &FileName) {
    clear();
    index_.deserialize(reader);
    size_t sz = 0;
    reader.read((char*)&sz, sizeof(sz));
    data_.resize(sz);
    reader.read((char*)&data_[0], sz * sizeof(data_[0]));
    reader.read((char*)&sz, sizeof(sz));
    push_back_buffer_.resize(sz);
    reader.read((char*)&push_back_buffer_[0], sz * sizeof(push_back_buffer_[0]));
    for (size_t i = 0; i < sz; ++i) {
      KMer s(K_);
      size_t idx;

      s.BinRead(reader);
      reader.read((char*)&idx, sizeof(idx));

      push_back_index_.insert(typename KMerPushBackIndexType::value_type(s, idx));
    }

    kmers = traits::raw_deserialize(reader, FileName);
  }

  const std::string &workdir() const {
    return workdir_;
  }

  friend class DeBruijnKMerIndexBuilder<Seq>;
 protected:
  bool contains(KMerIdx idx, const KMer &k, bool check_push_back = true) const {
    // Sanity check
    if (idx == InvalidKMerIdx || idx >= size())
      return false;

    if (idx < data_.size()) {
      auto it = kmers->begin() + idx;
      return (typename traits::raw_equal_to()(k, *it));
    }

    if (check_push_back) {
      auto it = push_back_index_.right.find(idx - data_.size());
      return (it != push_back_index_.right.end() &&
              it -> second == k);
    }

    return false;
  }

  size_t raw_seq_idx(const typename KMerIndexT::KMerRawReference s) const {
    return index_.raw_seq_idx(s);
  }
};


template<class IdType, class Seq = runtime_k::RtSeq,
    class traits = kmer_index_traits<Seq> >
class DeBruijnEdgeIndex : public DeBruijnKMerIndex<EdgeInfo<IdType>, Seq, traits> {
  typedef DeBruijnKMerIndex<EdgeInfo<IdType>, Seq, traits> base;
 public:
  typedef Seq                      KMer;
  typedef KMerIndex<KMer, traits>  KMerIndexT;

  DeBruijnEdgeIndex(unsigned K, const std::string &workdir)
      : base(K, workdir) {}

  ~DeBruijnEdgeIndex() {}

  bool ContainsInIndex(typename base::KMerIdx idx) const {
    if (idx == base::InvalidKMerIdx)
      return false;

    const typename base::KMerIndexValueType &entry = base::operator[](idx);
    return (entry.offset_ != -1);
  }

  bool ContainsInIndex(const KMer& kmer) const {
    typename base::KMerIdx idx = base::seq_idx(kmer);

    // Early exit if kmer has not been seen at all
    if (idx == base::InvalidKMerIdx)
      return false;

    // Otherwise, check, whether it's attached to any edge
    const typename base::KMerIndexValueType &entry = base::operator[](idx);
    return (entry.offset_ != -1);
  }

  /**
   * Number of edges coming into param edge's end
   */
  unsigned RivalEdgeCount(const KMer &kmer) const {
    KMer kmer2 = kmer << 'A';
    unsigned res = 0;
    for (char c = 0; c < 4; ++c)
      if (base::contains(kmer2 >> c))
        res += 1;

    return res;
  }

  unsigned RivalEdgeCount(typename base::KMerIdx idx) const {
    KMer kmer2 = kmer(idx) << 'A';
    unsigned res = 0;
    for (char c = 0; c < 4; ++c)
      if (base::contains(kmer2 >> c))
        res += 1;

    return res;
  }

  /**
   * Number of edges going out of the param edge's end
   */
  unsigned NextEdgeCount(const KMer &kmer) const {
    unsigned res = 0;
    for (char c = 0; c < 4; ++c)
      if (base::contains(kmer << c))
        res += 1;

    return res;
  }

  unsigned NextEdgeCount(typename base::KMerIdx idx) const {
    KMer kmer = this->kmer(idx);

    unsigned res = 0;
    for (char c = 0; c < 4; ++c)
      if (contains(kmer << c))
        res += 1;

    return res;
  }

  KMer NextEdge(const KMer &kmer) const { // returns any next edge
    for (char c = 0; c < 4; ++c) {
      KMer s = kmer << c;
      typename base::KMerIdx idx = base::seq_idx(s);
      if (base::contains(idx))
        return this->kmer(idx);
    }

    VERIFY_MSG(false, "Couldn't find requested edge!");
    return KMer(base::K());
    // no next edges (we should request one here).
  }

  KMer NextEdge(typename base::KMerIdx idx) const { // returns any next edge
    KMer kmer = this->kmer(idx);

    for (char c = 0; c < 4; ++c) {
      KMer s = kmer << c;
      if (base::contains(s))
        return s;
    }

    VERIFY_MSG(false, "Couldn't find requested edge!");
    return KMer(base::K());
    // no next edges (we should request one here).
  }

  std::pair<IdType, size_t> get(const KMer &kmer) const {
    typename base::KMerIdx idx = seq_idx(kmer);
    VERIFY(contains(idx, kmer));

    const EdgeInfo<IdType> &entry = base::operator[](idx);
    return std::make_pair(entry.edge_id, entry.offset);
  }

  std::pair<IdType, size_t> get(typename base::KMerIdx idx) const {
    const EdgeInfo<IdType> &entry = base::operator[](idx);
    return std::make_pair(entry.edge_id, entry.offset);
  }

  bool DeleteIfEqual(const KMer &kmer, IdType id) {
    typename base::KMerIdx idx = seq_idx(kmer);

    if (!contains(idx, kmer))
      return false;

    EdgeInfo<IdType> &entry = base::operator[](idx);

    if (entry.edge_id == id) {
      entry.offset = -1u;
      return true;
    }

    return false;
  }

  template<class Writer>
  void BinWrite(Writer &writer) const {
    base::index_.serialize(writer);
    size_t sz = base::data_.size();
    writer.write((char*)&sz, sizeof(sz));
    for (size_t i = 0; i < sz; ++i)
      writer.write((char*)&(base::data_[i].count_), sizeof(base::data_[0].count_));
    sz = base::push_back_buffer_.size();
    writer.write((char*)&sz, sizeof(sz));
    for (size_t i = 0; i < sz; ++i)
      writer.write((char*)&(base::push_back_buffer_[i].count_), sizeof(base::push_back_buffer_[0].count_));
    for (auto it = base::push_back_index_.left.begin(), e = base::push_back_index_.left.end(); it != e; ++it) {
      size_t idx = it->second;
      KMer::BinWrite(writer, it->first);
      writer.write((char*)&idx, sizeof(idx));
      sz -= 1;
    }
    VERIFY(sz == 0);
    traits::raw_serialize(writer, base::kmers);
  }

  template<class Reader>
  void BinRead(Reader &reader, const std::string &FileName) {
    base::clear();
    base::index_.deserialize(reader);
    size_t sz = 0;
    reader.read((char*)&sz, sizeof(sz));
    base::data_.resize(sz);
    for (size_t i = 0; i < sz; ++i)
      reader.read((char*)&(base::data_[i].count_), sizeof(base::data_[0].count_));
    reader.read((char*)&sz, sizeof(sz));
    base::push_back_buffer_.resize(sz);
    for (size_t i = 0; i < sz; ++i)
      reader.read((char*)&(base::push_back_buffer_[i].count_), sizeof(base::push_back_buffer_[0].count_));
    for (size_t i = 0; i < sz; ++i) {
      KMer s(base::K_);
      size_t idx;

      s.BinRead(reader);
      reader.read((char*)&idx, sizeof(idx));

      base::push_back_index_.insert(typename base::KMerPushBackIndexType::value_type(s, idx));
    }
    base::kmers = traits::raw_deserialize(reader, FileName);
  }

  void PutInIndex(const KMer &kmer, IdType id, int offset, bool ignore_new_kmer = false) {
    size_t idx = base::seq_idx(kmer);
    if (base::contains(idx, kmer)) {
      EdgeInfo<IdType> &entry = base::operator[](idx);
      entry.edge_id = id;
      entry.offset = offset;
    }
  }
};

//todo remove after checking that it is not needed
//template <>
//class OldDeBruijnEdgeIndexBuilder<runtime_k::RtSeq> :
//            public DeBruijnKMerIndexBuilder<kmer_index_traits<runtime_k::RtSeq> > {
//    typedef DeBruijnKMerIndexBuilder<kmer_index_traits<runtime_k::RtSeq> > base;
//
//    template <class ReadStream, class Graph>
//    size_t FillCoverageFromStream(ReadStream &stream,
//                                  DeBruijnEdgeIndex<Graph, runtime_k::RtSeq> &index) const {
//        unsigned K = index.K();
//        size_t rl = 0;
//
//        while (!stream.eof()) {
//            typename ReadStream::read_type r;
//            stream >> r;
//            rl = std::max(rl, r.size());
//
//            const Sequence &seq = r.sequence();
//            if (seq.size() < K)
//                continue;
//
//            runtime_k::RtSeq kmer = seq.start<runtime_k::RtSeq>(K);
//
//            size_t idx = index.seq_idx(kmer);
//            if (index.contains(idx, kmer)) {
//#   pragma omp atomic
//                index.data_[idx].count_ += 1;
//            }
//            for (size_t j = K; j < seq.size(); ++j) {
//                kmer <<= seq[j];
//                idx = index.seq_idx(kmer);
//                if (index.contains(idx, kmer)) {
//#     pragma omp atomic
//                    index.data_[idx].count_ += 1;
//                }
//            }
//        }
//
//        return rl;
//    }
//
//  public:
//
//    template<class Graph, class Read>
//    size_t ParallelFillCoverage(DeBruijnEdgeIndex<Graph, runtime_k::RtSeq> &index,
//                                io::ReadStreamVector<io::IReader<Read> > &streams,
//                                SingleReadStream* contigs_stream = 0) const {
//        INFO("Collecting k-mer coverage information from reads, this takes a while.");
//
//        unsigned nthreads = streams.size();
//        size_t rl = 0;
//        streams.reset();
//#pragma omp parallel for num_threads(nthreads) shared(rl)
//        for (size_t i = 0; i < nthreads; ++i) {
//            size_t crl = FillCoverageFromStream(streams[i], index);
//
//            // There is no max reduction in C/C++ OpenMP... Only in FORTRAN :(
//#pragma omp flush(rl)
//            if (crl > rl)
//#pragma omp critical
//            {
//                rl = std::max(rl, crl);
//            }
//        }
//
//        // Contigs have zero coverage!
//#if 0
//        if (contigs_stream) {
//            contigs_stream->reset();
//            FillCoverageFromStream(*contigs_stream, index);
//        }
//#endif
//
//#ifndef NDEBUG
//        for (auto idx = index.kmer_idx_begin(), eidx = index.kmer_idx_end();
//             idx != eidx; ++idx) {
//
//            runtime_k::RtSeq k = index.kmer(idx);
//
//            VERIFY(index[k].count_ == index[!k].count_);
//        }
//#endif
//        return rl;
//    }
//
//    template<class Graph, class Read>
//    size_t BuildIndexFromStream(DeBruijnEdgeIndex<Graph, runtime_k::RtSeq> &index,
//                                io::ReadStreamVector<io::IReader<Read> > &streams,
//                                SingleReadStream* contigs_stream = 0) const {
//        base::BuildIndexFromStream(index, streams, contigs_stream);
//
//        // Now use the index to fill the coverage and EdgeId's
//        return ParallelFillCoverage(index, streams, contigs_stream);
//    }
//
////    template<class Graph, class Read>
////    size_t BuildIndexWithCoverageFromGraph(Graph &graph,
////                                           DeBruijnEdgeIndex<Graph, runtime_k::RtSeq> &index,
////                                           io::ReadStreamVector<io::IReader<Read> > &streams,
////                                           SingleReadStream* contigs_stream = 0) const {
////        BuildIndexFromGraph(index, graph);
////
////        // Now use the index to fill the coverage and EdgeId's
////        return ParallelFillCoverage(index, streams, contigs_stream);
////    }
//
//    template <class Graph>
//    void BuildIndexFromGraph(DeBruijnEdgeIndex<Graph, runtime_k::RtSeq> &index,
//                             const Graph &g) const {
//        base::BuildIndexFromGraph(index, g);
//
//        // Now use the index to fill the coverage and EdgeId's
//        INFO("Collecting k-mer coverage information from graph, this takes a while.");
//    }
//
//  protected:
//    DECL_LOGGER("Edge Index Building");
//};

}
