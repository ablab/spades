#pragma once
/*
 * kmer_extension_index.hpp
 *
 *  Created on: May 24, 2013
 *      Author: anton
 */
#include "debruijn_kmer_index.hpp"
#include "kmer_splitters.hpp"

namespace debruijn_graph {

template<class Seq>
struct slim_kmer_index_traits : public kmer_index_traits<Seq> {
  typedef kmer_index_traits<Seq> __super;

  typedef MMappedRecordReader<typename Seq::DataType> FinalKMerStorage;

  template<class Writer>
  static void raw_serialize(Writer&, typename __super::RawKMerStorage*) {
    VERIFY(false && "Cannot save extension index");
  }

  template<class Reader>
  static typename __super::RawKMerStorage *raw_deserialize(Reader&, const std::string &) {
    VERIFY(false && "Cannot load extension index");
    return NULL;
  }

};

template<class Seq>
struct KmerWithHash {
    typedef size_t KMerIdx;
    Seq kmer;
    KMerIdx idx;

    KmerWithHash(Seq kmer_, KMerIdx idx_) :
            kmer(kmer_), idx(idx_) { }
};

template<class traits = slim_kmer_index_traits<runtime_k::RtSeq>>
class DeBruijnExtensionIndex : public KmerFreeIndex<uint8_t, traits> {
    typedef KmerFreeIndex<uint8_t, traits> base;

    bool CheckUnique(uint8_t mask) const {
        static bool unique[] = {0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0};
        return unique[mask];
    }

    char GetUnique(uint8_t mask) const {
        static char next[] = {-1, 0, 1, -1, 2, -1, -1, -1, 3, -1, -1, -1, -1, -1, -1, -1};
        return next[mask];
    }

    size_t Count(uint8_t mask) const {
        static char count[] = {0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4};
        return count[mask];
    }

  public:
    typedef typename base::traits_t traits_t;
    typedef typename base::KeyType KMer;
    typedef typename base::IdxType KMerIdx;

    DeBruijnExtensionIndex(unsigned K, const std::string &workdir)
            : base(K, workdir) {}

    void AddOutgoing(KMerIdx idx, char nnucl) {
        unsigned nmask = (1 << nnucl);
        if (!(this->operator [](idx) & nmask)) {
#           pragma omp atomic
            this->operator [](idx) |= nmask;
        }
    }

    void AddIncoming(KMerIdx idx, char pnucl) {
        unsigned pmask = (1 << (pnucl + 4));

        if (!(this->operator [](idx) & pmask)) {
#           pragma omp atomic
            this->operator [](idx) |= pmask;
        }
    }

    void DeleteOutgoing(KMerIdx idx, char nnucl) {
        unsigned nmask = (1 << nnucl);
        if (this->operator [](idx) & nmask) {
#           pragma omp atomic
            this->operator [](idx) &= ~nmask;
        }
    }

    void DeleteIncoming(KMerIdx idx, char pnucl) {
        unsigned pmask = (1 << (pnucl + 4));

        if (this->operator [](idx) & pmask) {
#         pragma omp atomic
            this->operator [](idx) &= ~pmask;
        }
    }

    void IsolateVertex(KMerIdx idx) {
        this->operator [](idx) = 0;
    }

    void AddOutgoing(const runtime_k::RtSeq &kmer, char nnucl) {
        AddOutgoing(this->seq_idx(kmer), nnucl);
    }

    void AddIncoming(const runtime_k::RtSeq &kmer, char pnucl) {
        AddIncoming(this->seq_idx(kmer), pnucl);
    }

    void DeleteOutgoing(const runtime_k::RtSeq &kmer, char nnucl) {
        DeleteOutgoing(this->seq_idx(kmer), nnucl);
    }

    void DeleteIncoming(const runtime_k::RtSeq &kmer, char pnucl) {
        DeleteIncoming(this->seq_idx(kmer), pnucl);
    }

    bool CheckOutgoing(KMerIdx idx, size_t nucl) {
        return (this->operator [](idx)) & (1 << nucl);
    }

    bool CheckIncoming(KMerIdx idx, size_t nucl) {
        return (this->operator [](idx)) & (16 << nucl);
    }

    bool IsDeadEnd(KMerIdx idx) const {
        return !(this->operator [](idx) & 15);
    }

    bool IsDeadStart(KMerIdx idx) const {
        return !(this->operator [](idx) >> 4);
    }

    bool CheckUniqueOutgoing(KMerIdx idx) const {
        return CheckUnique(this->operator [](idx) & 15);
    }

    char GetUniqueOutgoing(KMerIdx idx) const {
        return GetUnique(this->operator [](idx) & 15);
    }

    bool CheckUniqueIncoming(KMerIdx idx) const {
        return CheckUnique(this->operator [](idx) >> 4);
    }

    char GetUniqueIncoming(KMerIdx idx) const {
        return GetUnique(this->operator [](idx) >> 4);
    }

    size_t OutgoingEdgeCount(KMerIdx idx) const {
        return Count(this->operator [](idx) & 15);
    }

    size_t IncomingEdgeCount(KMerIdx idx) const {
        return Count(this->operator [](idx) >> 4);
    }

    ~DeBruijnExtensionIndex() {}

    KmerWithHash<KMer> CreateKmerWithHash(KMer kmer) const {
        return KmerWithHash<KMer>(kmer, this->seq_idx(kmer));
    }

};

template<class Builder>
class DeBruijnExtensionIndexBuilder : public Builder {
    typedef Builder base;
 public:
    typedef typename Builder::IndexT IndexT;

    template <class ReadStream>
    size_t FillExtensionsFromStream(ReadStream &stream,
                                    IndexT &index) const {
        unsigned k = index.k();
        size_t rl = 0;

        while (!stream.eof()) {
            typename ReadStream::read_type r;
            stream >> r;
            rl = std::max(rl, r.size());

            const Sequence &seq = r.sequence();
            if (seq.size() < k + 1)
                continue;

            runtime_k::RtSeq kmer = seq.start<runtime_k::RtSeq>(k);
            for (size_t j = k; j < seq.size(); ++j) {
                char nnucl = seq[j], pnucl = kmer[0];
                index.AddOutgoing(kmer, nnucl);
                kmer <<= nnucl;
                index.AddIncoming(kmer, pnucl);
            }
        }

        return rl;
    }

    void FillExtensionsFromIndex(const std::string &KPlusOneMersFilename,
                                 IndexT &index) const {
        unsigned KPlusOne = index.k() + 1;

        typename IndexT::kmer_iterator it(KPlusOneMersFilename,
                                          runtime_k::RtSeq::GetDataSize(KPlusOne));
        for (; it.good(); ++it) {
            runtime_k::RtSeq kpomer(KPlusOne, *it);

            char pnucl = kpomer[0], nnucl = kpomer[KPlusOne - 1];
            index.AddOutgoing(runtime_k::RtSeq(KPlusOne - 1, kpomer), nnucl);
            // FIXME: This is extremely ugly. Needs to add start / end methods to extract first / last N symbols...
            index.AddIncoming(runtime_k::RtSeq(KPlusOne - 1, kpomer << 0), pnucl);
        }
    }


  public:
    template <class Streams>
    size_t BuildExtensionIndexFromStream(IndexT &index,
                                         Streams &streams,
                                         SingleReadStream* contigs_stream = 0) const {
        unsigned nthreads = streams.size();

        // First, build a k+1-mer index
        DeBruijnReadKMerSplitter<typename Streams::ReaderType::read_type>
                splitter(index.workdir(), index.k() + 1, 0xDEADBEEF,
                         streams, contigs_stream);
        KMerDiskCounter<runtime_k::RtSeq> counter(index.workdir(), splitter);
        counter.CountAll(nthreads, nthreads, /* merge */ false);

        // Now, count unique k-mers from k+1-mers
        DeBruijnKMerKMerSplitter splitter2(index.workdir(),
                                           index.k(), index.k() + 1);
        for (size_t i = 0; i < nthreads; ++i)
          splitter2.AddKMers(counter.GetMergedKMersFname(i));
        KMerDiskCounter<runtime_k::RtSeq> counter2(index.workdir(), splitter2);

        index.BuildIndex(counter2, 16, nthreads);

        // Build the kmer extensions
        INFO("Building k-mer extensions from k+1-mers");
#       pragma omp parallel for num_threads(nthreads)
        for (size_t i = 0; i < nthreads; ++i)
          FillExtensionsFromIndex(counter.GetMergedKMersFname(i), index);
        INFO("Building k-mer extensions from k+1-mers finished.");

        return splitter.read_length();
    }

  private:
    DECL_LOGGER("DeBruijnExtensionIndexBuilder");
};

template<class Index>
struct ExtensionIndexHelper {
    typedef Index IndexT;
    typedef typename IndexT::traits_t traits_t;
    typedef typename IndexT::KMer Kmer;
    typedef typename IndexT::KMerIdx KMerIdx;
    typedef DeBruijnStreamKMerIndexBuilder<Kmer, IndexT> DeBruijnStreamKMerIndexBuilderT;
    typedef DeBruijnExtensionIndexBuilder<DeBruijnStreamKMerIndexBuilderT> DeBruijnExtensionIndexBuilderT;
};

}
