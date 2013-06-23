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
class DeBruijnExtensionIndexBuilder;

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

template<class Seq = runtime_k::RtSeq, class traits = slim_kmer_index_traits<Seq> >
class DeBruijnExtensionIndex : public DeBruijnKMerIndex<uint8_t, traits> {
  private:
    typedef DeBruijnKMerIndex<uint8_t, traits> base;

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

    std::string KMersFilename_;

  public:
    typedef typename traits::SeqType KMer;
    typedef traits                   kmer_index_traits;
    typedef KMerIndex<traits>        KMerIndexT;

    typedef MMappedFileRecordArrayIterator<typename KMer::DataType> kmer_iterator;

    kmer_iterator kmer_begin() const {
        return kmer_iterator(KMersFilename_, KMer::GetDataSize(base::K()));
    }

    DeBruijnExtensionIndex(unsigned K, const std::string &workdir)
            : base(K, workdir), KMersFilename_("") {}

    void AddOutgoing(size_t idx, char nnucl) {
        unsigned nmask = (1 << nnucl);
        if (!(this->operator [](idx) & nmask)) {
#           pragma omp atomic
            this->operator [](idx) |= nmask;
        }
    }

    void AddIncoming(size_t idx, char pnucl) {
        unsigned pmask = (1 << (pnucl + 4));

        if (!(this->operator [](idx) & pmask)) {
#           pragma omp atomic
            this->operator [](idx) |= pmask;
        }
    }

    void DeleteOutgoing(size_t idx, char nnucl) {
        unsigned nmask = (1 << nnucl);
        if (this->operator [](idx) & nmask) {
#           pragma omp atomic
            this->operator [](idx) &= ~nmask;
        }
    }

    void DeleteIncoming(size_t idx, char pnucl) {
        unsigned pmask = (1 << (pnucl + 4));

        if (this->operator [](idx) & pmask) {
#         pragma omp atomic
            this->operator [](idx) &= ~pmask;
        }
    }

    void IsolateVertex(size_t idx) {
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

    bool CheckOutgoing(size_t idx, size_t nucl) {
        return (this->operator [](idx)) & (1 << nucl);
    }

    bool CheckIncoming(size_t idx, size_t nucl) {
        return (this->operator [](idx)) & (16 << nucl);
    }

    bool IsDeadEnd(size_t idx) const {
        return !(this->operator [](idx) & 15);
    }

    bool IsDeadStart(size_t idx) const {
        return !(this->operator [](idx) >> 4);
    }

    bool CheckUniqueOutgoing(size_t idx) const {
        return CheckUnique(this->operator [](idx) & 15);
    }

    char GetUniqueOutgoing(size_t idx) const {
        return GetUnique(this->operator [](idx) & 15);
    }

    bool CheckUniqueIncoming(size_t idx) const {
        return CheckUnique(this->operator [](idx) >> 4);
    }

    char GetUniqueIncoming(size_t idx) const {
        return GetUnique(this->operator [](idx) >> 4);
    }

    size_t OutgoingEdgeCount(size_t idx) const {
        return Count(this->operator [](idx) & 15);
    }

    size_t IncomingEdgeCount(size_t idx) const {
        return Count(this->operator [](idx) >> 4);
    }

    ~DeBruijnExtensionIndex() {}


    struct KmerWithHash {
        KMer kmer;
        size_t idx;

        KmerWithHash(KMer _kmer,
                     const DeBruijnExtensionIndex<Seq, traits> &index) :
                kmer(_kmer), idx(index.seq_idx(kmer)) { }
    };

    KmerWithHash CreateKmerWithHash(KMer kmer) const {
        return KmerWithHash(kmer, *this);
    }

    friend class DeBruijnExtensionIndexBuilder<Seq>;
};

template <>
class DeBruijnKMerIndexBuilder<slim_kmer_index_traits<runtime_k::RtSeq>> {
 public:
  template <class IdType, class Read>
  std::string BuildIndexFromStream(DeBruijnKMerIndex<IdType, slim_kmer_index_traits<runtime_k::RtSeq>> &index,
                                   io::ReadStreamVector<io::IReader<Read> > &streams,
                                   SingleReadStream* contigs_stream = 0) const {
    DeBruijnReadKMerSplitter<Read> splitter(index.workdir(),
                                            index.K(),
                                            streams, contigs_stream);
    KMerDiskCounter<runtime_k::RtSeq> counter(index.workdir(), splitter);
    KMerIndexBuilder<typename DeBruijnKMerIndex<IdType, slim_kmer_index_traits<runtime_k::RtSeq>>::KMerIndexT> builder(index.workdir(), 16, streams.size());

    size_t sz = builder.BuildIndex(index.index_, counter, /* save final */ true);
    index.data_.resize(sz);
    index.kmers = NULL;

    return counter.GetFinalKMersFname();
  }

 protected:
  DECL_LOGGER("K-mer Index Building");
};

template <class Seq>
class DeBruijnExtensionIndexBuilder : public DeBruijnKMerIndexBuilder<typename DeBruijnExtensionIndex<Seq>::kmer_index_traits> {
  public:
    template <class Read>
    size_t BuildIndexFromStream(DeBruijnExtensionIndex<Seq> &index,
                                io::ReadStreamVector<io::IReader<Read> > &streams,
                                SingleReadStream* contigs_stream = 0) const;

  protected:
    DECL_LOGGER("Extension Index Building");
};

template <>
class DeBruijnExtensionIndexBuilder<runtime_k::RtSeq> :
            public DeBruijnKMerIndexBuilder<DeBruijnExtensionIndex<runtime_k::RtSeq>::kmer_index_traits> {
    typedef DeBruijnKMerIndexBuilder<DeBruijnExtensionIndex<runtime_k::RtSeq>::kmer_index_traits> base;

    template <class ReadStream>
    size_t FillExtensionsFromStream(ReadStream &stream,
                                    DeBruijnExtensionIndex<runtime_k::RtSeq> &index) const {
        unsigned K = index.K();
        size_t rl = 0;

        while (!stream.eof()) {
            typename ReadStream::read_type r;
            stream >> r;
            rl = std::max(rl, r.size());

            const Sequence &seq = r.sequence();
            if (seq.size() < K + 1)
                continue;

            runtime_k::RtSeq kmer = seq.start<runtime_k::RtSeq>(K);
            for (size_t j = K; j < seq.size(); ++j) {
                char nnucl = seq[j], pnucl = kmer[0];
                index.AddOutgoing(kmer, nnucl);
                kmer <<= nnucl;
                index.AddIncoming(kmer, pnucl);
            }
        }

        return rl;
    }

    void FillExtensionsFromIndex(std::string KPlusOneMersFilename,
                                 DeBruijnExtensionIndex<runtime_k::RtSeq> &index) const {
        unsigned KPlusOne = index.K() + 1;

        DeBruijnExtensionIndex<runtime_k::RtSeq>::kmer_iterator it(KPlusOneMersFilename,
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
    template <class Read>
    size_t BuildIndexFromStream(DeBruijnExtensionIndex<runtime_k::RtSeq> &index,
                                io::ReadStreamVector<io::IReader<Read> > &streams,
                                SingleReadStream* contigs_stream = 0) const {
        // First, build a k+1-mer index
        DeBruijnReadKMerSplitter<Read> splitter(index.workdir(),
                                                index.K() + 1,
                                                streams, contigs_stream);
        KMerDiskCounter<runtime_k::RtSeq> counter(index.workdir(), splitter);
        counter.CountAll(16, streams.size());

        // Now, count unique k-mers from k+1-mers
        DeBruijnKMerKMerSplitter splitter2(index.workdir(),
                                           index.K(),
                                           counter.GetFinalKMersFname(), index.K() + 1);
        KMerDiskCounter<runtime_k::RtSeq> counter2(index.workdir(), splitter2);

        KMerIndexBuilder<typename DeBruijnKMerIndex<uint8_t, slim_kmer_index_traits<runtime_k::RtSeq>>::KMerIndexT> builder(index.workdir(), 16, 1);

        size_t sz = builder.BuildIndex(index.index_, counter2, /* save final */ true);
        index.data_.resize(sz);
        index.kmers = NULL;
        index.KMersFilename_ = counter2.GetFinalKMersFname();

        // Build the kmer extensions
        INFO("Building k-mer extensions from k+1-mers");
        FillExtensionsFromIndex(counter.GetFinalKMersFname(), index);
        INFO("Building k-mer extensions from k+1-mers finished.");

        return splitter.read_length();
    }

  protected:
    DECL_LOGGER("Extension Index Building");
};

}
