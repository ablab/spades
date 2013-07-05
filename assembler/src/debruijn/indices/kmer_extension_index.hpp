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

//class DeBruijnExtensionIndexBuilder;

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
class DeBruijnExtensionIndex : public DeBruijnKMerIndex<uint8_t, traits> {
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

  public:
    typedef typename base::traits_t traits_t;
    typedef typename base::KMer KMer;
    typedef typename base::KMerIdx KMerIdx;
    typedef InnerDeBruijnTotallyKMerFreeIndexBuilder<DeBruijnExtensionIndex> BuilderT;

    //    typedef KMerIndex<traits>        KMerIndexT;

    typedef MMappedFileRecordArrayIterator<typename KMer::DataType> kmer_iterator;

    kmer_iterator kmer_begin() const {
        return kmer_iterator(this->KMersFilename_, KMer::GetDataSize(base::K()));
    }

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
        return KmerWithHash<KMer>(kmer, seq_idx(kmer));
    }

};

//template <>
//class DeBruijnKMerIndexBuilder<slim_kmer_index_traits<runtime_k::RtSeq>> {
// public:
//  template <class IdType, class Read>
//  std::string BuildIndexFromStream(DeBruijnKMerIndex<IdType, slim_kmer_index_traits<runtime_k::RtSeq>> &index,
//                                   io::ReadStreamVector<io::IReader<Read> > &streams,
//                                   SingleReadStream* contigs_stream = 0) const {
//    DeBruijnReadKMerSplitter<Read> splitter(index.workdir(),
//                                            index.K(),
//                                            streams, contigs_stream);
//    KMerDiskCounter<runtime_k::RtSeq> counter(index.workdir(), splitter);
//    KMerIndexBuilder<typename DeBruijnKMerIndex<IdType, slim_kmer_index_traits<runtime_k::RtSeq>>::KMerIndexT> builder(index.workdir(), 16, streams.size());
//
//    size_t sz = builder.BuildIndex(index.index_, counter, /* save final */ true);
//    index.data_.resize(sz);
//    index.kmers = NULL;
//
//    return counter.GetFinalKMersFname();
//  }
//
// protected:
//  DECL_LOGGER("K-mer Index Building");
//};

//class DeBruijnExtensionIndexBuilder :
//            public DeBruijnKMerIndexBuilder<DeBruijnExtensionIndex<>::kmer_index_traits> {
//    typedef DeBruijnKMerIndexBuilder<DeBruijnExtensionIndex<>::kmer_index_traits> base;
//
//    template <class ReadStream>
//    size_t FillExtensionsFromStream(ReadStream &stream,
//                                    DeBruijnExtensionIndex<> &index) const {
//        unsigned K = index.K();
//        size_t rl = 0;
//
//        while (!stream.eof()) {
//            typename ReadStream::read_type r;
//            stream >> r;
//            rl = std::max(rl, r.size());
//
//            const Sequence &seq = r.sequence();
//            if (seq.size() < K + 1)
//                continue;
//
//            runtime_k::RtSeq kmer = seq.start<runtime_k::RtSeq>(K);
//            for (size_t j = K; j < seq.size(); ++j) {
//                char nnucl = seq[j], pnucl = kmer[0];
//                index.AddOutgoing(kmer, nnucl);
//                kmer <<= nnucl;
//                index.AddIncoming(kmer, pnucl);
//            }
//        }
//
//        return rl;
//    }
//
//  public:
//    template <class Read>
//    size_t BuildIndexFromStream(DeBruijnExtensionIndex<> &index,
//                                io::ReadStreamVector<io::IReader<Read> > &streams,
//                                SingleReadStream* contigs_stream = 0) const {
//        unsigned nthreads = streams.size();
//
//        index.KMersFilename_ = base::BuildIndexFromStream(index, streams, contigs_stream);
//
//        // Now use the index to fill the coverage and EdgeId's
//        INFO("Building k-mer extensions from reads, this takes a while.");
//
//        size_t rl = 0;
//        streams.reset();
//# pragma omp parallel for num_threads(nthreads) shared(rl)
//        for (size_t i = 0; i < nthreads; ++i) {
//            size_t crl = FillExtensionsFromStream(streams[i], index);
//
//            // There is no max reduction in C/C++ OpenMP... Only in FORTRAN :(
//#   pragma omp flush(rl)
//            if (crl > rl)
//#     pragma omp critical
//            {
//                rl = std::max(rl, crl);
//            }
//        }
//
//        if (contigs_stream) {
//            contigs_stream->reset();
//            FillExtensionsFromStream(*contigs_stream, index);
//        }
//        INFO("Building k-mer extensions from reads finished.");
//
//        return rl;
//    }
//
//  protected:
//    DECL_LOGGER("Extension Index Building");
//};

template<class Builder>
class DeBruijnExtensionIndexBuilder : public Builder {
    typedef Builder base;
 public:
    typedef typename Builder::IndexT IndexT;

    template <class ReadStream>
    size_t FillExtensionsFromStream(ReadStream &stream,
                                    IndexT &index) const {
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

  public:
    template <class Streams>
    size_t BuildExtensionIndexFromStream(IndexT &index,
                                Streams &streams,
                                SingleReadStream* contigs_stream = 0) const {
        unsigned nthreads = streams.size();

        BuildIndexFromStream(index, streams, contigs_stream);

        // Now use the index to fill the coverage and EdgeId's
        INFO("Building k-mer extensions from reads, this takes a while.");

        size_t rl = 0;
        streams.reset();
# pragma omp parallel for num_threads(nthreads) shared(rl)
        for (size_t i = 0; i < nthreads; ++i) {
            size_t crl = FillExtensionsFromStream(streams[i], index);

            // There is no max reduction in C/C++ OpenMP... Only in FORTRAN :(
#   pragma omp flush(rl)
            if (crl > rl)
#     pragma omp critical
            {
                rl = std::max(rl, crl);
            }
        }

        if (contigs_stream) {
            contigs_stream->reset();
            FillExtensionsFromStream(*contigs_stream, index);
        }
        INFO("Building k-mer extensions from reads finished.");

        return rl;
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
    typedef typename IndexT::BuilderT BuilderT;
    typedef DeBruijnKMerIndexBuilder<Kmer, BuilderT> DeBruijnKMerIndexBuilderT;
    typedef DeBruijnExtensionIndexBuilder<DeBruijnKMerIndexBuilderT> DeBruijnExtensionIndexBuilderT;
};

}
