#pragma once
/*
 * kmer_extension_index.hpp
 *
 *  Created on: May 24, 2013
 *      Author: anton
 */
#include "perfect_hash_map.hpp"
#include "kmer_splitters.hpp"
#include "simple_tools.hpp"
#include "storing_traits.hpp"

namespace debruijn_graph {

inline uint8_t invert_byte_slow(uint8_t a) {
    size_t res = 0;
    for(size_t i = 0; i < 8; i++) {
        res <<= 1;
        res += a & 1;
        a = uint8_t(a >> 1);
    }
    return uint8_t(res);
}

inline vector<uint8_t> count_invert_byte() {
    vector<uint8_t> result;
    for(size_t a = 0; a < 256; a++) {
        result.push_back(invert_byte_slow((uint8_t)a));
    }
    return result;
}

inline uint8_t invert_byte(uint8_t a) {
    static vector<uint8_t> precalc = count_invert_byte();
    return precalc[a];
}

class InOutMask {
private:
	uint8_t mask_;

    bool CheckUnique(uint8_t mask) const {
        static bool unique[] =
                { 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 };
        return unique[mask];
    }

    char GetUnique(uint8_t mask) const {
        static char next[] = { -1, 0, 1, -1, 2, -1, -1, -1, 3, -1, -1, -1, -1,
                -1, -1, -1 };
        VERIFY(next[mask] != -1)
        return next[mask];
    }

    size_t Count(uint8_t mask) const {
        static char count[] = { 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4 };
        return count[mask];
    }


    char inv_position(char nucl, bool as_is) const {
        if(as_is)
            return nucl;
        else
            return char(7 - nucl);
    }

public:
	explicit InOutMask(uint8_t mask = 0) : mask_(mask){
	}

	uint8_t get_mask() const {
	    return mask_;
	}

	template<class Key>
	InOutMask conjugate(const Key & /*k*/) const {
		return InOutMask(invert_byte(mask_));
	}

    void AddOutgoing(char nnucl, bool as_is) {
        unsigned nmask = (unsigned) (1 << inv_position(nnucl, as_is));
        if (!(mask_ & nmask)) {
#           pragma omp atomic
            mask_ |= (unsigned char) nmask;
        }
    }

    void AddIncoming(char pnucl, bool as_is) {
        unsigned pmask = (unsigned) (1 << inv_position(char(pnucl + 4), as_is));
        if (!(mask_ & pmask)) {
#           pragma omp atomic
        	mask_|= (unsigned char) pmask;
        }
    }

    void DeleteOutgoing(char nnucl, bool as_is) {
        unsigned nmask = (1 << inv_position(nnucl, as_is));
        if (mask_ & nmask) {
#           pragma omp atomic
        	mask_ &= (unsigned char) ~nmask;
        }
    }

    void DeleteIncoming(char pnucl, bool as_is) {
        unsigned pmask = (1 << inv_position(char(pnucl + 4), as_is));
        if (mask_ & pmask) {
#           pragma omp atomic
        	mask_ &= (unsigned char) ~pmask;
        }
    }

    void IsolateVertex() {
    	mask_ = 0;
    }

    bool CheckOutgoing(char nucl) const {
        return mask_ & (1 << nucl);
    }

    bool CheckIncoming(char nucl) const {
        return mask_ & (1 << (4 + nucl));
    }

    bool IsDeadEnd() const {
        return !(mask_ & 15);
    }

    bool IsDeadStart() const {
        return !(mask_ >> 4);
    }

    bool CheckUniqueOutgoing() const {
        return CheckUnique(mask_ & 15);
    }

    bool CheckUniqueIncoming() const {
        return CheckUnique(uint8_t(mask_ >> 4));
    }

    char GetUniqueOutgoing() const {
        return GetUnique(mask_ & 15);
    }

    char GetUniqueIncoming() const {
        return GetUnique(uint8_t(mask_ >> 4));
    }

    size_t OutgoingEdgeCount() const {
        return Count(mask_ & 15);
    }

    size_t IncomingEdgeCount() const {
        return Count(uint8_t(mask_ >> 4));
    }
};

template<class Stream>
Stream &operator<<(Stream& stream, const InOutMask &mask) {
    return stream << std::bitset<8>(mask.get_mask());
}

template<class Seq>
struct slim_kmer_index_traits : public kmer_index_traits<Seq> {
    typedef kmer_index_traits<Seq> __super;

    typedef MMappedRecordReader<typename Seq::DataType> FinalKMerStorage;

    template<class Writer>
    static void raw_serialize(Writer&, typename __super::RawKMerStorage*) {
        VERIFY(false && "Cannot save extension index");
    }

    template<class Reader>
    static typename __super::RawKMerStorage *raw_deserialize(
            Reader&, const std::string &) {
        VERIFY(false && "Cannot load extension index");
        return NULL;
    }

};

template<typename KeyWithHash>
struct AbstractDeEdge {
    KeyWithHash start;
    KeyWithHash end;
    AbstractDeEdge(KeyWithHash _start, KeyWithHash _end) : start(_start), end(_end) {
    }

    AbstractDeEdge<KeyWithHash> &operator=(const AbstractDeEdge<KeyWithHash> &that) {
    	this->start = that.start;
    	this->end = that.end;
    	return *this;
    }

    bool operator==(const AbstractDeEdge &other) {
        return start.idx() == other.start.idx() && end.idx() == other.end.idx();
    }

    bool operator!=(const AbstractDeEdge &other) {
        return !(*this == other);
    }
};

template<class stream, class KWH>
stream &operator<<(stream &s, const AbstractDeEdge<KWH> de_edge) {
    return s << "DeEdge[" << de_edge.start << ", " << de_edge.end << "]";
}

template<class traits = slim_kmer_index_traits<runtime_k::RtSeq>, class StoringType = DefaultStoring>
class DeBruijnExtensionIndex : public KeyIteratingMap<typename traits::SeqType, InOutMask, traits, StoringType> {
    typedef KeyIteratingMap<typename traits::SeqType, InOutMask, traits, StoringType> base;

public:
    typedef typename base::traits_t traits_t;
    typedef StoringType storing_type;
    typedef typename base::KeyType KMer;
    typedef typename base::IdxType KMerIdx;
    typedef typename base::KeyWithHash KeyWithHash;
    typedef AbstractDeEdge<KeyWithHash> DeEdge;
    using base::ConstructKWH;

    DeBruijnExtensionIndex(unsigned K, const std::string &workdir)
            : base((size_t) K, workdir) {
    }

    void AddOutgoing(const KeyWithHash &kwh, char nucl) {
        TRACE("Add outgoing " << kwh << " " << size_t(nucl) << " " << kwh.is_minimal());
    	this->get_raw_value_reference(kwh).AddOutgoing(nucl, kwh.is_minimal());
    }

    void AddIncoming(const KeyWithHash &kwh, char nucl) {
        TRACE("Add incoming " << kwh << " " << size_t(nucl) << " " << kwh.is_minimal());
    	this->get_raw_value_reference(kwh).AddIncoming(nucl, kwh.is_minimal());
    }

    void DeleteOutgoing(const KeyWithHash &kwh, char nucl) {
        TRACE("Delete outgoing " << kwh << " " << size_t(nucl) << " " << kwh.is_minimal());
    	this->get_raw_value_reference(kwh).DeleteOutgoing(nucl, kwh.is_minimal());
    }

    void DeleteIncoming(const KeyWithHash &kwh, char nucl) {
        TRACE("Delete incoming " << kwh << " " << size_t(nucl) << " " << kwh.is_minimal());
    	this->get_raw_value_reference(kwh).DeleteIncoming(nucl, kwh.is_minimal());
    }

    void IsolateVertex(const KeyWithHash &kwh) {
        TRACE("Isolate vertex " << kwh);
        this->get_raw_value_reference(kwh).IsolateVertex();
    }

    bool CheckOutgoing(const KeyWithHash &kwh, char nucl) const {
    	return this->get_value(kwh).CheckOutgoing(nucl);
    }

    KeyWithHash GetOutgoing(const KeyWithHash &kwh, char nucl) const {
        return kwh << nucl;
    }

    bool CheckIncoming(const KeyWithHash &kwh, char nucl) const {
    	return this->get_value(kwh).CheckIncoming(nucl);
    }

    KeyWithHash GetIncoming(const KeyWithHash &kwh, char nucl) const {
        return kwh >> nucl;
    }

    bool IsDeadEnd(const KeyWithHash &kwh) const {
    	return this->get_value(kwh).IsDeadEnd();
    }

    bool IsDeadStart(const KeyWithHash &kwh) const {
    	return this->get_value(kwh).IsDeadStart();
    }

    bool CheckUniqueOutgoing(const KeyWithHash &kwh) const {
    	return this->get_value(kwh).CheckUniqueOutgoing();
    }

    KeyWithHash GetUniqueOutgoing(const KeyWithHash &kwh) const {
    	return GetOutgoing(kwh, this->get_value(kwh).GetUniqueOutgoing());
    }

    bool CheckUniqueIncoming(const KeyWithHash &kwh) const {
    	return this->get_value(kwh).CheckUniqueIncoming();
    }

    KeyWithHash GetUniqueIncoming(const KeyWithHash &kwh) const {
    	return GetIncoming(kwh, this->get_value(kwh).GetUniqueIncoming());
    }

    size_t OutgoingEdgeCount(const KeyWithHash &kwh) const {
    	return this->get_value(kwh).OutgoingEdgeCount();
    }

    size_t IncomingEdgeCount(const KeyWithHash &kwh) const {
    	return this->get_value(kwh).IncomingEdgeCount();
    }

    ~DeBruijnExtensionIndex() {
    }

private:
   DECL_LOGGER("ExtentionIndex");
};

template<class Builder>
class DeBruijnExtensionIndexBuilder : public Builder {
    typedef Builder base;
public:
    typedef typename Builder::IndexT IndexT;

    template<class ReadStream>
    size_t FillExtensionsFromStream(ReadStream &stream, IndexT &index) const {
        unsigned k = index.k();
        size_t rl = 0;

        while (!stream.eof()) {
            typename ReadStream::read_type r;
            stream >> r;
            rl = std::max(rl, r.size());

            const Sequence &seq = r.sequence();
            if (seq.size() < k + 1)
                continue;

            typename IndexT::KeyWithHash kwh = index.ConstructKWH(seq.start<runtime_k::RtSeq>(k));
            for (size_t j = k; j < seq.size(); ++j) {
                char nnucl = seq[j], pnucl = kwh[0];
                index.AddOutgoing(kwh, nnucl);
                kwh <<= nnucl;
                index.AddIncoming(kwh, pnucl);
            }
        }

        return rl;
    }

    void FillExtensionsFromIndex(const std::string &KPlusOneMersFilename,
                                 IndexT &index) const {
        unsigned KPlusOne = index.k() + 1;

        typename IndexT::kmer_iterator it(
                KPlusOneMersFilename, runtime_k::RtSeq::GetDataSize(KPlusOne));
        for (; it.good(); ++it) {
            runtime_k::RtSeq kpomer(KPlusOne, *it);

            char pnucl = kpomer[0], nnucl = kpomer[KPlusOne - 1];
            TRACE("processing k+1-mer " << kpomer);
            index.AddOutgoing(index.ConstructKWH(runtime_k::RtSeq(KPlusOne - 1, kpomer)),
                              nnucl);
            // FIXME: This is extremely ugly. Needs to add start / end methods to extract first / last N symbols...
            index.AddIncoming(index.ConstructKWH(runtime_k::RtSeq(KPlusOne - 1, kpomer << 0)),
                              pnucl);
        }
    }

public:
    template<class Streams>
    ReadStatistics BuildExtensionIndexFromStream(
            IndexT &index, Streams &streams, io::SingleStream* contigs_stream = 0,
            size_t read_buffer_size = 0) const {
        unsigned nthreads = (unsigned) streams.size();

        // First, build a k+1-mer index
        DeBruijnReadKMerSplitter<typename Streams::ReadT, StoringTypeFilter<typename IndexT::storing_type>> splitter(
                index.workdir(), index.k() + 1, 0xDEADBEEF, streams,
                contigs_stream, read_buffer_size);
        KMerDiskCounter<runtime_k::RtSeq> counter(index.workdir(), splitter);
        counter.CountAll(nthreads, nthreads, /* merge */false);

        // Now, count unique k-mers from k+1-mers
        DeBruijnKMerKMerSplitter<StoringTypeFilter<typename IndexT::storing_type> > splitter2(index.workdir(), index.k(),
                                           index.k() + 1, IndexT::storing_type::IsInvertable(), read_buffer_size);
        for (unsigned i = 0; i < nthreads; ++i)
            splitter2.AddKMers(counter.GetMergedKMersFname(i));
        KMerDiskCounter<runtime_k::RtSeq> counter2(index.workdir(), splitter2);

        index.BuildIndex(counter2, 16, nthreads);

        // Build the kmer extensions
        INFO("Building k-mer extensions from k+1-mers");
#       pragma omp parallel for num_threads(nthreads)
        for (unsigned i = 0; i < nthreads; ++i)
            FillExtensionsFromIndex(counter.GetMergedKMersFname(i), index);
        INFO("Building k-mer extensions from k+1-mers finished.");

        return splitter.stats();
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
