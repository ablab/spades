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
template<class Seq = runtime_k::RtSeq,
    class traits = kmer_index_traits<Seq> >
class DeBruijnExtensionIndex : public DeBruijnKMerIndex<uint8_t, Seq, traits> {
private:
  typedef DeBruijnKMerIndex<uint8_t, Seq, traits> base;

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
  typedef Seq                      KMer;
  typedef KMerIndex<KMer, traits>  KMerIndexT;

  DeBruijnExtensionIndex(unsigned K, const std::string &workdir)
      : base(K, workdir) {}

  void AddOutgoing(size_t idx, char nnucl) {
	  unsigned nmask = (1 << nnucl);
	  if (!(this->operator [](idx) & nmask)) {
#         pragma omp atomic
		  this->operator [](idx) |= nmask;
	  }
  }

  void AddIncoming(size_t idx, char pnucl) {
	  unsigned pmask = (1 << (pnucl + 4));

	  if (!(this->operator [](idx) & pmask)) {
#         pragma omp atomic
		  this->operator [](idx) |= pmask;
	  }
  }

  void DeleteOutgoing(size_t idx, char nnucl) {
	  unsigned nmask = (1 << nnucl);
	  if (this->operator [](idx) & nmask) {
#         pragma omp atomic
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
  	Seq kmer;
  	size_t idx;

  	KmerWithHash(
  			Seq _kmer,
  			const DeBruijnExtensionIndex<Seq, traits> &index) :
  			kmer(_kmer), idx(index.seq_idx(kmer)) {
  	}
  };

  KmerWithHash CreateKmerWithHash(Seq kmer) const {
  		return KmerWithHash(kmer, *this);
  }
};

template <class Seq>
class DeBruijnExtensionIndexBuilder : public DeBruijnKMerIndexBuilder<Seq> {
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
      public DeBruijnKMerIndexBuilder<runtime_k::RtSeq> {
  typedef DeBruijnKMerIndexBuilder<runtime_k::RtSeq> base;

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
//      size_t idx = index.seq_idx(kmer);

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
  template <class Read>
  size_t BuildIndexFromStream(DeBruijnExtensionIndex<runtime_k::RtSeq> &index,
                              io::ReadStreamVector<io::IReader<Read> > &streams,
                              SingleReadStream* contigs_stream = 0) const {
    unsigned nthreads = streams.size();

    base::BuildIndexFromStream(index, streams, contigs_stream);

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

 protected:
  DECL_LOGGER("Extension Index Building");
};

}
