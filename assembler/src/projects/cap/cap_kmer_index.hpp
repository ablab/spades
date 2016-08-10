//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "compare_standard.hpp"
#include "longseq.hpp"
#include "polynomial_hash.hpp"
#include "common/adt/kmer_map.hpp"
#include "utils/indices/edge_position_index.hpp"

#include "io/reads/sequence_reader.hpp"
#include "utils/mph_index/base_hash.hpp"

template<>
struct kmer_index_traits<cap::LSeq> {
    typedef cap::LSeq SeqType;
    typedef std::vector<cap::LSeq> RawKMerStorage;
    typedef std::vector<cap::LSeq> FinalKMerStorage;

    typedef RawKMerStorage::iterator             raw_data_iterator;
    typedef RawKMerStorage::const_iterator       raw_data_const_iterator;
    typedef RawKMerStorage::iterator::value_type KMerRawData;
    typedef RawKMerStorage::iterator::reference  KMerRawReference;

    struct raw_equal_to {
        inline bool operator()(const SeqType &lhs, const KMerRawReference rhs) {
            // Using fast_equal_to, which relies only on hash:
            // 1. True comparison leads to poor performance on large k
            // 2. Hashes are to be different (in other case MPH is impossible)
            return SeqType::fast_equal_to()(lhs, rhs);
        }
    };

    struct hash_function {
        inline uint64_t operator()(const SeqType &k) const {
            return k.GetHash().get<2>();
        }
        inline uint64_t operator()(const KMerRawReference k) const {
            return k.GetHash().get<2>();
        }
    };

    struct seeded_hash_function {
        static cxxmph::h128 hash128(const KMerRawReference k, uint32_t seed) {
            SeqType::HashType hash = k.GetHash();
            //      uint64_t salt = hash.get<2>();
            cxxmph::h128 h;
            MurmurHash3_x64_128(reinterpret_cast<const void*>(&hash), sizeof(hash),
                    seed, &h);
            //      h.set64(hash.get<0>(), 0);
            //      h.set64(hash.get<1>(), 1);

            //      INFO("SEE MAN:: hash=" << hash.get<0>() << ", seed=" << seed << ", result = " << h.get64(0) << " " << h.get64(1) << " " << k.str());
            return h;
        }

        static cxxmph::h128 hash128(const SeqType &k, uint32_t seed) {
            SeqType::HashType hash = k.GetHash();
            //      uint64_t salt = hash.get<2>();
            cxxmph::h128 h;
            MurmurHash3_x64_128(reinterpret_cast<const void*>(&hash), sizeof(hash),
                    seed, &h);
            //      h.set64(hash.get<0>(), 0);
            //      h.set64(hash.get<1>(), 1);

            //      INFO("SEE MAN:: hash=" << hash.get<0>() << ", seed=" << seed << ", result = " << h.get64(0) << " " << h.get64(1) << " " << k.str());

            return h;
        }
    };

    struct raw_create {
        inline SeqType operator()(unsigned /*K*/, const KMerRawReference kmer) {
            return SeqType(kmer);
        }
    };

    template <class Reader>
        static RawKMerStorage *raw_deserialize(Reader &/*reader*/, const std::string &/*FileName*/) {
            VERIFY(false);
            return NULL;
        }

};

namespace cap {

    template <class ReadType>
        class CapKMerCounter : public ::KMerCounter<LSeq> {
            typedef KMerCounter<LSeq> __super;
            typedef typename __super::RawKMerStorage RawKMerStorage;

            unsigned k_;
            io::ReadStreamList<ReadType> streams_;
            std::unordered_set<LSeq, LSeq::hash, LSeq::equal_to> storage_;
            RawKMerStorage *bucket;

            bool has_counted_;

            public:
            CapKMerCounter(const unsigned k, io::ReadStreamList<ReadType> streams)
                : k_(k),
                streams_(streams),
                storage_(),
                bucket(NULL),
                has_counted_(false) {
                }

            CapKMerCounter(const unsigned k) : CapKMerCounter(k, NULL) {
            }

            virtual ~CapKMerCounter() {
                ReleaseBucket(0);
            }

            virtual size_t KMerSize() const {
                return LSeq::GetDataSize(k_) * sizeof(typename LSeq::DataType);
            }

            virtual size_t Count(unsigned /* num_buckets */, unsigned /* num_threads */) {
                if (!has_counted_) {
                    Init();
                    has_counted_ = true;
                }
                INFO("K-mer counting done. There are " << storage_.size() << " kmers in total. ");
                return storage_.size();
            }

            virtual size_t CountAll(unsigned /* num_buckets */, unsigned /* num_threads */, bool /* merge  */= true) {
                if (!has_counted_) {
                    Init();
                    has_counted_ = true;
                }
                INFO("K-mer counting done. There are " << storage_.size() << " kmers in total. ");
                return storage_.size();
            }

            virtual void MergeBuckets(unsigned /* num_buckets */) {
                VERIFY(bucket == NULL);
            }

            virtual void OpenBucket(size_t /* idx */, bool /* unlink  */= true) {
                VERIFY(bucket == NULL);

                if (!has_counted_) {
                    Init();
                    has_counted_ = true;
                }
                TRACE("BUCKET OPEN");
                bucket = new RawKMerStorage();
                bucket->reserve(storage_.size());
                for (auto it = storage_.begin(); it != storage_.end(); ++it) {
                    bucket->push_back(*it);
                }
            }

            virtual void ReleaseBucket(size_t /* idx */) {
                TRACE("RELEASE BUCKET");
                delete bucket;
                bucket = NULL;
            }

            virtual RawKMerStorage* GetFinalKMers() {
                OpenBucket(0);
                VERIFY(bucket != NULL);

                RawKMerStorage *ret = bucket;
                bucket = NULL;

                return ret;
            }

            virtual typename __super::iterator bucket_begin(size_t /* idx */) {
                return bucket->begin();
            }
            virtual typename __super::iterator bucket_end(size_t /* idx */) {
                return bucket->end();
            }

            protected:
            virtual void Init() {
                VERIFY(streams_.size() > 0);
                for (size_t i = 0; i < streams_.size(); ++i) {
                    while (!streams_[i].eof()) {
                        ReadType r;
                        streams_[i] >> r;
                        const Sequence &seq = r.sequence();
                        if (seq.size() == 0) {
                            continue;
                        }
                        if (seq.size() < k_) {
                            INFO("WARNING: too small sequence!!");
                            continue;
                        }

                        LSeq kmer(k_, seq);
                        do {
                            storage_.insert(kmer);
                            kmer.Shift();
                        } while (kmer.IsValid());

                    }
                }
                streams_.clear();
            }

            void SetStreams(io::ReadStreamList<ReadType>& streams) {
                streams_ = streams;
            }

        };

    template <class Graph>
        class CapKMerGraphCounter : public CapKMerCounter<io::SingleRead> {

            public:
            CapKMerGraphCounter(const unsigned k, const Graph &g)
                : CapKMerCounter<io::SingleRead>(k),
                g_(g) {

                }

            protected:
            virtual void Init() {
                io::ReadStreamList<io::SingleRead> stream_vector;
                //fixme create reasonable reader from the graph
                for (auto it = g_.ConstEdgeBegin(); !it.IsEnd(); ++it) {
                    stream_vector.push_back(make_shared<io::SequenceReadStream<io::SingleRead>>(g_.EdgeNucls(*it)));
                }

                CapKMerCounter<io::SingleRead>::SetStreams(stream_vector);
                CapKMerCounter<io::SingleRead>::Init();
            }

            private:
            const Graph &g_;
        };

}

namespace debruijn_graph {

    template<class Index>
        class DeBruijnStreamKMerIndexBuilder<cap::LSeq, Index> {
            public:
                typedef Index IndexT;

                template <class Streams>
                    size_t BuildIndexFromStream(Index &index,
                            Streams &streams,
                            io::SingleStream* /* contigs_stream  */= 0) const {
                        /*
                           std::vector<io::IReader<io::SingleRead> *> stream_vec(streams.size());
                           for (size_t i = 0; i < streams.size(); ++i) {
                           stream_vec.push_back(&streams[i]);
                           }
                           auto streams_ptr = std::make_shared<Streams>(
                           new io::ReadStreamVector<io::IReader<io::SingleRead>>(stream_vec, false));
                           */
                        cap::CapKMerCounter<typename Streams::ReadT> counter(
                                index.k(), streams);

                        index.BuildIndex(counter, 1, 1);
                        return 0;
                    }

        };

    template <class Index>
        class DeBruijnGraphKMerIndexBuilder<Index, cap::LSeq> {
                  public:
                      typedef Index IndexT;

                      template<class Graph>
                          void BuildIndexFromGraph(IndexT &index, const Graph &g) const {
                              cap::CapKMerGraphCounter<Graph> counter(index.k(), g);
                              index.BuildIndex(counter, 16, 1);
                          }
              };

}

namespace runtime_k {

    //todo review this class
    template <typename Value>
        class KmerMap<Value, cap::LSeq> {
            public:
                typedef typename cap::LSeq key_type;
                typedef typename std::pair<const key_type, Value> value_type;
                typedef KmerMap<Value, cap::LSeq> map_type;

            private:
                // Note using equal_to which maintains special 'transitions' inside LongSeqs
                typedef std::unordered_map<key_type, Value, typename key_type::hash, typename key_type::equal_to> int_map_type;
                typedef typename std::pair<const key_type, const Value> const_value_type;
                int_map_type *data_;

                class InnerIterator {
                    friend class KmerMap<Value, cap::LSeq>;
                    typedef typename int_map_type::iterator base;
                    base base_;

                    public:

                    InnerIterator(const base &iter): base_(iter) {
                    }

                    InnerIterator &operator++() {
                        ++base_;
                        return *this;
                    }
                    InnerIterator operator++(int) {
                        InnerIterator stored = *this;
                        ++base_;
                        return stored;
                    }

                    value_type operator*() {
                        return base_->operator *();
                    }

                    const key_type first() {
                        return base_->first();
                    }

                    Value& second() {
                        return base_->second();
                    }

                    bool operator==(const InnerIterator& iter) const {
                        return base_ == iter.base_;
                    }
                    bool operator!=(const InnerIterator& iter) const {
                        return !operator==(iter);
                    }

                };

                class InnerConstIterator {
                    friend class KmerMap<Value, cap::LSeq>;
                    typedef typename int_map_type::const_iterator base;
                    base base_;

                    public:

                    InnerConstIterator(const base &iter): base_(iter) {
                    }

                    InnerConstIterator &operator++() {
                        ++base_;
                        return *this;
                    }
                    InnerConstIterator operator++(int) {
                        InnerConstIterator stored = *this;
                        ++base_;
                        return stored;
                    }

                    const value_type operator*() const {
                        return base_->operator *();
                    }

                    key_type first() const {
                        return base_->first;
                    }

                    const Value& second() const {
                        return base_->second;
                    }

                    bool operator==(const InnerConstIterator& iter) const {
                        return base_ == iter.base_;
                    }
                    bool operator!=(const InnerConstIterator& iter) const {
                        return !operator==(iter);
                    }

                };

            public:
                typedef InnerIterator iterator;
                typedef InnerConstIterator const_iterator;

                KmerMap(size_t /* k */) {
                    data_ = new int_map_type();
                }

                KmerMap(int_map_type *map) : data_(map) {
                }

                KmerMap(const map_type& map) {
                    data_ = new int_map_type(*map.data);
                }

                map_type& operator=(const map_type& map) {
                    if (map.data_ != data_) {
                        delete data_;
                        data_ = new int_map_type(*map.data_);
                    }

                    return *this;
                }

                ~KmerMap() {
                    delete data_;
                }

                bool empty() const {
                    return data_->empty();
                }

                size_t size() const {
                    return data_->size();
                }

                size_t max_size() const {
                    return data_->max_size();
                }

                const_iterator begin() const {
                    return InnerConstIterator(data_->begin());
                }

                iterator begin() {
                    return InnerIterator(data_->begin());
                }

                const_iterator end() const {
                    return InnerConstIterator(data_->end());
                }

                iterator end() {
                    return InnerIterator(data_->end());
                }

                Value& operator[](const key_type& kmer_seq) {
                    return data_->operator [](kmer_seq);
                }

                const_iterator find(const key_type& kmer_seq) const {
                    return InnerConstIterator(data_->find(kmer_seq));
                }

                iterator find(const key_type& kmer_seq) {
                    return InnerIterator(data_->find(kmer_seq));
                }

                size_t count(const key_type& kmer_seq) const {
                    return data_->count(kmer_seq);
                }

                pair<iterator, bool> insert(const value_type& val) {
                    return data_->insert(val);
                }

                size_t erase(const key_type& kmer_seq) {
                    return data_->erase(kmer_seq);
                }

                //    iterator erase(const const_iterator& iter) {
                //        return iterator(data_->erase(iter.get_data()));
                //    }

                iterator erase(const iterator& iter) {
                    return data_->erase(iter.base_);
                }

                void clear() {
                    data_->clear();
                }

                /*
                   size_t bucket_count() const {
                   return data_->bucket_count();
                   }

                   size_t max_bucket_count() const {
                   return data_->max_bucket_count();
                   }

                   size_t bucket_size(size_t n) const {
                   return data_->bucket_size(n);
                   }

                   size_t bucket(const RtSeq& kmer_seq) const {
                   return data_->bucket(kmer_seq);
                   }

                   float load_factor() const {
                   return data_->load_factor();
                   }

                   float max_load_factor() const {
                   return data_->max_load_factor();
                   }

                   void max_load_factor(float z) {
                   data_->max_load_factor(z);
                   }

                   void rehash(size_t n) {
                   data_->rehash(n);
                   }

                   size_t get_k() const {
                   return data_->get_k();
                   }

                   int_map_type& get_data() {
                   return *data_;
                   }
                   */



        };

}
