#pragma once
//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "kmer_index/kmer_mph/kmer_index_builder.hpp"
#include "kmer_index/kmer_mph/kmer_buckets.hpp"

#include "pipeline/partask_mpi.hpp"

namespace kmers {
template<class Index>
class KMerIndexBuilderMPI {
    typedef typename Index::KMerSeq Seq;
    typedef typename Index::kmer_index_traits kmer_index_traits;

 public:
    KMerIndexBuilderMPI() = default;
    size_t BuildIndexMPI(Index &out, KMerDiskStorage<Seq>& kmerstorage, bool save_final = false);

 private:
    DECL_LOGGER("K-mer Index Building MPI");

    class BuildKMerIndexTask {
        BuildKMerIndexTask() = default;
     public:
        BuildKMerIndexTask(KMerDiskStorage<Seq> kmerstorage, unsigned k) : storage_{std::move(kmerstorage)}, k_{k} {};
        BuildKMerIndexTask(std::istream &is) { deserialize(is); }

        std::ostream &serialize(std::ostream &os) const {
            storage_.BinWrite(os);
            return os;
        }

        std::istream &deserialize(std::istream &is) {
            storage_.BinRead(is);
            return is;
        }

        auto make_splitter(size_t, Index &) {
            return partask::make_seq_generator(storage_.num_buckets());
        }

        void process(std::istream &is, std::ostream &os, Index &) {
            std::vector<size_t> residuals;
            while (is.get() && is) {
                size_t i;
                io::binary::BinRead(is, i);
                DEBUG("BuildKMerIndexTask: process, i = " << i);
                residuals.push_back(i);
            }

            size_t num_buckets = storage_.num_buckets();

            std::vector<typename Index::KMerDataIndex *> indices(num_buckets);
            std::vector<size_t> sizes(num_buckets);
            DEBUG("NumBuckets: " << num_buckets);
#pragma omp parallel for
            for (size_t ii = 0; ii < residuals.size(); ++ii) {
                size_t i = residuals[ii];
                const auto &bucket_range = storage_.bucket(i);
                sizes[i] = storage_.bucket_size(i);
                DEBUG("Bucket size: " << sizes[i]);
                indices[i] = new typename Index::KMerDataIndex(sizes[i], Index::KMerDataIndex::ConflictPolicy::Ignore, /* gamma */ 4.0);
                indices[i]->build(boomphf::range(storage_.bucket_begin(i), storage_.bucket_end(i)));
                DEBUG("Index created");
            }

            for (size_t i = 0; i < indices.size(); ++i) {
                if (!indices[i]) continue;
                DEBUG("Sending index " << i);
                os.put('\1');
                io::binary::BinWrite(os, i, sizes[i]);
                if (partask::master()) {
                    io::binary::BinWrite(os, indices[i]);  // Pass the pointer
                } else {
                    indices[i]->save(os);
                    delete indices[i];
                }
            }
            os.put('\0');
        }

        void merge(const std::vector<std::istream *> &piss, Index &index) {
            INFO("Index merge started");
            index.clear();

            auto segment_policy = storage_.segment_policy();
            size_t segments = segment_policy.num_segments();
            index.segment_starts_.resize(segments + 1, 0);
            index.index_.resize(storage_.num_buckets(),
                                typename Index::KMerDataIndex(0, Index::KMerDataIndex::ConflictPolicy::Ignore, /* gamma */ 4.0));
            index.segment_policy_ = segment_policy;
            index.num_segments_ = segments;

            DEBUG("Index initialized with empty subindices");

            for (size_t node = 0; node < piss.size(); ++node) {
                DEBUG("Collecting stream " << node);
                std::istream &is = *piss[node];
                while (static_cast<char>(is.get())) {
                    size_t i, size;
                    io::binary::BinRead(is, i, size);
                    DEBUG("Load index " << i << " " << size << " from stream " << node);
                    if (node == 0) {
                        typename KMerIndex<kmer_index_traits>::KMerDataIndex *p;
                        io::binary::BinRead(is, p);
                        std::swap(index.index_[i], *p);
                        delete p;
                    } else {
                        index.index_[i].load(is);
                    }
                    DEBUG("Index loaded " << i << " " << size);
                    index.segment_starts_[i + 1] = size;
                }
            }

            // Finally, record the sizes of buckets.
            for (size_t i = 1; i <= segments; ++i) {
                index.segment_starts_[i] += index.segment_starts_[i - 1];
            }

            INFO("Index merge done");
        }

        void sync(Index &index) {
            INFO("Index broadcasting started");
            partask::broadcast(index,
                               [](std::ostream &os, const auto &index) { return index.serialize(os); },
                               [](std::istream &is, auto &index) { return index.deserialize(is); });
            INFO("Index broadcasting done");
        }

     private:
        KMerDiskStorage<Seq> storage_;
        unsigned k_;
    };
};

template<class Index>
size_t KMerIndexBuilderMPI<Index>::BuildIndexMPI(Index &index, KMerDiskStorage<Seq>& storage, bool save_final) {
    index.clear();

    INFO("Building kmer index ");

    size_t buckets = storage.num_buckets();
    unsigned k = storage.k();
    DEBUG("K: " << k);
    partask::TaskRegistry treg;
    auto build_index = treg.add<BuildKMerIndexTask>(std::ref(index));

    treg.listen();
    if (partask::master()) {
        build_index(storage, k);
    }
    treg.stop_listening();

    size_t kmers = index.segment_starts_[buckets];
    double bits_per_kmer = 8.0 * (double) index.mem_size() / (double) kmers;
    INFO("Index built. Total " << kmers << " k-mers, " << index.mem_size() << " bytes occupied (" << bits_per_kmer
                               << " bits per kmer), " << buckets << " buckets");

    if (partask::master() && save_final) {
        storage.merge();
    }

    return kmers;
}
}
