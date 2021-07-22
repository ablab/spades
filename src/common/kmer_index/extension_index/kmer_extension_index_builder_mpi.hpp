//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "kmer_extension_index.hpp"
#include "kmer_extension_index_builder.hpp"

#include "kmer_index/kmer_mph/kmer_index_builder.hpp"
#include "kmer_index/kmer_mph/kmer_index_builder_mpi.hpp"
#include "kmer_index/kmer_mph/kmer_splitters.hpp"
#include "kmer_index/kmer_counting.hpp"
#include "kmer_index/ph_map/perfect_hash_map_builder.hpp"
#include "io/reads/multifile_reader.hpp"

#include "pipeline/partask_mpi.hpp"

namespace kmers {
class DeBruijnExtensionIndexBuilderMPI : public DeBruijnExtensionIndexBuilder {
 public:
    template<class Index, class Streams>
    void BuildExtensionIndexFromStream(fs::TmpDir workdir, Index &index,
                                       Streams &streams,
                                       size_t read_buffer_size = 0) const {
        unsigned nthreads = omp_get_max_threads();
        using KmerFilter = StoringTypeFilter<typename Index::storing_type>;

        // First, build a k+1-mer index
        DeBruijnReadKMerSplitter<typename Streams::ReadT, KmerFilter>
            splitter(workdir, index.k() + 1, 0xDEADBEEF, streams, read_buffer_size);
        kmers::KMerDiskCounter<RtSeq> counter(workdir, splitter);
        auto storage = counter.CountAll(nthreads, /* merge */ false);

        BuildExtensionIndexFromKPOMers(workdir, index, storage,
                                       nthreads, read_buffer_size);
    }

    template<class Index, class KMerStorage>
    void BuildExtensionIndexFromKPOMersMPI(fs::TmpDir workdir,
                                           Index &index, KMerStorage &kmerstorage,
                                           unsigned nbuckets, size_t read_buffer_size = 0) const;

 private:
    DECL_LOGGER("DeBruijnExtensionIndexBuilderMPI");
};

template<class Index, class Seq>
class FillIndexTask {
 public:
    FillIndexTask() = default;

    FillIndexTask(kmers::KMerDiskStorage<Seq> &storage) : storage_{storage} {};
    FillIndexTask(std::istream &is) {
        storage_.BinRead(is);
    }

    std::ostream &serialize(std::ostream &os) const {
        io::binary::BinWrite(os, storage_);
        return os;
    }

    auto make_splitter(size_t, Index &) {
        return partask::make_seq_generator(storage_.num_buckets());
    }

    void process(std::istream &is, std::ostream & /*os*/, Index &index) {
        auto file_ids = partask::get_seq(is);
#       pragma omp parallel for
        for (size_t i = 0; i < file_ids.size(); ++i) {
            size_t idx = file_ids[i];
            builder_.FillExtensionsFromIndex(storage_.bucket_begin(idx), storage_.bucket_end(idx), index);
        }

        // Send nothing
    }

    void sync(Index &index) {
        INFO("FillIndexTask::sync started");
        partask::allreduce(index.raw_data(), index.raw_size(), MPI_BOR);
        INFO("FillIndexTask::sync finished");
    }

 private:
    kmers::KMerDiskStorage<Seq> storage_;
    DeBruijnExtensionIndexBuilder builder_;
};

template<class Index>
class SplitKPOMersTask {
    typedef typename Index::traits_t::SeqType Seq;
    SplitKPOMersTask() = default;
 public:
    SplitKPOMersTask(kmers::KMerDiskStorage<Seq> &storage,
                     unsigned k,
                     unsigned nthreads,
                     size_t read_buffer_size,
                     const std::string &dir)
        : storage_{storage}, k_{k}, nthreads_{nthreads}, read_buffer_size_{read_buffer_size}, dir_{dir} {};
    SplitKPOMersTask(std::istream &is) {
        io::binary::BinRead(is, storage_, k_, nthreads_, read_buffer_size_, dir_);

    }

    std::ostream &serialize(std::ostream &os) const {
        io::binary::BinWrite(os, storage_, k_, nthreads_, read_buffer_size_, dir_);
        return os;
    }

    auto make_splitter(size_t) {
        return partask::make_seq_generator(storage_.num_buckets());
    }

    void process(std::istream &is, std::ostream &os) {
        auto workdir = fs::tmp::acquire_temp_dir(dir_);
        workdir->release();
        DeBruijnKMerKMerSplitter<StoringTypeFilter<typename Index::storing_type>,
                                 typename kmers::KMerDiskStorage<Seq>::kmer_iterator> splitter(
            workdir, k_, k_ + 1, Index::storing_type::IsInvertable(), read_buffer_size_);

        policy_ = splitter.bucket_policy();
        auto file_ids = partask::get_seq(is);
        for (size_t i : file_ids) {
            splitter.AddKMers(storage_.bucket(i));
        }

        kmers::KMerDiskCounter<RtSeq> counter2(workdir, splitter);
        auto storage2 = counter2.CountAll(storage_.num_buckets(), nthreads_,/* merge */ false);
        storage2.BinWrite(os);
        storage2.release_all();
    }

    auto merge(const std::vector<std::istream *> &piss) {
        auto workdir = fs::tmp::acquire_temp_dir(dir_);
        workdir->release();
        std::vector<kmers::KMerDiskStorage<Seq>> storages;
        for (size_t i = 0; i < piss.size(); ++i) {
            auto &is = *piss[i];
            kmers::KMerDiskStorage<Seq> kmerstorage(workdir, k_, policy_);
            kmerstorage.BinRead(is);
            storages.push_back(std::move(kmerstorage));
        }

        return storages;
    }

 private:
    kmers::KMerDiskStorage<Seq> storage_;
    typename kmers::KMerDiskStorage<Seq>::KMerSegmentPolicy policy_;
    unsigned k_;
    unsigned nthreads_;
    size_t read_buffer_size_;
    std::string dir_;
};

template<class Index, class KMerStorage>
inline void DeBruijnExtensionIndexBuilderMPI::BuildExtensionIndexFromKPOMersMPI(fs::TmpDir workdir,
                                                                                Index &index,
                                                                                KMerStorage &kmerstorage,
                                                                                unsigned nthreads,
                                                                                size_t read_buffer_size) const {
    typedef typename Index::traits_t::SeqType Seq;
    VERIFY(kmerstorage.k() == index.k() + 1);

    KMerStorage kmerfiles2(workdir, index.k(), kmerstorage.segment_policy());

    INFO("DeBruijnExtensionIndexBuilder started nthreads = " << nthreads);
    partask::TaskRegistry treg;
    auto merge_kmer_files = treg.add<kmers::MergeKMerFilesTask<Seq>>();
    auto split_kpo_mers = treg.add<SplitKPOMersTask<Index>>();
    treg.listen();
    DEBUG("Listening started");

    if (partask::master()) {
        std::vector<std::string> outputfiles;
        DEBUG("Split_kpo_mers started");
        auto unmerged_kmerfiles2 = split_kpo_mers(kmerstorage, index.k(), nthreads, read_buffer_size, workdir->dir());
        DEBUG("Split_kpo_mers finished");

        for (unsigned i = 0; i < kmerfiles2.num_buckets(); ++i) {
            outputfiles.push_back(kmerfiles2.create(i)->file());
        }

        merge_kmer_files(std::move(unmerged_kmerfiles2), outputfiles, index.k(), workdir->dir());
        DEBUG("Merge_kmer_files finished");
    }
    treg.stop_listening();

    partask::broadcast(kmerfiles2);
    INFO("Total kmers=" << kmerfiles2.total_kmers());
    BuildIndexMPI(index, kmerfiles2, /* save_final */ true);

    auto fill_index = treg.add<FillIndexTask<Index, Seq>>(std::ref(index));
    treg.listen();
    DEBUG("Listening started");

    if (partask::master()) {
        fill_index(kmerstorage);
    }
    treg.stop_listening();

    INFO("Building k-mer extensions from k+1-mers finished.");
}
} //utils
