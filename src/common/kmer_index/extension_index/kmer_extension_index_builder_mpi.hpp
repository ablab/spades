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

    /*
     * Build extension index from k+1-mers.
     *
     * kpostorage -- storage with k+1 mers. This storage is unmerged.
     * */
    template<class Index, class KMerStorage>
    void BuildExtensionIndexFromKPOMersMPI(fs::TmpDir workdir,
                                           Index &index, const KMerStorage &kpostorage,
                                           unsigned nbuckets, size_t read_buffer_size = 0) const;

 private:
    DECL_LOGGER("DeBruijnExtensionIndexBuilderMPI");
};

template<class Index, class Seq>
class FillIndexTask {
 public:
    FillIndexTask() = default;

    FillIndexTask(const kmers::KMerDiskStorage<Seq> &kpostorage) : kpostorage_{kpostorage} {};
    FillIndexTask(std::istream &is) {
        kpostorage_.BinRead(is);
    }

    std::ostream &serialize(std::ostream &os) const {
        io::binary::BinWrite(os, kpostorage_);
        return os;
    }

    auto make_splitter(size_t, Index &) {
        return partask::make_seq_generator(kpostorage_.num_buckets());
    }

    void process(std::istream &is, std::ostream & /*os*/, Index &index) {
        auto file_ids = partask::get_seq(is);
#       pragma omp parallel for
        for (size_t i = 0; i < file_ids.size(); ++i) {
            size_t idx = file_ids[i];
            builder_.FillExtensionsFromIndex(kpostorage_.bucket_begin(idx), kpostorage_.bucket_end(idx), index);
        }

        // Send nothing
    }

    void sync(Index &index) {
        INFO("FillIndexTask::sync started");
        partask::allreduce(index.raw_data(), index.raw_size(), MPI_BOR);
        INFO("FillIndexTask::sync finished");
    }

 private:
    kmers::KMerDiskStorage<Seq> kpostorage_;
    DeBruijnExtensionIndexBuilder builder_;
};


/*
 * Build kmer storages from k+1-mer storage
 *
 * Return the vector of kmer storages one for each node. Each kmer storage contain unique sorted kmers.
 * The segmentation policy for each returned kmer storage the same as policy for k+1-mer storage.
 * All returned storages are unmerged.
 */
template<class Index>
class SplitKPOMersTask {
    typedef typename Index::traits_t::SeqType Seq;
    SplitKPOMersTask() = default;
 public:
    SplitKPOMersTask(const kmers::KMerDiskStorage<Seq> &kpostorage,
                     unsigned k,
                     unsigned nthreads,
                     size_t read_buffer_size,
                     const std::string &dir)
        : kpostorage_{kpostorage}, k_{k}, nthreads_{nthreads}, read_buffer_size_{read_buffer_size}, dir_{dir} {};

    SplitKPOMersTask(std::istream &is) {
        io::binary::BinRead(is, kpostorage_, k_, nthreads_, read_buffer_size_, dir_);

    }

    std::ostream &serialize(std::ostream &os) const {
        io::binary::BinWrite(os, kpostorage_, k_, nthreads_, read_buffer_size_, dir_);
        return os;
    }

    auto make_splitter(size_t) {
        return partask::make_seq_generator(kpostorage_.num_buckets());
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
            splitter.AddKMers(kpostorage_.bucket(i));
        }

        kmers::KMerDiskCounter<RtSeq> counter2(workdir, splitter);
        auto storage2 = counter2.CountAll(kpostorage_.num_buckets(), nthreads_,/* merge */ false);
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
    kmers::KMerDiskStorage<Seq> kpostorage_;
    typename kmers::KMerDiskStorage<Seq>::KMerSegmentPolicy policy_;
    unsigned k_;
    unsigned nthreads_;
    size_t read_buffer_size_;
    std::string dir_;
};


/* Merge KMerDiskStorages from a few nodes
 *
 * storages -- KMerDiskStorages which were received from different nodes after SplitKPOMersTask.
 * Each storage must be valid and unmerged. It is mean that each storage contain kmers splitted on files consist with
 * segmentation_policy, in each file all kmers are unique and sorted. However kmer can be present in different storage,
 * but in one particular storage only one time.
 *
 * MergeKMerFiles merge buckets from different Storages into one storage. The policy of splitting kmers should be the
 * same for all storages and should stay the same for final storage. In MergeKMerFilesTask all kmers from 0 buckets
 * from all storages should be merge into 0 bucket in final storage, 1 buckets should be merge into 1 bucket etc.
 *
 * Instead of final storage MergeKMerFilesTask take as a input the vector of files name(ofiles) correspondent to final bucket
 * in new storage. After the work of MergeKmerFilesTask files from ofiles should contain merged kmers. Kmers
 * in each file should be sorted and unique. 0 file should contain all kmers from all 0 buckets from all storages,
 * the 1 file should contain all kmers from 1 buckets from all storages and etc.
 */
template <typename Seq>
class MergeKMerFilesTask {
    MergeKMerFilesTask() = default;
 public:
    MergeKMerFilesTask(std::vector<kmers::KMerDiskStorage<Seq>> storages, std::vector<std::string>& ofiles) : storages_{std::move(storages)}, ofiles_{ofiles} {};
    MergeKMerFilesTask(std::istream &is) {
        io::binary::BinRead(is, storages_, ofiles_);
    }

    void serialize(std::ostream &os) const {
        io::binary::BinWrite(os, storages_, ofiles_);
    }

    auto make_splitter(size_t) {
        return partask::make_seq_generator(storages_[0].num_buckets());
    }

    void process(std::istream &is, std::ostream &) {
        std::vector<size_t> residuals;
        while (is.get() && is) {
            size_t i;
            io::binary::BinRead(is, i);
            DEBUG("MergeKMerFilesTask: process, i = " << i);
            residuals.push_back(i);
        }

        /*size_t num_open_files = omp_get_max_threads() * (2 * partask::world_size());
        INFO("Setting open file limit to " << num_open_files);
        utils::limit_file(num_open_files);*/

        int totalsum = 0;
#pragma omp parallel for
        for (size_t idx = 0; idx < residuals.size(); ++idx) {
            size_t i = residuals[idx];  // TODO rename var i -> residual

            MMappedRecordArrayWriter<typename Seq::DataType> os(ofiles_[i], Seq::GetDataSize(storages_[0].k()));
            auto elcnt = Seq::GetDataSize(storages_[0].k());
            std::vector<MMappedRecordArrayReader<typename Seq::DataType>> ins;
            std::vector<int> strs(storages_.size(), 0);
            std::vector<std::pair<int, int>> oids;

            bool notEmpty = true;
            size_t prevId = -1ULL;
            unsigned sumsize = 0;
            for (size_t sid = 0; sid < storages_.size(); ++sid) {
                ins.push_back(MMappedRecordArrayReader<typename Seq::DataType>(*storages_[sid].bucket_file(i), Seq::GetDataSize(storages_[0].k()), /* unlink */ false));
                sumsize += ins.back().size();
            }

            int total = 0;
            while (notEmpty) {
                size_t bstpos = -1ULL;
                for (size_t sid = 0; sid < storages_.size(); ++sid) {
                    if (ins[sid].size() == strs[sid]) {
                        continue;
                    }
                    if (bstpos == -1ULL || adt::array_less<typename Seq::DataType>()(*(ins[sid].begin() + strs[sid]),
                                                                                     *(ins[bstpos].begin() + strs[bstpos]))) {
                        bstpos = sid;
                    }
                }
                if (bstpos != -1ULL) {
                    if (prevId == -1ULL || adt::array_less<typename Seq::DataType>()(*(ins[prevId].begin() + (strs[prevId] - 1)), *(ins[bstpos].begin() + (strs[bstpos])))) {
                        oids.push_back(std::make_pair(bstpos, strs[bstpos]));
                        total += 1;
                    }
                    prevId = bstpos;
                    strs[bstpos] += 1;
                } else {
                    notEmpty = false;
                }
            }
            os.resize(total);
            for (auto oid : oids) {
                os.write(ins[oid.first].data() + oid.second * elcnt, 1);
            }

            totalsum += total;
        }
        DEBUG("Total kmers writen= " << totalsum);

    }

    void merge(const std::vector<std::istream *> &) {}

 private:
    std::vector<kmers::KMerDiskStorage<Seq>> storages_;
    std::vector<std::string> ofiles_;
};

/*
 * Build extension index from k+1-mers.
 *
 * kpostorage -- storage with k+1 mers. This storage is unmerged.
 */
template<class Index, class KMerStorage>
inline void DeBruijnExtensionIndexBuilderMPI::BuildExtensionIndexFromKPOMersMPI(fs::TmpDir workdir,
                                                                                Index &index,
                                                                                const KMerStorage &kpostorage,
                                                                                unsigned nthreads,
                                                                                size_t read_buffer_size) const {
    typedef typename Index::traits_t::SeqType Seq;
    VERIFY(kpostorage.k() == index.k() + 1);

    KMerStorage kmerfiles2(workdir, index.k(), kpostorage.segment_policy());

    INFO("DeBruijnExtensionIndexBuilder started nthreads = " << nthreads);
    partask::TaskRegistry treg;
    auto merge_kmer_files = treg.add<MergeKMerFilesTask<Seq>>();
    auto split_kpo_mers = treg.add<SplitKPOMersTask<Index>>();
    treg.listen();
    DEBUG("Listening started");

    if (partask::master()) {
        std::vector<std::string> outputfiles;
        INFO("Split_kpo_mers started");
        auto unmerged_kmerfiles2 = split_kpo_mers(kpostorage, index.k(), nthreads, read_buffer_size, workdir->dir());

        //VERIFY that number of buckets in each splitted storage the same
        for (size_t i = 0; i < unmerged_kmerfiles2.size(); ++i) {
            VERIFY(unmerged_kmerfiles2[i].num_buckets() == kpostorage.num_buckets());
        }

        INFO("Split_kpo_mers finished");

        for (unsigned i = 0; i < kmerfiles2.num_buckets(); ++i) {
            outputfiles.push_back(kmerfiles2.create(i)->file());
        }

        merge_kmer_files(std::move(unmerged_kmerfiles2), outputfiles);
        INFO("Merge_kmer_files finished");
    }
    treg.stop_listening();

    partask::broadcast(kmerfiles2);
    INFO("Total kmers= " << kmerfiles2.total_kmers());
    BuildIndexMPI(index, kmerfiles2, /* save_final */ true);

    auto fill_index = treg.add<FillIndexTask<Index, Seq>>(std::ref(index));
    treg.listen();
    DEBUG("Listening started");

    if (partask::master()) {
        fill_index(kpostorage);
    }
    treg.stop_listening();

    INFO("Building k-mer extensions from k+1-mers finished.");
}
} //utils
