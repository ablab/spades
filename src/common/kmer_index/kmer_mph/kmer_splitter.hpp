//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2020-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "kmer_buckets.hpp"

#include "adt/kmer_vector.hpp"
#include "utils/filesystem/file_limit.hpp"
#include "utils/filesystem/temporary.hpp"
#include "utils/memory_limit.hpp"
#include "utils/logger/logger.hpp"

#include <pdqsort/pdqsort_pod.h>
#include <string>
#include <cstdio>

namespace kmers {

template<class Seq>
class KMerSplitter {
public:
    typedef typename kmer::KMerSegmentPolicy<Seq> KMerBuckets;
    typedef std::vector<fs::DependentTmpFile> RawKMers;

    KMerSplitter(const std::filesystem::path &work_dir, unsigned K)
            : KMerSplitter(fs::tmp::make_temp_dir(work_dir, "kmer_splitter"), K) {}

    KMerSplitter(fs::TmpDir work_dir, unsigned K)
            : work_dir_(work_dir), K_(K) {}

    virtual ~KMerSplitter() {}

    virtual RawKMers Split(size_t num_files, unsigned nthreads) = 0;

    size_t kmer_size() const {
        return Seq::GetDataSize(K_) * sizeof(typename Seq::DataType);
    }

    unsigned K() const { return K_; }
    KMerBuckets bucket_policy() const { return bucket_; }

protected:
    fs::TmpDir work_dir_;
    unsigned K_;
    KMerBuckets bucket_;

    DECL_LOGGER("K-mer Splitting");
};

template<class Seq>
class KMerSortingSplitter : public KMerSplitter<Seq> {
public:
    using typename KMerSplitter<Seq>::RawKMers;

    KMerSortingSplitter(const std::filesystem::path &work_dir, unsigned K)
            : KMerSplitter<Seq>(work_dir, K), cell_size_(0), num_files_(0) {}

    KMerSortingSplitter(fs::TmpDir work_dir, unsigned K)
            : KMerSplitter<Seq>(work_dir, K), cell_size_(0), num_files_(0) {}

protected:
    using SeqKMerVector = adt::KMerVector<Seq>;
    using KMerBuffer = std::vector<SeqKMerVector>;

    std::vector<KMerBuffer> kmer_buffers_;
    size_t cell_size_;
    size_t num_files_;

    RawKMers PrepareBuffers(size_t num_files, unsigned nthreads, size_t reads_buffer_size) {
        num_files_ = num_files;
        this->bucket_.reset(num_files);

        // Determine the set of output files
        RawKMers out;
        auto tmp_prefix = this->work_dir_->tmp_file("kmers_raw");
        for (unsigned i = 0; i < num_files_; ++i) {
            out.emplace_back(tmp_prefix->CreateDep(std::to_string(i)));
            fclose(fopen(out.back()->file().c_str(), "ab"));
        }

        size_t file_limit = num_files_ + 2*nthreads;
        size_t res = utils::limit_file(file_limit);
        if (res < file_limit) {
            WARN("Failed to setup necessary limit for number of open files. The process might crash later on.");
            WARN("Do 'ulimit -n " << file_limit << "' in the console to overcome the limit");
        }

        if (reads_buffer_size == 0) {
            reads_buffer_size = 536870912ull;
            size_t mem_limit =  (size_t)((double)(utils::get_free_memory()) / (nthreads * 3));
            INFO("Memory available for splitting buffers: " << (double)mem_limit / 1024.0 / 1024.0 / 1024.0 << " Gb");
            reads_buffer_size = std::min(reads_buffer_size, mem_limit);
        }
        cell_size_ = reads_buffer_size / (num_files_ * this->kmer_size());
        // Set sane minimum cell size
        if (cell_size_ < 16384)
            cell_size_ = 16384;

        INFO("Using cell size of " << cell_size_);
        kmer_buffers_.resize(nthreads);
        for (unsigned i = 0; i < nthreads; ++i) {
            KMerBuffer &entry = kmer_buffers_[i];
            entry.resize(num_files_, adt::KMerVector<Seq>(this->K_, (size_t) (1.1 * (double) cell_size_)));
        }

        return out;
    }

    bool push_back_internal(const Seq &seq, unsigned thread_id) {
        VERIFY(thread_id < kmer_buffers_.size());
        KMerBuffer &entry = kmer_buffers_[thread_id];

        size_t idx = this->bucket_(seq);
        entry[idx].push_back(seq);
        return entry[idx].size() > cell_size_;
    }

    void DumpBuffers(const RawKMers &ostreams) {
        VERIFY(ostreams.size() == num_files_ && kmer_buffers_[0].size() == num_files_);

#   pragma omp parallel for
        for (size_t k = 0; k < num_files_; ++k) {
            // Below k is thread id!

            size_t sz = 0;
            for (size_t i = 0; i < kmer_buffers_.size(); ++i)
                sz += kmer_buffers_[i][k].size();

            adt::KMerVector<Seq> SortBuffer(this->K_, sz);
            for (auto & entry : kmer_buffers_) {
                const auto &buffer = entry[k];
                for (size_t j = 0; j < buffer.size(); ++j)
                    SortBuffer.push_back(buffer[j]);
            }
            pdqsort_pod(SortBuffer.data(), SortBuffer.data() + SortBuffer.size() * SortBuffer.el_size(), SortBuffer.el_size());
            auto it = std::unique(SortBuffer.begin(), SortBuffer.end(), typename adt::KMerVector<Seq>::equal_to());

#     pragma omp critical
            {
                size_t cnt =  it - SortBuffer.begin();

                // Write k-mers
                FILE *f = fopen(ostreams[k]->file().c_str(), "ab");
                if (!f)
                    FATAL_ERROR("Cannot open temporary file " << ostreams[k]->file() << " for writing");
                size_t res = fwrite(SortBuffer.data(), SortBuffer.el_data_size(), cnt, f);
                if (res != cnt)
                    FATAL_ERROR("I/O error! Incomplete write! Reason: " << strerror(errno) << ". Error code: " << errno);
                fclose(f);

                // Write index
                f = fopen((ostreams[k]->file().native() + ".idx").c_str(), "ab");
                if (!f)
                    FATAL_ERROR("Cannot open temporary file " << ostreams[k]->file() << " for writing");
                res = fwrite(&cnt, sizeof(cnt), 1, f);
                if (res != 1)
                    FATAL_ERROR("I/O error! Incomplete write! Reason: " << strerror(errno) << ". Error code: " << errno);
                fclose(f);
            }
        }

        for (auto & entry : kmer_buffers_)
            for (auto & eentry : entry)
                eentry.clear();
    }

    void ClearBuffers() {
        for (auto & entry : kmer_buffers_)
            for (auto & eentry : entry) {
                eentry.clear();
                eentry.shrink_to_fit();
            }
    }
};

}
