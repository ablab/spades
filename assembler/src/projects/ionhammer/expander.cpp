//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "expander.hpp"

#include "config_struct.hpp"
#include "kmer_data.hpp"
#include "valid_hkmer_generator.hpp"

#include "io/reads/file_reader.hpp"

#include <vector>
#include <cstring>

bool Expander::operator()(const io::SingleRead &r) {
    size_t sz = r.size();

    std::vector<unsigned> covered_by_solid(sz, false);
    std::vector<size_t> kmer_indices(sz, -1ull);

    ValidHKMerGenerator<hammer::K> gen(r);
    while (gen.HasMore()) {
        hammer::HKMer kmer = gen.kmer();
        size_t idx = data_.seq_idx(kmer), kl = kmer.size();
        size_t read_pos = gen.pos() - kl;

        kmer_indices[read_pos] = idx;
        if (data_[idx].changeto == idx &&
            data_[idx].qual < cfg::get().center_qual_threshold) {
            for (size_t j = read_pos; j < read_pos + kl; ++j) {
                VERIFY_MSG(j < sz, "read_pos == " << read_pos << ", r.size() == " << r.size() << ", kmer: " << kmer << ", read: " << r.GetSequenceString());
                covered_by_solid[j] = true;
            }
        }

        gen.Next();
    }

    for (size_t j = 0; j < sz; ++j) {
        if (!covered_by_solid[j] || kmer_indices[j] == -1ull)
            continue;

        size_t idx = kmer_indices[j];
        auto &kmer_data = data_[idx];
        if (kmer_data.changeto != idx) {
#           pragma omp atomic
            changed_ += 1;

            kmer_data.lock();
            kmer_data.changeto = static_cast<unsigned>(idx);
            kmer_data.unlock();
        }
    }

    return false;
}
