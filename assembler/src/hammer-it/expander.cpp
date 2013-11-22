#include "expander.hpp"

#include "config_struct.hpp"
#include "kmer_data.hpp"
#include "valid_hkmer_generator.hpp"

#include "io/reader.hpp"

#include <vector>
#include <cstring>

bool Expander::operator()(const io::SingleRead &r) {
    size_t sz = r.size();

    std::vector<unsigned> covered_by_solid(sz, false);
    std::vector<size_t> kmer_indices(sz, -1ull);

    ValidHKMerGenerator<hammer::K> gen(r);
    while (gen.HasMore()) {
        hammer::HKMer kmer = gen.kmer();
        size_t idx = data_.seq_idx(kmer);
        size_t read_pos = gen.pos() - gen.nlen();

        kmer_indices[read_pos] = idx;
        if (data_[idx].changeto == idx &&
            data_[idx].qual < cfg::get().center_qual_threshold) {
            for (size_t j = read_pos, l = kmer.size(); j < read_pos + l; ++j)
                covered_by_solid[j] = true;
        }

        gen.Next();
    }

    for (size_t j = 0; j < sz; ++j) {
        if (!covered_by_solid[j] || kmer_indices[j] == -1ull)
            continue;

        // FIXME: Do not lock everything
        size_t idx = kmer_indices[j];
        auto &kmer_data = data_[idx];
        if (kmer_data.changeto != idx) {
#     pragma omp atomic
            changed_ += 1;

            kmer_data.lock();
            kmer_data.changeto = idx;
            kmer_data.unlock();
        }
    }

    return false;
}
