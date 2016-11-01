#pragma once

#include "pipeline/graph_pack.hpp"
#include "utils/indices/perfect_hash_map_builder.hpp"

namespace debruijn_graph {

typedef uint16_t Mpl;
typedef std::size_t Offset;
static const Mpl INVALID_MPL = Mpl(-1);

typedef typename std::vector<Mpl> MplVector;
typedef typename std::vector<double> AbundanceVector;

void SetSampleCount(size_t sample_cnt);
size_t SampleCount();

class KmerProfile {

public:
    typedef Mpl value_type;

    KmerProfile(const value_type* ptr = nullptr):
        ptr_(ptr) {
    }

    KmerProfile(const MplVector& vec):
        ptr_(&vec.front()) {
    }

    size_t size() const {
        return SampleCount();
    }

    Mpl operator[](size_t i) const {
        VERIFY(i < size());
        return ptr_[i];
    }

    const value_type* begin() const {
        return ptr_;
    }

    const value_type* end() const {
        return ptr_ + size();
    }

private:
    const value_type* ptr_;
};

typedef std::vector<KmerProfile> KmerProfiles;

template<class CovVecs>
AbundanceVector MeanVector(const CovVecs& cov_vecs) {
    VERIFY(cov_vecs.size() != 0);
    size_t sample_cnt = cov_vecs.front().size();
    AbundanceVector answer(sample_cnt, 0.);

    for (const auto& cov_vec : cov_vecs) {
        for (size_t i = 0; i < sample_cnt; ++i) {
            answer[i] += double(cov_vec[i]);
        }
    }

    for (size_t i = 0; i < sample_cnt; ++i) {
        answer[i] /= double(cov_vecs.size());
    }
    return answer;
}

template<class AbVector>
std::string PrintVector(const AbVector& mpl_vector) {
    stringstream ss;
    copy(mpl_vector.begin(), mpl_vector.end(),
         ostream_iterator<typename AbVector::value_type>(ss, " "));
    return ss.str();
}

class SingleClusterAnalyzer {
    static const uint MAX_IT = 10;

    double coord_vise_proximity_;
    double central_clust_share_;

    MplVector SampleMpls(const KmerProfiles& kmer_mpls, size_t sample) const;
    Mpl SampleMedian(const KmerProfiles& kmer_mpls, size_t sample) const;
    MplVector MedianVector(const KmerProfiles& kmer_mpls) const;
    bool AreClose(const KmerProfile& c, const KmerProfile& v) const;
    KmerProfiles CloseKmerMpls(const KmerProfiles& kmer_mpls, const KmerProfile& center) const;

public:
    SingleClusterAnalyzer(double coord_vise_proximity = 0.7,
                          double central_clust_share = 0.7) :
        coord_vise_proximity_(coord_vise_proximity),
        central_clust_share_(central_clust_share) {
    }

    boost::optional<AbundanceVector> operator()(const KmerProfiles& kmer_mpls) const;

private:
    DECL_LOGGER("SingleClusterAnalyzer");
};

class ContigAbundanceCounter {
    typedef typename InvertableStoring::trivial_inverter<Offset> InverterT;

    typedef KeyStoringMap<conj_graph_pack::seq_t,
                          Offset,
                          kmer_index_traits<conj_graph_pack::seq_t>,
                          InvertableStoring> IndexT;

    unsigned k_;
    SingleClusterAnalyzer cluster_analyzer_;
    double min_earmark_share_;
    IndexT kmer_mpl_;
    InverterT inverter_;
    std::vector<Mpl> mpl_data_;

    void FillMplMap(const std::string& kmers_mpl_file);

    vector<std::string> SplitOnNs(const std::string& seq) const;

public:
    ContigAbundanceCounter(unsigned k,
                           const SingleClusterAnalyzer& cluster_analyzer,
                           const std::string& work_dir,
                           double min_earmark_share = 0.7) :
        k_(k),
        cluster_analyzer_(cluster_analyzer),
        min_earmark_share_(min_earmark_share),
        kmer_mpl_(k_, work_dir) {
    }

    void Init(const std::string& kmer_mpl_file);

    boost::optional<AbundanceVector> operator()(const std::string& s, const std::string& /*name*/ = "") const;

private:
    DECL_LOGGER("ContigAbundanceCounter");
};

}
