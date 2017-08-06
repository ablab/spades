#pragma once

#include "pipeline/graph_pack.hpp"
#include "utils/ph_map/perfect_hash_map_builder.hpp"
#include "io/kmers/mmapped_reader.hpp"

namespace debruijn_graph {

typedef uint16_t Mpl;
typedef std::size_t Offset;
static const Mpl INVALID_MPL = Mpl(-1);

typedef typename std::vector<Mpl> MplVector;
typedef typename std::vector<double> AbundanceVector;

class KmerProfileIndex {
private:
    typedef typename utils::InvertableStoring::trivial_inverter<Offset> InverterT;
    typedef utils::KeyStoringMap<conj_graph_pack::seq_t,
        Offset,
        utils::kmer_index_traits<conj_graph_pack::seq_t>,
        utils::InvertableStoring> IndexT;

    static size_t sample_cnt_;

public:
    static void SetSampleCount(size_t sample_cnt) { sample_cnt_ = sample_cnt; }
    static size_t SampleCount() { return sample_cnt_; }

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

    typedef typename IndexT::KeyType KeyType;
    typedef typename IndexT::KeyWithHash KeyWithHash;

    KmerProfileIndex(unsigned k, const std::string& index_prefix, const std::string& work_dir);
    KeyWithHash Construct(const KeyType& seq) const;
    boost::optional<KmerProfile> operator[](const KeyWithHash& kwh) const;

private:
    typedef MMappedRecordArrayReader<Mpl> ProfilesT;
    InverterT inverter_;
    IndexT index_;
    boost::optional<ProfilesT> profiles_;
};

using KmerProfile = KmerProfileIndex::KmerProfile;
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

MplVector SampleMpls(const KmerProfiles& kmer_mpls, size_t sample);
Mpl SampleMedian(const KmerProfiles& kmer_mpls, size_t sample);
MplVector MedianVector(const KmerProfiles& kmer_mpls);

class ClusterAnalyzer {
public:
    virtual boost::optional<AbundanceVector> operator()(const KmerProfiles& kmer_mpls) const = 0;
    virtual ~ClusterAnalyzer() {};
};

class TrivialClusterAnalyzer : public ClusterAnalyzer {
public:
    TrivialClusterAnalyzer() {}

    boost::optional<AbundanceVector> operator()(const KmerProfiles& kmer_mpls) const override;

private:
    DECL_LOGGER("TrivialClusterAnalyzer");
};

class SingleClusterAnalyzer : public ClusterAnalyzer {
    static const uint MAX_IT = 10;

    double coord_vise_proximity_;
    double central_clust_share_;

    bool AreClose(const KmerProfile& c, const KmerProfile& v) const;
    KmerProfiles CloseKmerMpls(const KmerProfiles& kmer_mpls, const KmerProfile& center) const;

public:
    SingleClusterAnalyzer(double coord_vise_proximity = 0.7,
                          double central_clust_share = 0.7) :
        coord_vise_proximity_(coord_vise_proximity),
        central_clust_share_(central_clust_share) {
    }

    boost::optional<AbundanceVector> operator()(const KmerProfiles& kmer_mpls) const override;

private:
    DECL_LOGGER("SingleClusterAnalyzer");
};

class ContigAbundanceCounter {
    unsigned k_;
    std::unique_ptr<ClusterAnalyzer> cluster_analyzer_;
    double min_earmark_share_;

    vector<std::string> SplitOnNs(const std::string& seq) const;

    KmerProfileIndex profile_index_;

public:
    ContigAbundanceCounter(unsigned k,
                           std::unique_ptr<ClusterAnalyzer> cluster_analyzer,
                           const std::string& index_prefix,
                           const std::string& work_dir,
                           double min_earmark_share = 0.7) :
        k_(k),
        cluster_analyzer_(cluster_analyzer),
        min_earmark_share_(min_earmark_share),
        kmer_mpl_(k_) {
    }

    void Init(const std::string& kmer_mpl_file);

    boost::optional<AbundanceVector> operator()(const std::string& s, const std::string& /*name*/ = "") const;

private:
    DECL_LOGGER("ContigAbundanceCounter");
};

}
