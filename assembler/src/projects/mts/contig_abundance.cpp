#include "contig_abundance.hpp"
#include "utils/kmer_mph/kmer_splitters.hpp"

namespace debruijn_graph {

MplVector SampleMpls(const KmerProfiles& kmer_mpls, size_t sample) {
    MplVector answer;
    answer.reserve(kmer_mpls.size());
    for (const auto& kmer_mpl : kmer_mpls) {
        answer.push_back(kmer_mpl[sample]);
    }
    return answer;
}

//------------------------------------------------------------------------------

template<typename T>
T SampleMedian(MplVector &&sample_mpls);

template<>
Abundance SampleMedian(MplVector &&sample_mpls) {
    std::nth_element(sample_mpls.begin(), sample_mpls.begin() + sample_mpls.size() / 2, sample_mpls.end());
    return (Abundance)sample_mpls[sample_mpls.size() / 2];
}

template<>
AbVar SampleMedian(MplVector &&sample_mpls) {
    auto median = SampleMedian<Abundance>(MplVector(sample_mpls));

    //Median absolute deviation
    for (auto& i: sample_mpls)
        i = (Mpl)std::abs(int(i) - int(median));

    std::nth_element(sample_mpls.begin(), sample_mpls.begin() + sample_mpls.size() / 2, sample_mpls.end());
    Var mad = (Var) sample_mpls[sample_mpls.size() / 2]; // * 1.4826; //constant scale factor for normal distibution

    return {median, mad};
}

//------------------------------------------------------------------------------

template<typename T>
T WinsoredMean(std::vector<Mpl>& v);

template<>
Abundance WinsoredMean(std::vector<Mpl>& v) {
    const float frac = 0.05f;
    size_t offset = (size_t)std::ceil(float(v.size()) * frac);
    //std::sort(v.begin(), v.end());
    std::nth_element(v.begin(), v.end() - offset, v.end());
    for (size_t i = 0; i < offset; ++i) {
        v[i] = v[offset];
        v[v.size() - i - 1] = v[v.size() - offset - 1];
    }
    size_t sum = 0;
    for (Mpl i: v)
        sum += i;
    return (Abundance)sum / (Abundance)v.size();
}

Var Variance(const std::vector<Mpl>& v, Abundance mean) {
    size_t sum = 0;
    for (Mpl i: v)
        sum += (i * i);
    return (Abundance)sum / (Abundance)v.size() - mean * mean;
}

template<>
AbVar WinsoredMean(std::vector<Mpl>& v) {
    Abundance mean = WinsoredMean<Abundance>(v);
    Var var = Variance(v, mean);
    return {mean, var};
}

std::ostream& operator<<(std::ostream& str, AbVar ab_var) {
    str << ab_var.ab << "\t" << ab_var.var;
    return str;
}

//------------------------------------------------------------------------------

template<typename T> using PointEstimator = std::function<T(MplVector &&)>;

template<typename T>
Profile<T> CountProfile(const KmerProfiles& kmer_mpls, PointEstimator<T> point_estimator) {
    VERIFY(kmer_mpls.size() != 0);
    std::vector<T> res;
    res.reserve(KmerProfileIndex::SampleCount());
    for (size_t i = 0; i < KmerProfileIndex::SampleCount(); ++i) {
        auto sample_mpls = SampleMpls(kmer_mpls, i);
        res.push_back(point_estimator(std::move(sample_mpls)));
        TRACE(sample_mpls << ": " << res.back());
    }
    return res;
}

template<typename T>
typename ClusterAnalyzer<T>::Result TrivialClusterAnalyzer<T>::operator()(const KmerProfiles& kmer_mpls) const {
    return typename ClusterAnalyzer<T>::Result(CountProfile(kmer_mpls, PointEstimator<T>(SampleMedian<T>)));
}

/*
bool SingleClusterAnalyzer::AreClose(const KmerProfileIndex::KmerProfile& c,
                                     const KmerProfileIndex::KmerProfile& v) const {
    //VERIFY(c.size() == v.size());
    double sum = 0;
    size_t non_zero_cnt = 0;
    for (size_t i = 0; i < c.size(); ++i) {
        double norm = 1.;
        if (c[i] != 0) {
            //norm = std::sqrt(double(c[i]));
            norm = double(c[i]);
            ++non_zero_cnt;
        }
        sum += std::abs(double(c[i]) - double(v[i])) / norm;
    }
    return math::ls(sum, coord_vise_proximity_ * double(non_zero_cnt));
}

KmerProfiles SingleClusterAnalyzer::CloseKmerMpls(const KmerProfiles& kmer_mpls, const KmerProfile& center) const {
    KmerProfiles answer;
    for (const auto& kmer_mpl : kmer_mpls) {
        if (AreClose(center, kmer_mpl)) {
            answer.push_back(kmer_mpl);
        } else {
            TRACE("Far kmer mpl " << PrintVector(kmer_mpl));
        }
    }
    return answer;
}

ClusterAnalyzer::Result SingleClusterAnalyzer::operator()(const KmerProfiles& kmer_mpls) const {
    MplVector center = MedianVector(kmer_mpls);
    auto locality = CloseKmerMpls(kmer_mpls, KmerProfile(center));

    for (size_t it_cnt = 0; it_cnt < MAX_IT; ++it_cnt) {
        DEBUG("Iteration " << it_cnt);
        DEBUG("Center is " << PrintVector(center));

        DEBUG("Locality size is " << locality.size()
                  << " making " << (double(locality.size()) / double(kmer_mpls.size()))
                  << " of total # points");

        double center_share = double(locality.size()) / double(kmer_mpls.size());
        if (math::ls(center_share, central_clust_share_)) {
            DEBUG("Detected central area contains too few k-mers: share " << center_share
                      << " ; center size " << locality.size()
                      << " ; total size " << kmer_mpls.size());
            return ClusterAnalyzer::Result();
        }

        MplVector update = MedianVector(locality);
        DEBUG("Center update is " << PrintVector(update));

        if (center == update) {
            DEBUG("Old and new centers matched on iteration " << it_cnt);
            break;
        }

        center = update;
        locality = CloseKmerMpls(kmer_mpls, center);
    }

    return ClusterAnalyzer::Result(MeanVector(locality));
}*/

vector<std::string> SplitOnNs(const std::string& seq) {
    vector<std::string> answer;
    for (size_t i = 0; i < seq.size(); i++) {
        size_t j = i;
        while (j < seq.size() && is_nucl(seq[j])) {
            j++;
        }
        if (j > i) {
            answer.push_back(seq.substr(i, j - i));
            i = j;
        }
    }
    return answer;
}

size_t KmerProfileIndex::sample_cnt_ = 0;

KmerProfileIndex::KmerProfileIndex(unsigned k,
                                   const std::string& index_prefix):
    index_(k) {
    std::string index_file = index_prefix + ".kmm";
    INFO("Loading kmer index from " << index_file);
    std::ifstream kmers_in(index_file, std::ios::binary);
    index_.BinRead(kmers_in, index_file);

    const size_t data_size = SampleCount() * index_.size();
    std::string profiles_file = index_prefix + ".bpr";
    INFO("Loading profiles data of " << data_size << " elements from " << profiles_file);
    profiles_ = ProfilesT(profiles_file, data_size, false);
    //Prefetch mmapped data by force
    Mpl checksum = 0;
    for (const Mpl* ptr = profiles_->data(); ptr < profiles_->data() + data_size; ptr += 4096 / sizeof(Mpl))
        checksum = Mpl(checksum + *ptr);
    INFO("Kmer index loaded; checksum: " << checksum);
}

KmerProfileIndex::KmerProfileIndex(KmerProfileIndex&& other):
    inverter_(std::move(other.inverter_)),
    index_(std::move(other.index_)),
    profiles_(std::move(other.profiles_)) {
}

KmerProfileIndex::KeyWithHash KmerProfileIndex::Construct(const KmerProfileIndex::KeyType& key) const {
    return index_.ConstructKWH(key);
}

boost::optional<KmerProfile> KmerProfileIndex::operator[](const KmerProfileIndex::KeyWithHash& kwh) const {
    if (index_.valid(kwh)) {
        return KmerProfile(profiles_->data() + index_.get_value(kwh, inverter_));
    } else
        return boost::none;
}

template<typename T>
ProfileCounter<T>::ProfileCounter(unsigned k,
                                  std::unique_ptr<ClusterAnalyzer<T>> cluster_analyzer,
                                  const std::string& index_prefix,
                                  double min_earmark_share) :
    k_(k),
    cluster_analyzer_(std::move(cluster_analyzer)),
    min_earmark_share_(min_earmark_share),
    profile_index_(k_, index_prefix) {
}

template<typename T>
ProfileCounter<T>::ProfileCounter(ProfileCounter&& other):
    k_(other.k_),
    cluster_analyzer_(std::move(other.cluster_analyzer_)),
    min_earmark_share_(other.min_earmark_share_),
    profile_index_(std::move(other.profile_index_)) {
}

template<typename T>
typename ClusterAnalyzer<T>::Result ProfileCounter<T>::operator()(
        const std::string& s,
        const std::string& /*name*/) const {
    KmerProfiles kmer_mpls;

    for (const auto& seq : SplitOnNs(s)) {
        if (seq.size() < k_)
            continue;

        auto kwh = profile_index_.Construct(RtSeq(k_, seq));
        kwh >>= 'A';

        for (size_t j = k_ - 1; j < seq.size(); ++j) {
            kwh <<= seq[j];
            TRACE("Processing kmer " << kwh.key().str());
            auto prof = profile_index_[kwh];
            if (prof) {
                TRACE("Valid");
                kmer_mpls.push_back(*prof);
                //if (!name.empty()) {
                //    os << PrintVector(kmer_mpl_.get_value(kwh, inverter_), sample_cnt_) << std::endl;
                //}
                TRACE(PrintVector(*prof));
            } else {
                TRACE("Invalid");
            }
        }
    }

    double earmark_share = double(kmer_mpls.size()) / double(s.size() - k_ + 1);
    DEBUG("Earmark k-mers: share " << earmark_share
              << " # earmarks " << kmer_mpls.size()
              << " ; total # " << (s.size() - k_ + 1));
    if (math::ls(earmark_share, min_earmark_share_)) {
        DEBUG("Too few earmarks");
        return boost::none;
    }

    return (*cluster_analyzer_)(kmer_mpls);
}

//Explicit instantiations for what we can estimate
#define INSTANTIATE(type) \
template class TrivialClusterAnalyzer<type>; \
template class ProfileCounter<type>;

INSTANTIATE(Abundance);
INSTANTIATE(AbVar);

}
