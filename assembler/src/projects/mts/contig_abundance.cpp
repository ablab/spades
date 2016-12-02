#include "contig_abundance.hpp"
#include "utils/indices/kmer_splitters.hpp"

namespace debruijn_graph {

size_t sample_cnt_ = 0;

void SetSampleCount(size_t sample_cnt) {
    sample_cnt_ = sample_cnt;
}

size_t SampleCount() {
    return sample_cnt_;
}

MplVector SingleClusterAnalyzer::SampleMpls(const KmerProfiles& kmer_mpls, size_t sample) const {
    MplVector answer;
    answer.reserve(kmer_mpls.size());
    for (const auto& kmer_mpl : kmer_mpls) {
        answer.push_back(kmer_mpl[sample]);
    }
    return answer;
}

Mpl SingleClusterAnalyzer::SampleMedian(const KmerProfiles& kmer_mpls, size_t sample) const {
    std::vector<Mpl> sample_mpls = SampleMpls(kmer_mpls, sample);

    std::nth_element(sample_mpls.begin(), sample_mpls.begin() + sample_mpls.size()/2, sample_mpls.end());
    return sample_mpls[sample_mpls.size()/2];
}

MplVector SingleClusterAnalyzer::MedianVector(const KmerProfiles& kmer_mpls) const {
    VERIFY(kmer_mpls.size() != 0);
    MplVector answer(SampleCount(), 0);
    for (size_t i = 0; i < SampleCount(); ++i) {
        answer[i] = SampleMedian(kmer_mpls, i);
    }
    return answer;
}

bool SingleClusterAnalyzer::AreClose(const KmerProfile& c, const KmerProfile& v) const {
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

boost::optional<AbundanceVector> SingleClusterAnalyzer::operator()(const KmerProfiles& kmer_mpls) const {
    auto med = MedianVector(kmer_mpls);
    return AbundanceVector(med.begin(), med.end());
    //return boost::optional<AbundanceVector>(answer);
    //MplVector center = MedianVector(kmer_mpls);
    //auto locality = CloseKmerMpls(kmer_mpls, KmerProfile(center));

    //for (size_t it_cnt = 0; it_cnt < MAX_IT; ++it_cnt) {
    //    DEBUG("Iteration " << it_cnt);
    //    DEBUG("Center is " << PrintVector(center));

    //    DEBUG("Locality size is " << locality.size()
    //              << " making " << (double(locality.size()) / double(kmer_mpls.size()))
    //              << " of total # points");

    //    double center_share = double(locality.size()) / double(kmer_mpls.size());
    //    if (math::ls(center_share, central_clust_share_)) {
    //        DEBUG("Detected central area contains too few k-mers: share " << center_share
    //                  << " ; center size " << locality.size()
    //                  << " ; total size " << kmer_mpls.size());
    //        return boost::none;
    //    }

    //    MplVector update = MedianVector(locality);
    //    DEBUG("Center update is " << PrintVector(update));

    //    if (center == update) {
    //        DEBUG("Old and new centers matched on iteration " << it_cnt);
    //        break;
    //    }

    //    center = update;
    //    locality = CloseKmerMpls(kmer_mpls, center);
    //}

    //return boost::optional<AbundanceVector>(MeanVector(locality, sample_cnt_));
}

vector<std::string> ContigAbundanceCounter::SplitOnNs(const std::string& seq) const {
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

void ContigAbundanceCounter::Init(const std::string& file_prefix) {
    VERIFY(SampleCount() != 0);
    INFO("Loading kmer index");
    std::ifstream kmers_in(file_prefix + ".kmm", std::ios::binary);
    kmer_mpl_.BinRead(kmers_in, file_prefix + ".kmm");

    INFO("Loading kmer profiles data");
    const size_t data_size = SampleCount() * kmer_mpl_.size();
    mpl_data_.resize(data_size);
    std::ifstream mpls_in(file_prefix + ".bpr", std::ios::binary);
    mpls_in.read((char *)&mpl_data_[0], data_size * sizeof(Mpl));
}

boost::optional<AbundanceVector> ContigAbundanceCounter::operator()(
        const std::string& s,
        const std::string& /*name*/) const {
    KmerProfiles kmer_mpls;

    for (const auto& seq : SplitOnNs(s)) {
        if (seq.size() < k_)
            continue;

        auto kwh = kmer_mpl_.ConstructKWH(RtSeq(k_, seq));
        kwh >>= 'A';

        for (size_t j = k_ - 1; j < seq.size(); ++j) {
            kwh <<= seq[j];
            TRACE("Processing kmer " << kwh.key().str());
            if (kmer_mpl_.valid(kwh)) {
                TRACE("Valid");
                KmerProfile prof(&mpl_data_[kmer_mpl_.get_value(kwh, inverter_)]);
                kmer_mpls.push_back(prof);
                //if (!name.empty()) {
                //    os << PrintVector(kmer_mpl_.get_value(kwh, inverter_), sample_cnt_) << std::endl;
                //}
                TRACE(PrintVector(prof));
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

    return cluster_analyzer_(kmer_mpls);
}

}
