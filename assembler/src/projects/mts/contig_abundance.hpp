#pragma once

namespace debruijn_graph {

typedef uint16_t Mpl;
static const Mpl INVALID_MPL = Mpl(-1);

typedef typename std::vector<Mpl> MplVector;
typedef typename std::vector<double> AbundanceVector;

namespace settings {
    size_t sample_cnt = 0;
}

void SetSampleCount(size_t sample_cnt) {
    settings::sample_cnt = sample_cnt;
}

size_t SampleCount() {
    return settings::sample_cnt;
}

class KmerProfile {

public:
    typedef Mpl value_type;

    KmerProfile(Mpl* ptr = nullptr):
        ptr_(ptr)
    {}

    size_t size() const {
        return SampleCount();
    }

    Mpl operator[](size_t i) const {
        return ptr_[i];
    }

    Mpl* begin() {
        return ptr_;
    }

    const Mpl* begin() const {
        return ptr_;
    }

    Mpl* end() {
        return ptr_ + size();
    }

    const Mpl* end() const {
        return ptr_ + size();
    }

private:
    Mpl* ptr_;
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

    MplVector SampleMpls(const KmerProfiles& kmer_mpls, size_t sample) const {
        MplVector answer;
        answer.reserve(kmer_mpls.size());
        for (const auto& kmer_mpl : kmer_mpls) {
            answer.push_back(kmer_mpl[sample]);
        }
        return answer;
    }

    Mpl SampleMedian(const KmerProfiles& kmer_mpls, size_t sample) const {
        std::vector<Mpl> sample_mpls = SampleMpls(kmer_mpls, sample);

        std::nth_element(sample_mpls.begin(), sample_mpls.begin() + sample_mpls.size()/2, sample_mpls.end());
        return sample_mpls[sample_mpls.size()/2];
    }

    MplVector MedianVector(const KmerProfiles& kmer_mpls) const {
        VERIFY(kmer_mpls.size() != 0);
        size_t sample_cnt = kmer_mpls.front().size();
        MplVector answer(sample_cnt, 0);
        for (size_t i = 0; i < sample_cnt; ++i) {
            answer[i] = SampleMedian(kmer_mpls, i);
        }
        return answer;
    }

    bool AreClose(const KmerProfile& c, const KmerProfile& v) const {
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

    KmerProfiles CloseKmerMpls(const KmerProfiles& kmer_mpls, const KmerProfile& center) const {
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

public:
    SingleClusterAnalyzer(double coord_vise_proximity = 0.7,
                          double central_clust_share = 0.7) :
        coord_vise_proximity_(coord_vise_proximity),
        central_clust_share_(central_clust_share) {
    }

    boost::optional<AbundanceVector> operator()(const KmerProfiles& kmer_mpls) const {
        auto med = MedianVector(kmer_mpls);
        return AbundanceVector(med.begin(), med.end());
        //return boost::optional<AbundanceVector>(answer);
        //MplVector center = MedianVector(kmer_mpls);
        //auto locality = CloseKmerMpls(kmer_mpls, center);

        //for (size_t it_cnt = 0; it_cnt < MAX_IT; ++it_cnt) {
        //    DEBUG("Iteration " << it_cnt);
        //    DEBUG("Center is " << PrintVector(center, sample_cnt_));

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
        //    DEBUG("Center update is " << PrintVector(update, sample_cnt_));

        //    if (center == update) {
        //        DEBUG("Old and new centers matched on iteration " << it_cnt);
        //        break;
        //    }

        //    center = update;
        //    locality = CloseKmerMpls(kmer_mpls, center);
        //}

        //return boost::optional<AbundanceVector>(MeanVector(locality, sample_cnt_));
    }

private:
    DECL_LOGGER("SingleClusterAnalyzer");
};

class ContigAbundanceCounter {
    typedef typename InvertableStoring::immutant_inverter<KmerProfile> InverterT;

    unsigned k_;
    SingleClusterAnalyzer cluster_analyzer_;
    double min_earmark_share_;

    KeyStoringMap<conj_graph_pack::seq_t,
                  KmerProfile,
                  kmer_index_traits<conj_graph_pack::seq_t>,
                  InvertableStoring>
        kmer_mpl_;
    std::vector<Mpl> mpl_data_;

    InverterT inverter_;

    void FillMplMap(const std::string& kmers_mpl_file) {
        INFO("Kmer profile fill start");
        //We must allocate the whole buffer for all profiles at once
        //to avoid pointer invalidation after possible vector resize
        const size_t data_size = SampleCount() * kmer_mpl_.size();
        mpl_data_.reserve(data_size);
        INFO("Allocated buffer of " << data_size << " elements");
        std::ifstream kmers_in(kmers_mpl_file + ".kmer", std::ios::binary);
        std::ifstream kmers_mpl_in(kmers_mpl_file + ".mpl");
        while (true) {
            runtime_k::RtSeq kmer(k_);
            kmer.BinRead(kmers_in);
            if (kmers_in.fail()) {
                break;
            }

//            VERIFY(kmer_str.length() == k_);
//            conj_graph_pack::seq_t kmer(k_, kmer_str.c_str());
//            kmer = gp_.kmer_mapper.Substitute(kmer);

            KmerProfile mpls(&*mpl_data_.end());
            for (size_t i = 0; i < SampleCount(); ++i) {
                Mpl mpl;
                kmers_mpl_in >> mpl;
                VERIFY(!kmers_mpl_in.fail());
                mpl_data_.push_back(mpl);
            }
            //Double-check we haven't invalidated vector views
            VERIFY(mpl_data_.size() <= data_size);

            auto kwh = kmer_mpl_.ConstructKWH(kmer);
            VERIFY(kmer_mpl_.valid(kwh));
            kmer_mpl_.put_value(kwh, mpls, inverter_);
        }
        INFO("Kmer profile fill finish");
    }

    vector<std::string> SplitOnNs(const std::string& seq) const {
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

public:
    ContigAbundanceCounter(unsigned k,
                           size_t sample_cnt,
                           const SingleClusterAnalyzer& cluster_analyzer,
                           const std::string& work_dir,
                           double min_earmark_share = 0.7) :
        k_(k),
        cluster_analyzer_(cluster_analyzer),
        min_earmark_share_(min_earmark_share),
        kmer_mpl_(k_, work_dir) {
        SetSampleCount(sample_cnt);
    }

    void Init(const std::string& kmer_mpl_file,
              size_t /*fixme some buffer size*/read_buffer_size = 0) {
        INFO("Initializing kmer profile index");
        DeBruijnKMerKMerSplitter<StoringTypeFilter<InvertableStoring>>
            splitter(kmer_mpl_.workdir(), k_, k_, true, read_buffer_size);

        splitter.AddKMers(kmer_mpl_file + ".kmer");

        KMerDiskCounter<runtime_k::RtSeq> counter(kmer_mpl_.workdir(), splitter);

        kmer_mpl_.BuildIndex(counter, 16, /*nthreads*/ 1);

        FillMplMap(kmer_mpl_file);
    }

    boost::optional<AbundanceVector> operator()(const std::string& s, const std::string& /*name*/ = "") const {
        KmerProfiles kmer_mpls;

        for (const auto& seq : SplitOnNs(s)) {
            if (seq.size() < k_)
                continue;

            auto kwh = kmer_mpl_.ConstructKWH(runtime_k::RtSeq(k_, seq));
            kwh >>= 'A';

            for (size_t j = k_ - 1; j < seq.size(); ++j) {
                kwh <<= seq[j];
                TRACE("Processing kmer " << kwh.key().str());
                if (kmer_mpl_.valid(kwh)) {
                    TRACE("Valid");
                    kmer_mpls.push_back(kmer_mpl_.get_value(kwh, inverter_));
                    //if (!name.empty()) {
                    //    os << PrintVector(kmer_mpl_.get_value(kwh, inverter_), sample_cnt_) << std::endl;
                    //}
                    TRACE(PrintVector(kmer_mpl_.get_value(kwh, inverter_)));
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

private:
    DECL_LOGGER("ContigAbundanceCounter");
};

}
