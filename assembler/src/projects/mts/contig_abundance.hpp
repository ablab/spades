#pragma once

namespace debruijn_graph {

static const uint INVALID_MPL = uint(-1);
//TODO get rid of MAX_SAMPLE_CNT
static const uint MAX_SAMPLE_CNT = 20;
typedef uint Mpl;
typedef typename std::array<Mpl, MAX_SAMPLE_CNT> MplVector;
typedef typename std::array<double, MAX_SAMPLE_CNT> AbundanceVector;

template<class CovVecs>
AbundanceVector MeanVector(const CovVecs& cov_vecs, size_t sample_cnt) {
    AbundanceVector answer;

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
std::string PrintVector(const AbVector& mpl_vector, size_t sample_cnt) {
    stringstream ss;
    copy(mpl_vector.begin(), mpl_vector.begin() + sample_cnt,
         ostream_iterator<typename AbVector::value_type>(ss, " "));
    return ss.str();
}

class SingleClusterAnalyzer {
    static const uint MAX_IT = 10;

    size_t sample_cnt_;
    double coord_vise_proximity_;
    double central_clust_share_;

    std::vector<Mpl> SampleMpls(const std::vector<MplVector>& kmer_mpls, size_t sample) const {
        std::vector<Mpl> answer;
        answer.reserve(kmer_mpls.size());
        for (const auto& kmer_mpl : kmer_mpls) {
            answer.push_back(kmer_mpl[sample]);
        }
        return answer;
    }

    Mpl SampleMedian(const std::vector<MplVector>& kmer_mpls, size_t sample) const {
        std::vector<Mpl> sample_mpls = SampleMpls(kmer_mpls, sample);

        std::nth_element(sample_mpls.begin(), sample_mpls.begin() + sample_mpls.size()/2, sample_mpls.end());
        return sample_mpls[sample_mpls.size()/2];
    }

    MplVector MedianVector(const std::vector<MplVector>& kmer_mpls) const {
        MplVector answer;
        for (size_t i = 0; i < sample_cnt_; ++i) {
            answer[i] = SampleMedian(kmer_mpls, i);
        }
        return answer;
    }

    bool AreClose(const MplVector& c, const MplVector& v) const {
        double sum = 0;
        size_t non_zero_cnt = 0;
        for (size_t i = 0; i < sample_cnt_; ++i) {
            double norm = 1.;
            if (c[i] != 0) {
                norm = std::sqrt(double(c[i]));
                ++non_zero_cnt;
            }
            sum += std::abs(double(c[i]) - double(v[i])) / norm;
        }
        return math::ls(sum, coord_vise_proximity_ * double(non_zero_cnt));
    }

    std::vector<MplVector> CloseKmerMpls(const std::vector<MplVector>& kmer_mpls, const MplVector& center) const {
        std::vector<MplVector> answer;
        for (const auto& kmer_mpl : kmer_mpls) {
            if (AreClose(center, kmer_mpl)) {
                answer.push_back(kmer_mpl);
            } else {
                TRACE("Far kmer mpl " << PrintVector(kmer_mpl, sample_cnt_));
            }
        }
        return answer;
    }

public:
    SingleClusterAnalyzer(size_t sample_cnt,
                          double coord_vise_proximity = 2.,
                          double central_clust_share = 0.7) :
        sample_cnt_(sample_cnt),
        coord_vise_proximity_(coord_vise_proximity),
        central_clust_share_(central_clust_share) {

    }

    boost::optional<AbundanceVector> operator()(const std::vector<MplVector>& kmer_mpls) const {
        MplVector center = MedianVector(kmer_mpls);
        auto locality = CloseKmerMpls(kmer_mpls, center);

        for (size_t it_cnt = 0; it_cnt < MAX_IT; ++it_cnt) {
            DEBUG("Iteration " << it_cnt);
            DEBUG("Center is " << PrintVector(center, sample_cnt_));

            DEBUG("Locality size is " << locality.size()
                      << " making " << (double(locality.size()) / double(kmer_mpls.size()))
                      << " of total # points");

            double center_share = double(locality.size()) / double(kmer_mpls.size());
            if (math::ls(center_share, central_clust_share_)) {
                DEBUG("Detected central area contains too few k-mers: share " << center_share
                          << " ; center size " << locality.size()
                          << " ; total size " << kmer_mpls.size());
                return boost::none;
            }

            MplVector update = MedianVector(locality);
            DEBUG("Center update is " << PrintVector(update, sample_cnt_));

            if (center == update) {
                DEBUG("Old and new centers matched on iteration " << it_cnt);
                break;
            }

            center = update;
            locality = CloseKmerMpls(kmer_mpls, center);
        }

        return boost::optional<AbundanceVector>(MeanVector(locality, sample_cnt_));
    }

private:
    DECL_LOGGER("SingleClusterAnalyzer");
};

class ContigAbundanceCounter {
    typedef typename InvertableStoring::immutant_inverter<MplVector> InverterT;

    unsigned k_;
    size_t sample_cnt_;
    SingleClusterAnalyzer cluster_analyzer_;
    double min_earmark_share_;

    KeyStoringMap<conj_graph_pack::seq_t,
    MplVector,
    kmer_index_traits<conj_graph_pack::seq_t>,
    InvertableStoring> kmer_mpl_;

    InverterT inverter_;

    void FillMplMap(const std::string& kmers_mpl_file) {
        //INFO("Fill start");
        for (auto it = kmer_mpl_.value_begin(); it != kmer_mpl_.value_end(); ++it) {
            it->fill(INVALID_MPL);
        }
        std::ifstream kmers_in(kmers_mpl_file + ".kmer", std::ios::binary);
        //std::ifstream kmers_mpl_in(kmers_mpl_file + ".mpl", std::ios::binary);
        std::ifstream kmers_mpl_in(kmers_mpl_file + ".mpl");
        while (true) {
//            std::string kmer_str;
//            kmers_in >> kmer_str;
            runtime_k::RtSeq kmer(k_);
            kmer.BinRead(kmers_in);
            if (kmers_in.fail()) {
                break;
            }

//            VERIFY(kmer_str.length() == k_);
//            conj_graph_pack::seq_t kmer(k_, kmer_str.c_str());
//            kmer = gp_.kmer_mapper.Substitute(kmer);

            MplVector mpls;
            mpls.fill(INVALID_MPL);
            for (size_t i = 0; i < sample_cnt_; ++i) {
//                kmers_mpl_in.read((char*) &kmer_mpl[i], sizeof(uint));
                kmers_mpl_in >> mpls[i];
                VERIFY(!kmers_mpl_in.fail());
            }

            auto kwh = kmer_mpl_.ConstructKWH(kmer);
            VERIFY(kmer_mpl_.valid(kwh));
            kmer_mpl_.put_value(kwh, mpls, inverter_);
        }

        //check that everything was filled
        for (auto it = kmer_mpl_.value_begin(); it != kmer_mpl_.value_end(); ++it) {
            for (size_t i = 0; i < sample_cnt_; ++i) {
                VERIFY((*it)[i] != INVALID_MPL);
            }
        }
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
        sample_cnt_(sample_cnt),
        cluster_analyzer_(cluster_analyzer),
        min_earmark_share_(min_earmark_share),
        kmer_mpl_(k_, work_dir) {
    }

    void Init(const std::string& kmer_mpl_file,
              size_t /*fixme some buffer size*/read_buffer_size = 0) {
        INFO("Filling kmer multiplicities.");
        DeBruijnKMerKMerSplitter<StoringTypeFilter<InvertableStoring>>
            splitter(kmer_mpl_.workdir(), k_, k_, true, read_buffer_size);

        splitter.AddKMers(kmer_mpl_file + ".kmer");

        KMerDiskCounter<runtime_k::RtSeq> counter(kmer_mpl_.workdir(), splitter);

        kmer_mpl_.BuildIndex(counter, 16, /*nthreads*/ 1);

        FillMplMap(kmer_mpl_file);
    }

    boost::optional<AbundanceVector> operator()(const std::string& s) const {
        std::vector<MplVector> kmer_mpls;

        for (const auto& seq : SplitOnNs(s)) {
            if (seq.size() < k_)
                continue;

            auto kwh = kmer_mpl_.ConstructKWH(runtime_k::RtSeq(k_, seq));
            kwh >>= 'A';

            for (size_t j = k_ - 1; j < seq.size(); ++j) {
                kwh <<= seq[j];
                if (kmer_mpl_.valid(kwh)) {
                    kmer_mpls.push_back(kmer_mpl_.get_value(kwh, inverter_));
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
