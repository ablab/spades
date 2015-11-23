#include "standard_base.hpp"
#include "graph_pack.hpp"
#include "sequence_mapper.hpp"
#include "io/io_helper.hpp"
#include <regex>

namespace debruijn_graph {

typedef string bin_id;
typedef string contig_id;
typedef pair<contig_id, vector<bin_id>> ContigAnnotation;

inline contig_id GetId(const io::SingleRead& contig) {
     std::smatch m;
     std::regex e ("ID_(\\d+)$");
     bool success = std::regex_search(contig.name(), m, e);
     VERIFY(success);
     return m[1];
}

class AnnotationStream {
    std::ifstream inner_stream_;

    ContigAnnotation Parse(const std::string& s) const {
        ContigAnnotation annotation;
        stringstream ss(s);
        ss >> annotation.first;
        string delim;
        ss >> delim;
        VERIFY(delim == ":");
        while (true) {
            bin_id bin;
            ss >> bin;
            if (ss.fail())
                break;
            annotation.second.push_back(bin);
        }
        return annotation;
    }

public:

    AnnotationStream(const std::string& fn) : inner_stream_(fn) {

    }

    bool eof() const {
        return inner_stream_.eof();
    }

    AnnotationStream& operator >>(ContigAnnotation& annotation) {
        VERIFY(!inner_stream_.eof())
        std::string s;
        inner_stream_ >> s;

        annotation = Parse(s);
        return *this;
    }

    void close() {
        inner_stream_.close();
    }
};

class AnnotationOutStream {
    std::ofstream inner_stream_;
public:

    AnnotationOutStream(const std::string& fn) : inner_stream_(fn) {
    }

    AnnotationOutStream& operator <<(const ContigAnnotation& annotation) {
        inner_stream_ << annotation.first;
        string delim = " : ";
        for (bin_id bin : annotation.second) {
            inner_stream_ << delim << bin;
            delim = " ";
        }
        return *this;
    }

    void close() {
        inner_stream_.close();
    }
};

class EdgeAnnotation {
    const conj_graph_pack& gp_;
    set<bin_id> bins_of_interest_;
    shared_ptr<SequenceMapper<Graph>> mapper_;
    map<EdgeId, set<bin_id>> edge_annotation_;

    vector<bin_id> FilterInteresting(const vector<bin_id>& bins) const {
        vector<bin_id> answer;
        for (bin_id bin : bins) {
            if (bins_of_interest_.count(bin)) {
                answer.push_back(bin);
            }
        }
        return answer;
    }

    template<class BinCollection>
    void InnerStickAnnotation(EdgeId e, const BinCollection& bins) {
        auto& annotation = edge_annotation_[e];
        for (bin_id bin : bins) {
            annotation.insert(bin);
        }
    }

public:

    EdgeAnnotation(const conj_graph_pack& gp,
                   const vector<bin_id>& bins_of_interest) :
                       gp_(gp),
                       bins_of_interest_(bins_of_interest.begin(), bins_of_interest.end()) {
    }

    template<class BinCollection>
    void StickAnnotation(EdgeId e, const BinCollection& bins) {
        InnerStickAnnotation(e, bins);
        InnerStickAnnotation(gp_.g.conjugate(e), bins);
    }

    void StickAnnotation(EdgeId e, const bin_id& bin) {
        StickAnnotation(e, vector<bin_id>{bin});
    }

    template<class EdgeCollection>
    void StickAnnotation(const EdgeCollection& edges, const bin_id& bin) {
        for (EdgeId e : edges) {
            StickAnnotation(e, bin);
        }
    }

    void Fill(io::SingleStream& contigs, AnnotationStream& annotation_stream) {
        set<bin_id> all_bins;

        map<contig_id, vector<bin_id>> annotation_map;
        ContigAnnotation contig_annotation;
        while (!annotation_stream.eof()) {
            annotation_stream >> contig_annotation;
            auto bins = contig_annotation.second;
            if (!bins_of_interest_.empty()) {
                bins = FilterInteresting(bins);
            } else {
                insert_all(all_bins, bins);
            }
            if (!bins.empty()) {
                annotation_map[contig_annotation.first] = bins;
            }
        }

        io::SingleRead contig;
        while (!contigs.eof()) {
            contigs >> contig;
            contig_id id = GetId(contig);
            vector<bin_id> bins;
            if (annotation_map.count(id)) {
                bins = annotation_map[id];
            }
            if (!bins.empty()) {
                for (EdgeId e : mapper_->MapRead(contig).simple_path()) {
                    StickAnnotation(e, bins);
                }
            }
        }

        if (bins_of_interest_.empty()) {
            bins_of_interest_ = all_bins;
        }
    }

    vector<bin_id> Annotation(EdgeId e) const {
        if (!edge_annotation_.count(e)) {
            return {};
        }
        const auto& annotation = get(edge_annotation_, e);
        return vector<bin_id>(annotation.begin(), annotation.end());
    }

    set<bin_id> RelevantBins(const io::SingleRead& r) const {
        set<bin_id> answer;
        for (EdgeId e : mapper_->MapRead(r).simple_path()) {
            insert_all(answer, Annotation(e));
        }
        return answer;
    }

    set<EdgeId> EdgesOfBin(bin_id bin) const {
        set<EdgeId> answer;
        for (auto ann_pair : edge_annotation_) {
            if (ann_pair.second.count(bin)) {
                answer.insert(ann_pair.first);
            }
        }
        return answer;
    }

    const set<bin_id>& interesting_bins() const {
        return bins_of_interest_;
    }

};

}
