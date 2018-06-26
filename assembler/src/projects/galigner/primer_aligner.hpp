#pragma once

#include "rna16S_config.hpp"

namespace rna16S_mapping {

struct ReadMapping {
    void ExtractPositions() {
        start_pos_ = -1;
        end_pos_ = -1;
        string name = read_;
        int open_bracket_ind = -1;
        int close_bracket_ind = -1;
        int comma_ind = -1;
        for (int i = 0; i < name.size(); ++ i) {
            if (name[i] == '[') {
                open_bracket_ind = i;
            }
            if (name[i] == ']') {
                close_bracket_ind = i;
            }
            if (name[i] == ',') {
                comma_ind = i;
            }
        }
        start_pos_ = std::stoi(name.substr(open_bracket_ind + 1, comma_ind - open_bracket_ind));
        end_pos_ = std::stoi(name.substr(comma_ind + 1, close_bracket_ind - comma_ind));
        INFO("Name=" << name << " start_pos_=" << start_pos_ << " end_pos_=" << end_pos_);
    }

    ReadMapping(const std::string &read, const std::string &read_seq, const int &v, const std::vector<EdgeId> &e, const std::vector<MappingRange> &range)
        : read_(read), read_seq_(read_seq), v_(v), e_(e), range_(range) {
        ExtractPositions();
    }

    bool operator < (const ReadMapping &mapping) const {
        return (this->start_pos_ < mapping.start_pos_ || (this->start_pos_ == mapping.start_pos_ && this->end_pos_ < mapping.start_pos_) );
    }

    std::string read_;
    std::string read_seq_;
    int v_;
    std::vector<EdgeId> e_;
    std::vector<MappingRange> range_;
    int start_pos_;
    int end_pos_;
};

class PrimerAligner {
public:
    PrimerAligner(debruijn_graph::conj_graph_pack &gp, RnaAlignerConfig &cfg): gp_(gp), cfg_(cfg) {}
    void PreparePrimers(std::vector<io::SingleRead> &wrappedprimers, int threads);

    std::vector<ReadMapping>& primers() { return primers_; }

    std::map<VertexId, std::set<VertexId>>& dist() { return dist_; }

private:
    void AlignPrimer(const io::SingleRead &read, int &v, std::vector<EdgeId> &e, std::vector<MappingRange> &range);
    void FormDistanceMatrix();

    const debruijn_graph::conj_graph_pack &gp_;
    const rna16S_mapping::RnaAlignerConfig cfg_;
    
    std::vector<ReadMapping> primers_;
    std::map<VertexId, std::set<VertexId> > dist_;
};


} //namespace rna16S_mapping