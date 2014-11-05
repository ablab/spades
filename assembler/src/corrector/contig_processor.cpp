#include "contig_processor.hpp"
#include "config_struct.hpp"
#include "variants_table.hpp"

#include "io/ireader.hpp"
#include "io/osequencestream.hpp"
#include "io/file_reader.hpp"
#include "io/single_read.hpp"
#include "path_helper.hpp"

#include <boost/algorithm/string.hpp>

using namespace std;

namespace corrector {

void ContigProcessor::ReadContig() {
    io::FileReadStream frs(contig_file_);
    io::SingleRead cur_read;
    frs >> cur_read;
    if (!frs.eof()) {
#pragma omp critical
        {
            ERROR("Non unique sequnce in one contig fasta!");
        }
    }
    contig_name_ = cur_read.name();
    contig_ = cur_read.GetSequenceString();

    output_contig_file_ = path::append_path(path::parent_path(contig_file_), path::basename(contig_file_) + ".ref.fasta");
    charts_.resize(contig_.length());
}

void ContigProcessor::UpdateOneRead(const SingleSamRead &tmp, MappedSamStream &sm) {
    unordered_map<size_t, position_description> all_positions;
    if (tmp.get_contig_id() < 0) {
        return;
    }
    auto cur_s = sm.get_contig_name(tmp.get_contig_id());
    if (contig_name_.compare(cur_s) != 0) {
        return;
    }
    tmp.CountPositions(all_positions, contig_.length());
    size_t error_num = 0;

    for (auto &pos : all_positions) {
        charts_[pos.first].update(pos.second);
        if (pos.second.FoundOptimal(contig_[pos.first]) != var_to_pos[(int) contig_[pos.first]]) {
            error_num++;
        }
    }

    if (error_num >= error_counts_.size())
        error_counts_[error_counts_.size() - 1]++;
    else
        error_counts_[error_num]++;
}

//returns: number of changed nucleotides;
size_t ContigProcessor::UpdateOneBase(size_t i, stringstream &ss, const unordered_map<size_t, position_description> &interesting_positions) const{
    char old = (char) toupper(contig_[i]);
    auto strat = corr_cfg::get().strat;
    size_t maxi = charts_[i].FoundOptimal(contig_[i]);
    auto i_position = interesting_positions.find(i);
    if (i_position != interesting_positions.end()) {
        size_t maxj = i_position->second.FoundOptimal(contig_[i]);
        if (maxj != maxi) {
            DEBUG("Interesting positions differ with majority!");
            DEBUG("On position " << i << "  old: " << old << " majority: " << pos_to_var[maxi] << "interesting: " << pos_to_var[maxj]);
            if (strat != Strategy::MajorityOnly)
                maxi = maxj;
        }
    }
    if (old != pos_to_var[maxi]) {
        DEBUG("On position " << i << " changing " << old << " to " << pos_to_var[maxi]);
        DEBUG(charts_[i].str());
        if (maxi < Variants::Deletion) {
            ss << pos_to_var[maxi];
            return 1;
        } else if (maxi == Variants::Deletion) {
            return 1;
        } else if (maxi == Variants::Insertion) {
            string maxj = "";
            //first base before insertion;
            size_t new_maxi = var_to_pos[(int) contig_[i]];
            int new_maxx = charts_[i].votes[new_maxi];
            for (size_t k = 0; k < MAX_VARIANTS; k++) {
                if (new_maxx < charts_[i].votes[k] && (k != Variants::Insertion) && (k != Variants::Deletion)) {
                    new_maxx = charts_[i].votes[k];
                    new_maxi = k;
                }
            }
            ss << pos_to_var[new_maxi];
            int max_ins = 0;
            for (const auto &ic : charts_[i].insertions) {
                if (ic.second > max_ins) {
                    max_ins = ic.second;
                    maxj = ic.first;
                }
            }
            DEBUG("most popular insertion: " << maxj);
            ss << maxj;
            if (old == maxj[0]) {
                return (int) maxj.length() - 1;
            } else {
                return (int) maxj.length();
            }
        } else {
            //something strange happened
            WARN("While processing base " << i << " unknown decision was made");
            return 0;
        }
    } else {
        ss << old;
        return 0;
    }
}

size_t ContigProcessor::ProcessMultipleSamFiles() {
    error_counts_.resize(kMaxErrorNum);
    for (const auto &sf : sam_files_) {
        MappedSamStream sm(sf.first);
        while (!sm.eof()) {
            SingleSamRead tmp;
            sm >> tmp;

            UpdateOneRead(tmp, sm);
        }
        sm.close();
    }

    ipp_.FillInterestingPositions(charts_);
    for (const auto &sf : sam_files_) {
        MappedSamStream sm(sf.first);
        while (!sm.eof()) {
            unordered_map<size_t, position_description> ps;
            if (sf.second == io::LibraryType::PairedEnd ) {
                PairedSamRead tmp;
                sm >> tmp;
                tmp.CountPositions(ps, contig_.length());
            } else {
                SingleSamRead tmp;
                sm >> tmp;
                tmp.CountPositions(ps, contig_.length());
            }
            ipp_.UpdateInterestingRead(ps);
        }
        sm.close();
    }
    ipp_.UpdateInterestingPositions();
    unordered_map<size_t, position_description> interesting_positions = ipp_.get_weights();
    stringstream s_new_contig;
    size_t total_changes = 0;
    for (size_t i = 0; i < contig_.length(); i++) {
        total_changes += UpdateOneBase(i, s_new_contig, interesting_positions);
    }
    vector<string> contig_name_splitted;
    boost::split(contig_name_splitted, contig_name_, boost::is_any_of("_"));
    if (contig_name_splitted.size() >= 8) {
        io::osequencestream_with_manual_node_id oss(output_contig_file_);
        oss.setNodeID(std::stoi(contig_name_splitted[1]));
        oss.setCoverage(std::stod(contig_name_splitted[5]));
        oss.setID(std::stoi(contig_name_splitted[7]));
        oss << s_new_contig.str();
    } else {
        io::osequencestream oss(output_contig_file_);
        oss << io::SingleRead(contig_name_, s_new_contig.str());
    }
    return total_changes;
}

}
;
